"""`search`: build samples.tsv for one biome, starting from its survey records."""

import csv
import re
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path

import yaml

from . import score
from .clients import ena, gold, mgnify, ncbi, ols
from .constants import (
    CO_ASSEMBLY_SAMPLE,
    CURATION_FLAGS,
    EXCLUDED_TAXA,
    EXCLUDING_FLAGS,
    GENOME_FIELDS,
    HARD_FILTERS,
    KEPT_COLUMNS,
    KINGDOM_TAXA,
    LAYER_TYPES,
    LONG_READ_PLATFORMS,
    MIN_RUN_READS,
    REPORTED_BIOME,
    SAMPLE_ID,
    SAMPLES_COLUMNS,
    SURVEILLANCE_PACKAGES,
)
from .log import counts, info, warn
from .models import (
    EnaGenome,
    EnaRunLibrary,
    EnaSampleMetadata,
    GoldBiosample,
    MgnifyStudy,
    NcbiAssembly,
    NcbiTaxdump,
    NcbiTaxon,
    SampleRow,
    SurveyRecord,
    from_row,
)
from .score import Evidence, Terms
from .survey import Config, genome_samples, taxa_filter, write_tsv


def survey_records(survey_dir: Path, name: str) -> list[SurveyRecord]:
    """The biome's survey records (survey/biomes/<name>/records.tsv), minus assemblies and runs
    whose sample has no evidence for the config's kingdom (`kingdom_match` false). Survey finds
    assemblies and runs by their metagenome taxon (`tax_eq`), so one without its own taxon takes
    that one, with its tax ID from survey.tsv."""
    path = survey_dir / "biomes" / name / "records.tsv"
    if not path.exists():
        raise ValueError(f"{path} not found: run `biome-finder survey` on this biome first")
    with open(survey_dir / "survey.tsv", newline="") as f:
        taxon_ids = {
            r["metagenome_taxon"]: r["tax_id"] for r in csv.DictReader(f, delimiter="\t") if r["biome"] == name
        }
    with open(path, newline="") as f:
        everything = [from_row(SurveyRecord, r) for r in csv.DictReader(f, delimiter="\t")]
    records = [r for r in everything if r.kingdom_match is not False]
    info(
        f"{name}: {len(everything):,} survey records in {path}; {len(everything) - len(records):,} assemblies "
        f"and runs left out (no evidence for the kingdom); kept {counts(Counter(r.data_type for r in records))}"
    )
    return [
        (
            replace(r, tax_id=taxon_ids.get(r.metagenome_taxon, ""), scientific_name=r.metagenome_taxon)
            if not r.tax_id and r.data_type in ("contigs", "reads")
            else r
        )
        for r in records
    ]


def record_flags(r: SurveyRecord) -> list[str]:
    """Flags that the survey record alone settles (PROPOSAL.md §3.6 `flags`)."""
    out = []
    if r.in_catalogue:
        out.append("in_existing_catalogue")
    if r.data_type == "reads":
        if r.instrument_platform in LONG_READ_PLATFORMS:
            out.append("long_read")
        if r.library_layout == "SINGLE":
            out.append("single_end")
        if (r.read_count or 0) < MIN_RUN_READS:
            out.append("low_reads")
    return out


def data_type_evidence(r: SurveyRecord) -> str:
    """The fields that decided a survey record's data type (PROPOSAL.md §3.4). Survey found
    MAGs, bins and assemblies by their ENA `assembly_type`, and runs by their library; an
    assembly found by its sample can be untyped, and counts as a community assembly when its
    organism is a metagenome."""
    if r.data_type == "reads":
        return "ena.library_source=METAGENOMIC; ena.library_strategy=WGS or OTHER"
    if r.assembly_type:
        return f"ena.assembly_type={r.assembly_type.lower()}"
    return f"ena.assembly_type=(empty); ena.scientific_name={r.scientific_name}"


def to_row(r: SurveyRecord) -> SampleRow:
    """A survey record as a samples.tsv row. A genome's biome is its source samples' biome, so
    `sample_accession` holds those (several for co-assemblies) when ENA links them."""
    reads = r.data_type == "reads"
    return SampleRow(
        record_accession=r.accession,
        data_type=r.data_type,
        data_type_evidence=data_type_evidence(r),
        source_db="ENA",
        sample_accession=[a for a in (r.source_samples or r.sample_accession).split(";") if a],
        study_accession=r.study_accession,
        reported_tax_id=r.tax_id,
        reported_tax_name=r.scientific_name,
        flags=record_flags(r),
        read_stats=f"{r.read_count or 0} reads; {r.library_layout}; {r.instrument_platform}" if reads else "",
        first_public=r.first_public,
        read_count=r.read_count or 0,
    )


def genome_query(cfg: Config) -> str:
    """ENA `assembly` query for the single-organism genomes whose sample names one of the config's
    `keywords_include` in an environment field, in the config's kingdom minus the excluded taxa,
    leaving out the typed MAG, bin and SAG layers (survey has the first). The Portal's wildcards
    can't span a space, so a phrase matches when each of its words is in the same field."""
    words = [k.split() for k in cfg.get("keywords_include") or []]
    terms = [
        "(" + " AND ".join(f'{field}="*{w}*"' for w in phrase) + ")" for phrase in words for field in GENOME_FIELDS
    ]
    trees = " OR ".join(f"tax_tree({t})" for t in KINGDOM_TAXA[cfg["kingdom"]])
    layers = " AND ".join(f'NOT assembly_type="{LAYER_TYPES[t]}"' for t in ("mag", "bin", "sag"))
    return f"({trees}){taxa_filter(cfg)} AND {layers} AND ({' OR '.join(terms)})"


def genome_type(g: EnaGenome, ncbi: NcbiAssembly | None = None) -> tuple[str, str]:
    """A genome's data type and the fields that decided it (PROPOSAL.md §3.4, A.4): NCBI's
    `excluded_from_refseq` (derived from metagenome → mag, from a single cell → sag), the
    BioSample package (MIMAG → mag, MISAG → sag), `environmental_sample` (metagenome-derived →
    mag), and a strain or culture collection for an isolate candidate. With none of those it's
    `unknown`, for a curator."""
    package = g.ncbi_reporting_standard
    evidence = [f"ena.assembly_type={g.assembly_type.lower()}"] if g.assembly_type else []
    evidence += [f"biosample.package={package}"] if package else []
    excluded = ncbi.excluded_from_refseq if ncbi else ""
    if "derived from metagenome" in excluded:
        return "mag", "; ".join([*evidence, f"ncbi.excluded_from_refseq={excluded}"])
    if "derived from single cell" in excluded:
        return "sag", "; ".join([*evidence, f"ncbi.excluded_from_refseq={excluded}"])
    if package.upper().startswith("MIMAG"):
        return "mag", "; ".join(evidence)
    if package.upper().startswith("MISAG"):
        return "sag", "; ".join(evidence)
    if g.environmental_sample == "true":
        return "mag", "; ".join([*evidence, "ena.environmental_sample=true"])
    culture = [
        f"biosample.{k}={v}" for k, v in (("strain", g.strain), ("culture_collection", g.culture_collection)) if v
    ]
    if culture:
        return "isolate", "; ".join(evidence + culture)
    return "unknown", "; ".join(evidence + ([f"biosample.isolate={g.isolate}"] if g.isolate else []))


def genome_rows(
    genomes: list[EnaGenome], known_samples: set[str], ncbi: dict[str, NcbiAssembly] | None = None
) -> list[SampleRow]:
    """samples.tsv rows for the genomes `search` finds itself, minus those whose own sample is a
    survey genome's (the same genome, as a GCA). A MAG takes its source samples (`derived from`,
    as primary IDs) like the survey's MAGs; an isolate is its own sample. `ncbi` (NCBI's assembly
    summary rows) types the NCBI genomes ENA leaves untyped. A type strain (NCBI's type material,
    or "type strain" in the strain) is flagged `type_strain`, since its sample often records
    where the strain was first isolated, not where this genome came from (A.4); a public-health
    surveillance isolate (One Health Enteric or Pathogen package, NCBI Pathogen Detection) is
    flagged `surveillance_isolate`: a faecal or food-borne pathogen sampled from water, for a
    curator to judge; a SAG is flagged `single_cell`."""
    ncbi = ncbi or {}
    rows: list[SampleRow] = []
    for g in genomes:
        if g.sample_accession in known_samples:
            continue
        listed = ncbi.get(g.accession)
        data_type, evidence = genome_type(g, listed)
        type_material = listed and "type material" in listed.relation_to_type_material
        surveillance = (
            g.ncbi_reporting_standard.startswith(SURVEILLANCE_PACKAGES) or "Pathogen Detection" in g.description
        )
        row = SampleRow(
            record_accession=g.accession,
            data_type=data_type,
            data_type_evidence=evidence,
            assembly_level=g.assembly_level,
            source_db="ENA",
            sample_accession=[g.sample_accession],
            study_accession=g.study_accession,
            reported_tax_id=g.tax_id,
            reported_tax_name=g.scientific_name,
            first_public=g.first_created[:10],
        )
        if type_material or re.search(r"\btype strain\b", g.strain, re.I):
            row.flag("type_strain")
        if surveillance and data_type in ("isolate", "unknown"):
            row.flag("surveillance_isolate")
        if data_type == "sag":
            row.flag("single_cell")
        rows.append(row)
    mags = [r for r in rows if r.data_type == "mag"]
    origin = [
        SurveyRecord("mag", r.record_accession, r.study_accession, sample_accession=r.sample_accession[0]) for r in mags
    ]
    genome_samples(origin)
    for r, o in zip(mags, origin):
        r.sample_accession = [a for a in o.source_samples.split(";") if a] or r.sample_accession
    return rows


def lineage(taxon: NcbiTaxon) -> str:
    """The taxon's ranks as `d__Bacteria;k__Pseudomonadati;p__…;s__…`, one prefix per NCBI rank,
    with an empty name where the lineage has no taxon at that rank (NCBI gives metagenomes the
    species rank only: `d__;…;s__freshwater metagenome`); empty for a taxon with no ranks."""
    if not any(taxon.ranks.values()):
        return ""
    return ";".join(f"{rank[0]}__{name or ''}" for rank, name in taxon.ranks.items())


def resolve_taxonomy(rows: list[SampleRow], cfg: Config, dump: NcbiTaxdump) -> None:
    """Fill `tax_lineage`, `domain` and `taxdump_date` from the release `dump`, following a merged
    tax ID to the current one. ENA keeps the reported tax ID and name current, so the resolved
    ones aren't stored: in freshwater-water-prok (2026-11) no ID resolves to another, and names
    differ only by NCBI renames of the same taxon (`Binatia bacterium (SeqCode)`). Flags
    `out_of_scope_taxon` on genomes whose resolved lineage is outside the config's kingdom or
    under an excluded taxon (defaults plus `exclude_taxa_extra`). A tax ID that doesn't resolve
    is flagged `taxid_deleted` (in `delnodes.dmp`) or `taxid_unknown` (e.g. newer than the
    release), not out of scope: ENA can still use a tax ID NCBI deleted, so a curator decides.
    Assemblies and runs carry their metagenome taxon, so the scope doesn't apply to them."""
    excluded = {str(t) for t in EXCLUDED_TAXA + (cfg.get("exclude_taxa_extra") or [])}
    kingdoms = {k: {str(t) for t in ids} for k, ids in KINGDOM_TAXA.items()}
    for r in rows:
        reported = r.reported_tax_id
        current = dump.merged.get(reported, reported)
        taxon = dump.taxa.get(current)
        tree = taxon.ancestors | {current} if taxon else frozenset()
        r.domain = next((k for k, ids in kingdoms.items() if ids & tree), "")
        r.tax_lineage = lineage(taxon) if taxon else ""
        r.taxdump_date = dump.release
        if not taxon and reported:
            r.flag("taxid_deleted" if reported in dump.deleted else "taxid_unknown")
        genome = r.data_type not in ("contigs", "reads")
        if genome and taxon and (tree & excluded or r.domain != cfg["kingdom"]):
            r.flag("out_of_scope_taxon")


def reported_biome(samples: list[EnaSampleMetadata]) -> str:
    """The samples' biome fields as found, `key=value; …` (REPORTED_BIOME order, empty ones left
    out); a co-assembly's samples that differ are `|`-joined."""
    parts = []
    for key, field in REPORTED_BIOME:
        values = dict.fromkeys(v for s in samples if (v := getattr(s, field)))
        if values:
            parts.append(f"{key}={' | '.join(values)}")
    return "; ".join(parts)


def run_flags(rows: list[SampleRow], libraries: list[EnaRunLibrary]) -> None:
    """Flags from each run's library and its sample's runs: `amplicon_only` when the library is
    amplicon (strategy AMPLICON, or selection PCR: amplicon data the survey's WGS/OTHER query
    lets through), and `sample_low_reads` when none of the sample's runs reaches MIN_RUN_READS,
    so its reads can't be routed (CLAUDE.md). Both send the row to tier D. The library also goes
    to `data_type_evidence`."""
    lib = {x.run_accession: x for x in libraries}
    best: dict[str, int] = defaultdict(int)
    reads = [r for r in rows if r.data_type == "reads"]
    for r in reads:
        for s in r.sample_accession:
            best[s] = max(best[s], r.read_count)
    for r in reads:
        x = lib.get(r.record_accession)
        if x:
            r.data_type_evidence = (
                f"ena.library_source=METAGENOMIC; ena.library_strategy={x.library_strategy}; "
                f"ena.library_selection={x.library_selection}"
            )
            if x.library_strategy == "AMPLICON" or x.library_selection == "PCR":
                r.flag("amplicon_only")
        if max((best[s] for s in r.sample_accession), default=0) < MIN_RUN_READS:
            r.flag("sample_low_reads")


def expand_co_assemblies(rows: list[SampleRow], metadata: list[EnaSampleMetadata]) -> list[EnaSampleMetadata]:
    """Replace ENA's virtual co-assembly samples ("mock sample for coassembly … of samples:
    ERS…,ERS…") in `sample_accession` with the samples they list, as primary IDs, since a
    co-assembled genome's biome is its source samples'. Returns those samples' metadata."""
    members = {
        s.sample_accession: SAMPLE_ID.findall(f"{s.sample_title} {s.sample_description}")
        for s in metadata
        if CO_ASSEMBLY_SAMPLE.search(s.sample_title)
    }
    if not members:
        return []
    found = ena.sample_metadata({m for ms in members.values() for m in ms}).data
    primary = {i: s.sample_accession for s in found for i in s.ids}
    expanded = 0
    for r in rows:
        if any(s in members for s in r.sample_accession):
            r.sample_accession = list(
                dict.fromkeys(primary.get(m, m) for s in r.sample_accession for m in members.get(s, [s]))
            )
            expanded += 1
    info(
        f"{len(members):,} virtual co-assembly samples replaced by the {len(found):,} samples they list, "
        f"in {expanded:,} rows"
    )
    return found


def sample_evidence(rows: list[SampleRow], metadata: list[EnaSampleMetadata]) -> None:
    """Fill `reported_biome` from each row's samples (a genome's source samples), add the hard
    filters their text matches to `flags` (HARD_FILTERS), and put rows carrying an
    excluding flag in tier D. A genome is flagged if any of its source samples is. An isolate is
    a culture by definition, so `culture` and `enrichment` don't apply to it."""
    by_id = {i: s for s in metadata for i in s.ids}
    for r in rows:
        samples = [by_id[a] for a in r.sample_accession if a in by_id]
        r.reported_biome = reported_biome(samples)
        hard = list(dict.fromkeys(f for s in samples for f in score.hard_flags(s)))
        if r.data_type == "isolate":
            hard = [f for f in hard if f not in ("culture", "enrichment")]
        r.flag(*hard)
        if EXCLUDING_FLAGS & set(r.flags):
            r.tier = "D"


def envo_terms(ids: list[str]) -> Terms:
    """The config's ENVO terms with everything under them (OLS)."""
    return Terms.of(t for i in ids for t in ols.with_descendants(i).data)


def siblings(rows: list[SampleRow], items: dict[str, list[Evidence]]) -> dict[str, float]:
    """Action 3, per study: (samples whose own metadata points in − samples pointing out) / all
    its samples, from −1 to 1. A sample points in with an `in` metadata item and no opposing one."""
    samples_of: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        samples_of[r.study_accession].update(r.sample_accession)
    lean = {}
    for s, its in items.items():
        meta = {e.claim for e in its if e.action == "sample_metadata"}
        lean[s] = -1 if meta & {"out", "excluded_envo"} else 1 if "in" in meta else 0
    return {
        study: sum(lean.get(s, 0) for s in samples) / len(samples) for study, samples in samples_of.items() if samples
    }


def biome_evidence(
    rows: list[SampleRow],
    cfg: Config,
    metadata: list[EnaSampleMetadata],
    study_texts: dict[str, str],
    include: Terms,
    exclude: Terms,
    mgnify_studies: dict[str, MgnifyStudy] | None = None,
    gold_samples: dict[str, GoldBiosample] | None = None,
) -> None:
    """Score each row (actions 1-5, score.py) and fill `evidence`, `score`, `tier`,
    `needs_curation` and `curated_biome`. The sample is the unit: each of a row's samples is
    scored on its own metadata plus the row's study, and the row takes its worst sample, since a
    co-assembled genome is accepted only if all its source samples are. A sample with an env_*
    term under `envo_exclude` is flagged `excluded_envo` and goes to tier D; one whose own text
    has a `keywords_exclude` word is flagged `off_biome_keyword` (a conflict, so tier C unless
    another flag excludes it), so the reason shows next to the hard filters; a sample whose own
    evidence points only to another biome (score.off_biome) is flagged `off_biome`, tier D; tier D
    from the hard filters stays. `mgnify_studies` (by ENA study ID) gives the study's curated
    MGnify biome (action 1), also written to `third_party_biome`; `gold_samples` (by BioSample)
    each sample's GOLD ecosystem, likewise."""
    mgnify_studies, gold_samples = mgnify_studies or {}, gold_samples or {}
    kw_in, kw_out = score.keyword(cfg.get("keywords_include") or []), score.keyword(cfg.get("keywords_exclude") or [])
    taxa = {t.lower() for t in cfg.get("metagenome_taxa") or []}
    by_id = {i: s for s in metadata for i in s.ids}
    items = {}
    for a in {a for r in rows for a in r.sample_accession}:
        its = score.sample_items(by_id[a], include, exclude, kw_in, kw_out, taxa) if a in by_id else []
        if a in gold_samples:
            its = score.gold_items(gold_samples[a].gold_id, gold_samples[a].ecosystem_path, cfg) + its
        if its:
            items[a] = its
    lean = siblings(rows, items)
    study_items = {s: score.study_items(s, t, kw_in, kw_out, taxa) for s, t in study_texts.items()}
    for s in {r.study_accession for r in rows} & set(mgnify_studies):
        mg = mgnify_studies[s]
        study_items[s] = score.third_party_items(f"mgnify:{mg.accession}", mg.lineage, cfg) + study_items.get(s, [])
    for r in rows:
        study = study_items.get(r.study_accession, [])
        samples = [items.get(a, []) for a in r.sample_accession] or [[]]
        scored = [score.combine(its + study, lean.get(r.study_accession, 0.0)) for its in samples]
        value, conflict = min(s for s, _ in scored), any(c for _, c in scored)
        found = list(dict.fromkeys(e for its in samples for e in its + study))
        if any(e.claim == "excluded_envo" for e in found):
            r.flag("excluded_envo")
        if any(e.action == "sample_metadata" and e.claim == "out" for e in found):
            r.flag("off_biome_keyword")
        if any(score.off_biome(its) for its in samples if its):
            r.flag("off_biome")
        r.evidence = [str(e) for e in found]
        listed = mgnify_studies.get(r.study_accession)
        third = [f"mgnify:{listed.accession}:{listed.lineage}"] if listed and listed.lineage else []
        third += [f"gold:{g.gold_id}:{g.ecosystem_path}" for a in r.sample_accession if (g := gold_samples.get(a))]
        r.third_party_biome = "; ".join(dict.fromkeys(third))
        r.score = int(value)
        if EXCLUDING_FLAGS & set(r.flags):
            r.tier = "D"
        else:
            curate = CURATION_FLAGS & set(r.flags) or r.data_type == "unknown"
            r.tier = "C" if curate else score.tier(value, conflict)
        r.needs_curation = r.tier == "C"  # conflicts land there; tier D is excluded, not curated
        r.curated_biome = cfg.get("catalogue_biome", "") if r.tier in ("A", "B") else ""


def study_flags(rows: list[SampleRow], study_texts: dict[str, str]) -> None:
    """Flags from the study's title and description, for samples that carry only a default title:
    `virome` when the study sequenced the viral fraction (e.g. the Yangtze River virome,
    PRJNA672207; in freshwater-water-prok every study that names it is a virome or a cell-free
    fraction, and a study mixing viromes with cellular metagenomes would lose those too), and
    `manipulated` when it ran a micro- or mesocosm experiment (e.g. PRJEB56877 at Jiřická Pond),
    which only sends the rows to curation."""
    for flag in ("virome", "manipulated"):
        hits = {s for s, t in study_texts.items() if HARD_FILTERS[flag].search(t)}
        for r in rows:
            if r.study_accession in hits:
                r.flag(flag)


def merge(old: list[SampleRow], new: list[SampleRow]) -> list[SampleRow]:
    """New rows with the human and `prepare` columns of the old row of the same record. Old rows
    the new search no longer finds are kept unchanged, so no curation is lost."""
    by_acc = {r.record_accession: r for r in old}
    merged = []
    for r in new:
        prev = by_acc.pop(r.record_accession, None)
        kept = {k: getattr(prev, k) for k in KEPT_COLUMNS if prev and getattr(prev, k)}
        merged.append(replace(r, **kept))
    return merged + list(by_acc.values())


def read_samples(path: Path) -> list[SampleRow]:
    """An existing samples.tsv, as rows."""
    with open(path, newline="") as f:
        return [SampleRow.from_tsv(r) for r in csv.DictReader(f, delimiter="\t")]


def run(
    biome_path: str, survey_dir: str | Path, out: str | Path, since: str = "", taxdump: str = "", gold_export: str = ""
) -> None:
    """Write `out` (samples.tsv) for one biome config, merged into `out` if it exists. `taxdump`
    (YYYY-MM-DD) picks the monthly NCBI taxonomy release of that month; the default is the
    latest. `gold_export` is the GOLD data export the curator downloaded (docs/search.md)."""
    with open(biome_path) as f:
        cfg = yaml.safe_load(f)
    name = cfg["name"]
    info(f"{name}: search ({biome_path})")
    records = survey_records(Path(survey_dir), name)
    rows = [to_row(r) for r in records]
    info(f"{name}: genomes whose sample names the biome (isolates, and MAGs outside ENA's MAG layer)")
    genomes = ena.genomes(genome_query(cfg)).data
    known = {r.sample_accession for r in records if r.data_type in ("mag", "bin")}
    info(f"{name}: NCBI assembly summary for {len(genomes):,} genomes (one ~2 GB copy a day, cached)")
    summary = ncbi.assembly_summary(g.accession for g in genomes).data
    extra = genome_rows(genomes, known, summary)
    rows += extra
    info(
        f"{name}: {len(genomes):,} genomes, {len(genomes) - len(extra):,} already in the survey; added "
        f"{counts(Counter(r.data_type for r in extra))}"
    )
    if since:
        rows = [r for r in rows if r.first_public >= since]
        info(f"{name}: {len(rows):,} records public on or after {since}")

    release = f"{taxdump[:7]}-01" if taxdump else ncbi.latest_taxdump().data
    info(f"{name}: taxonomy from NCBI taxdump {release}")
    dump = ncbi.taxdump(release, {r.reported_tax_id for r in rows}).data
    if not dump:
        raise ValueError(f"NCBI taxdump release {release!r} not found in the archive")
    resolve_taxonomy(rows, cfg, dump)
    info(f"{name}: taxonomy flags: {flag_counts(rows, ('out_of_scope_taxon', 'taxid_deleted', 'taxid_unknown'))}")

    runs = [r.record_accession for r in rows if r.data_type == "reads"]
    info(f"{name}: libraries of {len(runs):,} runs from ENA")
    run_flags(rows, ena.run_libraries(runs).data)
    info(f"{name}: run flags: {flag_counts(rows, ('amplicon_only', 'sample_low_reads'))}")

    samples = {a for r in rows for a in r.sample_accession}
    info(f"{name}: sample metadata for {len(samples):,} samples from ENA")
    metadata = ena.sample_metadata(samples).data
    found = {i for s in metadata for i in s.ids}
    info(f"{name}: {len(samples & found):,} samples found, {len(samples - found):,} not returned by ENA")
    metadata += expand_co_assemblies(rows, metadata)
    sample_evidence(rows, metadata)
    info(
        f"{name}: hard filters (action 0; enrichment goes to curation, the rest to tier D): "
        f"{flag_counts(rows, tuple(HARD_FILTERS))}"
    )

    studies = {r.study_accession for r in rows}
    info(f"{name}: titles and descriptions of {len(studies):,} studies from ENA")
    texts = {s.study_accession: s.text for s in ena.studies(studies).data}
    study_flags(rows, texts)
    info(f"{name}: flags from study texts: {flag_counts(rows, ('virome', 'manipulated'))}")
    include, exclude = envo_terms(cfg.get("envo_include") or []), envo_terms(cfg.get("envo_exclude") or [])
    info(
        f"{name}: ENVO from OLS: {len(cfg.get('envo_include') or [])} include terms → {len(include.ids):,} with "
        f"descendants, {len(cfg.get('envo_exclude') or [])} exclude terms → {len(exclude.ids):,}"
    )
    info(f"{name}: MGnify study biomes (action 1)")
    by_ena = {a: m for m in mgnify.studies().data for a in m.ena_accessions}
    in_mgnify = studies & set(by_ena)
    info(
        f"{name}: {len(in_mgnify):,} of {len(studies):,} studies have an MGnify biome: "
        f"{counts(Counter(by_ena[s].lineage for s in in_mgnify))}"
    )
    gold_samples: dict[str, GoldBiosample] = {}
    if gold_export:
        listed = {g.biosample_accession: g for g in gold.biosamples(gold_export).data}
        gold_samples = {a: listed[a] for a in samples | {m for r in rows for m in r.sample_accession} if a in listed}
        info(
            f"{name}: {len(gold_samples):,} samples have a GOLD ecosystem: "
            f"{counts(Counter(' > '.join(g.ecosystem_path.split(' > ')[:4]) for g in gold_samples.values()))}"[:600]
        )
    elif cfg.get("gold_ecosystems"):
        warn(f"{name}: the config has gold_ecosystems but no --gold export was given, so no GOLD evidence")
    biome_evidence(rows, cfg, metadata, texts, include, exclude, by_ena, gold_samples)
    info(f"{name}: biome flags: {flag_counts(rows, ('excluded_envo', 'off_biome', 'off_biome_keyword'))}")
    for data_type, tiers in sorted(tier_counts(rows).items()):
        info(f"{name}: {data_type}: " + ", ".join(f"tier {t} {tiers.get(t, 0):,}" for t in "ABCD"))
    info(
        f"{name}: {sum(r.tier == 'C' for r in rows):,} of {len(rows):,} rows in tier C, "
        f"{sum(bool(r.needs_curation) for r in rows):,} need curation"
    )

    out = Path(out)
    old: list[SampleRow] = []
    if out.exists():
        old = read_samples(out)
        decided = sum(bool(r.decision) for r in old)
        info(f"{name}: merging into {out}: {len(old):,} existing rows, {decided:,} with a curator decision kept")
    merged = merge(old, rows)
    write_tsv(out, SAMPLES_COLUMNS, merged)
    info(f"{name}: wrote {len(merged):,} rows to {out}")


def flag_counts(rows: list[SampleRow], flags: tuple[str, ...]) -> str:
    """How many rows carry each of `flags`."""
    return counts(dict(Counter(f for r in rows for f in r.flags if f in flags)))


def tier_counts(rows: list[SampleRow]) -> dict[str, dict[str, int]]:
    """Rows per tier, per data type."""
    out: dict[str, Counter[str]] = defaultdict(Counter)
    for r in rows:
        out[r.data_type][r.tier] += 1
    return {k: dict(v) for k, v in out.items()}
