"""`search`: build samples.tsv for one biome, starting from its survey records."""

import csv
from pathlib import Path
from typing import Any

import yaml

from . import score
from collections import Counter, defaultdict

from .clients import ena, mgnify, ncbi, ols
from .log import counts, info
from .constants import EXCLUDED_TAXA, KINGDOM_TAXA, LONG_READ_PLATFORMS, SAMPLE_ID
from .models import EnaSampleMetadata, MgnifyStudy, NcbiTaxdump, NcbiTaxon
from .score import Evidence, Terms
from .survey import Config, write_tsv

Row = dict[str, Any]  # one samples.tsv row (SAMPLES_COLUMNS)

SAMPLES_COLUMNS = [  # PROPOSAL.md §3.6; columns are added, never renamed
    "record_accession",
    "data_type",
    "data_type_evidence",
    "assembly_level",
    "source_db",
    "sample_accession",
    "related_accessions",
    "study_accession",
    "secondary_study_accession",
    "bioproject",
    "reported_tax_id",
    "reported_tax_name",
    "tax_lineage",
    "domain",
    "taxdump_date",
    "publication",
    "paper_access",
    "paper_pdf",
    "reported_biome",
    "third_party_biome",
    "curated_biome",
    "score",
    "tier",
    "needs_curation",
    "evidence",
    "flags",
    "read_stats",
    "first_public",
    "curation_note",
    "decision",
    "destination",
    "local_path",
    "prepare_status",
]
KEPT_COLUMNS = ("curation_note", "decision", "destination", "local_path", "prepare_status")
MIN_RUN_READS = 200_000  # a sample's reads are routed only if one of its runs reaches this
REPORTED_BIOME = [  # (reported_biome key, ENA sample field)
    ("env_broad_scale", "broad_scale_environmental_context"),
    ("env_local_scale", "local_environmental_context"),
    ("env_medium", "environmental_medium"),
    ("isolation_source", "isolation_source"),
    ("organism", "scientific_name"),
    ("host", "host"),
    ("body_site", "host_body_site"),
]


def survey_records(survey_dir: Path, name: str) -> list[dict[str, str]]:
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
        everything = list(csv.DictReader(f, delimiter="\t"))
    records = [r for r in everything if r["kingdom_match"] != "False"]
    info(
        f"{name}: {len(everything):,} survey records in {path}; {len(everything) - len(records):,} assemblies "
        f"and runs left out (no evidence for the kingdom); kept {counts(Counter(r['data_type'] for r in records))}"
    )
    for r in records:
        if not r["tax_id"] and r["data_type"] in ("contigs", "reads"):
            r["tax_id"], r["scientific_name"] = taxon_ids.get(r["metagenome_taxon"], ""), r["metagenome_taxon"]
    return records


def flags(r: dict[str, str]) -> list[str]:
    """Flags that the survey record alone settles (PROPOSAL.md §3.6 `flags`)."""
    out = []
    if r["in_catalogue"] == "True":
        out.append("in_existing_catalogue")
    if r["data_type"] == "reads":
        if r["instrument_platform"] in LONG_READ_PLATFORMS:
            out.append("long_read")
        if r["library_layout"] == "SINGLE":
            out.append("single_end")
        if int(r["read_count"] or 0) < MIN_RUN_READS:
            out.append("low_reads")
    return out


def data_type_evidence(r: dict[str, str]) -> str:
    """The fields that decided a survey record's data type (PROPOSAL.md §3.4). Survey found
    MAGs, bins and assemblies by their ENA `assembly_type`, and runs by their library; an
    assembly found by its sample can be untyped, and counts as a community assembly when its
    organism is a metagenome."""
    if r["data_type"] == "reads":
        return "ena.library_source=METAGENOMIC; ena.library_strategy=WGS or OTHER"
    if r.get("assembly_type"):
        return f"ena.assembly_type={r['assembly_type'].lower()}"
    return f"ena.assembly_type=(empty); ena.scientific_name={r['scientific_name']}"


def to_row(r: dict[str, str]) -> Row:
    """A survey record as a samples.tsv row. A genome's biome is its source samples' biome, so
    `sample_accession` holds those (`;`-joined for co-assemblies) when ENA links them."""
    reads = r["data_type"] == "reads"
    return {
        "record_accession": r["accession"],
        "data_type": r["data_type"],
        "data_type_evidence": data_type_evidence(r),
        "source_db": "ENA",
        "sample_accession": r["source_samples"] or r["sample_accession"],
        "study_accession": r["study_accession"],
        "reported_tax_id": r["tax_id"],
        "reported_tax_name": r["scientific_name"],
        "flags": ";".join(flags(r)),
        "read_stats": f"{r['read_count']} reads; {r['library_layout']}; {r['instrument_platform']}" if reads else "",
        "first_public": r["first_public"],
    }


def lineage(taxon: NcbiTaxon) -> str:
    """The taxon's ranks as `d__Bacteria;k__Pseudomonadati;p__…;s__…`, one prefix per NCBI rank,
    with an empty name where the lineage has no taxon at that rank (NCBI gives metagenomes the
    species rank only: `d__;…;s__freshwater metagenome`); empty for a taxon with no ranks."""
    if not any(taxon.ranks.values()):
        return ""
    return ";".join(f"{rank[0]}__{name or ''}" for rank, name in taxon.ranks.items())


def resolve_taxonomy(rows: list[Row], cfg: Config, dump: NcbiTaxdump) -> None:
    """Fill `tax_lineage`, `domain` and `taxdump_date` from the release `dump`, following a merged
    tax ID to the current one. ENA keeps the reported tax ID and name current, so the resolved
    ones aren't stored: in freshwater-water-prok (2026-11) no ID resolves to another, and names
    differ only by NCBI renames of the same taxon (`Binatia bacterium (SeqCode)`). Flags `out_of_scope_taxon` on genomes whose resolved
    lineage is outside the config's kingdom or under an excluded taxon (defaults plus
    `exclude_taxa_extra`). A tax ID that doesn't resolve is flagged `taxid_deleted` (in
    `delnodes.dmp`) or `taxid_unknown` (e.g. newer than the release), not out of scope: ENA can
    still use a tax ID NCBI deleted, so a curator decides. Assemblies and runs carry their
    metagenome taxon, so the scope doesn't apply to them."""
    excluded = {str(t) for t in EXCLUDED_TAXA + (cfg.get("exclude_taxa_extra") or [])}
    kingdoms = {k: {str(t) for t in ids} for k, ids in KINGDOM_TAXA.items()}
    for r in rows:
        reported: str = r["reported_tax_id"]
        current = dump.merged.get(reported, reported)
        taxon = dump.taxa.get(current)
        tree = taxon.ancestors | {current} if taxon else frozenset()
        domain = next((k for k, ids in kingdoms.items() if ids & tree), "")
        new_flags = []
        if not taxon and reported:
            new_flags.append("taxid_deleted" if reported in dump.deleted else "taxid_unknown")
        genome = r["data_type"] not in ("contigs", "reads")
        if genome and taxon and (tree & excluded or domain != cfg["kingdom"]):
            new_flags.append("out_of_scope_taxon")
        r.update(
            tax_lineage=lineage(taxon) if taxon else "",
            domain=domain,
            taxdump_date=dump.release,
            flags=";".join(filter(None, [r["flags"], *new_flags])),
        )


def reported_biome(samples: list[EnaSampleMetadata]) -> str:
    """The samples' biome fields as found, `key=value; …` (REPORTED_BIOME order, empty ones left
    out); a co-assembly's samples that differ are `|`-joined."""
    parts = []
    for key, field in REPORTED_BIOME:
        values = dict.fromkeys(v for s in samples if (v := getattr(s, field)))
        if values:
            parts.append(f"{key}={' | '.join(values)}")
    return "; ".join(parts)


def expand_co_assemblies(rows: list[Row], metadata: list[EnaSampleMetadata]) -> list[EnaSampleMetadata]:
    """Replace ENA's virtual co-assembly samples ("mock sample for coassembly … of samples:
    ERS…,ERS…") in `sample_accession` with the samples they list, as primary IDs, since a
    co-assembled genome's biome is its source samples'. Returns those samples' metadata."""
    members = {
        s.sample_accession: SAMPLE_ID.findall(f"{s.sample_title} {s.sample_description}")
        for s in metadata
        if score.CO_ASSEMBLY_SAMPLE.search(s.sample_title)
    }
    if not members:
        return []
    found = ena.sample_metadata({m for ms in members.values() for m in ms}).data
    primary = {i: s.sample_accession for s in found for i in s.ids}
    expanded = 0
    for r in rows:
        samples = r["sample_accession"].split(";")
        if any(s in members for s in samples):
            r["sample_accession"] = ";".join(
                dict.fromkeys(primary.get(m, m) for s in samples for m in members.get(s, [s]))
            )
            expanded += 1
    info(
        f"{len(members):,} virtual co-assembly samples replaced by the {len(found):,} samples they list, "
        f"in {expanded:,} rows"
    )
    return found


def sample_evidence(rows: list[Row], metadata: list[EnaSampleMetadata]) -> None:
    """Fill `reported_biome` from each row's samples (a genome's source samples), add the hard
    filters their text matches to `flags` (score.HARD_FILTERS), and put rows carrying an
    excluding flag in tier D. A genome is flagged if any of its source samples is."""
    by_id = {i: s for s in metadata for i in s.ids}
    for r in rows:
        samples = [by_id[a] for a in r["sample_accession"].split(";") if a in by_id]
        r["reported_biome"] = reported_biome(samples)
        hard = dict.fromkeys(f for s in samples for f in score.hard_flags(s))
        flags = [f for f in r["flags"].split(";") if f] + list(hard)
        r["flags"] = ";".join(flags)
        if score.EXCLUDING_FLAGS & set(flags):
            r["tier"] = "D"


def envo_terms(ids: list[str]) -> Terms:
    """The config's ENVO terms with everything under them (OLS)."""
    return Terms.of(t for i in ids for t in ols.with_descendants(i).data)


def siblings(rows: list[Row], items: dict[str, list[Evidence]]) -> dict[str, float]:
    """Action 3, per study: (samples whose own metadata points in − samples pointing out) / all
    its samples, from −1 to 1. A sample points in with an `in` metadata item and no opposing one."""
    samples_of: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        samples_of[r["study_accession"]].update(a for a in r["sample_accession"].split(";") if a)
    lean = {}
    for s, its in items.items():
        meta = {e.claim for e in its if e.action == "sample_metadata"}
        lean[s] = -1 if meta & {"out", "excluded_envo"} else 1 if "in" in meta else 0
    return {study: sum(lean.get(s, 0) for s in samples) / len(samples) for study, samples in samples_of.items()}


def biome_evidence(
    rows: list[Row],
    cfg: Config,
    metadata: list[EnaSampleMetadata],
    study_texts: dict[str, str],
    include: Terms,
    exclude: Terms,
    mgnify_studies: dict[str, MgnifyStudy] | None = None,
) -> None:
    """Score each row (actions 1-5, score.py) and fill `evidence`, `score`, `tier`,
    `needs_curation` and `curated_biome`. The sample is the unit: each of a row's samples is
    scored on its own metadata plus the row's study, and the row takes its worst sample, since a
    co-assembled genome is accepted only if all its source samples are. A sample with an env_*
    term under `envo_exclude` is flagged `excluded_envo` and goes to tier D; tier D from the hard
    filters stays. `mgnify_studies` (by ENA study ID) gives the study's curated MGnify biome
    (action 1), also written to `third_party_biome`."""
    mgnify_studies = mgnify_studies or {}
    kw_in, kw_out = score.keyword(cfg.get("keywords_include") or []), score.keyword(cfg.get("keywords_exclude") or [])
    taxa = {t.lower() for t in cfg.get("metagenome_taxa") or []}
    by_id = {i: s for s in metadata for i in s.ids}
    items = {
        a: score.sample_items(by_id[a], include, exclude, kw_in, kw_out, taxa)
        for r in rows
        for a in r["sample_accession"].split(";")
        if a in by_id
    }
    lean = siblings(rows, items)
    study_items = {s: score.study_items(s, t, kw_in, kw_out, taxa) for s, t in study_texts.items()}
    for s in {r["study_accession"] for r in rows} & set(mgnify_studies):
        mg = mgnify_studies[s]
        study_items[s] = score.third_party_items(f"mgnify:{mg.accession}", mg.lineage, cfg) + study_items.get(s, [])
    for r in rows:
        study = study_items.get(r["study_accession"], [])
        samples = [items.get(a, []) for a in r["sample_accession"].split(";") if a] or [[]]
        scored = [score.combine(its + study, lean.get(r["study_accession"], 0.0)) for its in samples]
        value, conflict = min(s for s, _ in scored), any(c for _, c in scored)
        found = list(dict.fromkeys(e for its in samples for e in its + study))
        flags = [f for f in r["flags"].split(";") if f]
        if any(e.claim == "excluded_envo" for e in found) and "excluded_envo" not in flags:
            flags.append("excluded_envo")
        r["flags"] = ";".join(flags)
        r["evidence"] = "; ".join(str(e) for e in found)
        listed = mgnify_studies.get(r["study_accession"])
        r["third_party_biome"] = f"mgnify:{listed.accession}:{listed.lineage}" if listed and listed.lineage else ""
        r["score"] = value
        r["tier"] = "D" if score.EXCLUDING_FLAGS & set(flags) else score.tier(value, conflict)
        r["needs_curation"] = r["tier"] == "C"  # conflicts land there; tier D is excluded, not curated
        r["curated_biome"] = cfg.get("catalogue_biome", "") if r["tier"] in ("A", "B") else ""


def merge(old: list[Row], new: list[Row]) -> list[Row]:
    """New rows with the human and `prepare` columns of the old row of the same record. Old rows
    the new search no longer finds are kept unchanged, so no curation is lost."""
    by_acc = {r["record_accession"]: r for r in old}
    merged = []
    for r in new:
        prev = by_acc.pop(r["record_accession"], {})
        merged.append({**r, **{k: prev[k] for k in KEPT_COLUMNS if prev.get(k)}})
    return merged + list(by_acc.values())


def run(biome_path: str, survey_dir: str | Path, out: str | Path, since: str = "", taxdump: str = "") -> None:
    """Write `out` (samples.tsv) for one biome config, merged into `out` if it exists. `taxdump`
    (YYYY-MM-DD) picks the monthly NCBI taxonomy release of that month; the default is the
    latest."""
    with open(biome_path) as f:
        cfg = yaml.safe_load(f)
    name = cfg["name"]
    info(f"{name}: search ({biome_path})")
    rows = [to_row(r) for r in survey_records(Path(survey_dir), name)]
    if since:
        rows = [r for r in rows if r["first_public"] >= since]
        info(f"{name}: {len(rows):,} records public on or after {since}")

    release = f"{taxdump[:7]}-01" if taxdump else ncbi.latest_taxdump().data
    info(f"{name}: taxonomy from NCBI taxdump {release}")
    dump = ncbi.taxdump(release, {r["reported_tax_id"] for r in rows}).data
    if not dump:
        raise ValueError(f"NCBI taxdump release {release!r} not found in the archive")
    resolve_taxonomy(rows, cfg, dump)
    info(f"{name}: taxonomy flags: {flag_counts(rows, ('out_of_scope_taxon', 'taxid_deleted', 'taxid_unknown'))}")

    samples = {a for r in rows for a in r["sample_accession"].split(";") if a}
    info(f"{name}: sample metadata for {len(samples):,} samples from ENA")
    metadata = ena.sample_metadata(samples).data
    found = {i for s in metadata for i in s.ids}
    info(f"{name}: {len(samples & found):,} samples found, {len(samples - found):,} not returned by ENA")
    metadata += expand_co_assemblies(rows, metadata)
    sample_evidence(rows, metadata)
    info(f"{name}: hard filters (action 0): {flag_counts(rows, tuple(score.HARD_FILTERS))}")

    studies = {r["study_accession"] for r in rows}
    info(f"{name}: titles and descriptions of {len(studies):,} studies from ENA")
    texts = {s.study_accession: s.text for s in ena.studies(studies).data}
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
    biome_evidence(rows, cfg, metadata, texts, include, exclude, by_ena)
    for data_type, tiers in sorted(tier_counts(rows).items()):
        info(f"{name}: {data_type}: " + ", ".join(f"tier {t} {tiers.get(t, 0):,}" for t in "ABCD"))
    info(
        f"{name}: {sum(r['tier'] == 'C' for r in rows):,} of {len(rows):,} rows in tier C, "
        f"{sum(bool(r['needs_curation']) for r in rows):,} need curation"
    )

    out = Path(out)
    old: list[Row] = []
    if out.exists():
        with open(out, newline="") as f:
            old = list(csv.DictReader(f, delimiter="\t"))
        decided = sum(bool(r.get("decision")) for r in old)
        info(f"{name}: merging into {out}: {len(old):,} existing rows, {decided:,} with a curator decision kept")
    merged = merge(old, rows)
    write_tsv(out, SAMPLES_COLUMNS, merged)
    info(f"{name}: wrote {len(merged):,} rows to {out}")


def flag_counts(rows: list[Row], flags: tuple[str, ...]) -> str:
    """How many rows carry each of `flags`."""
    found = Counter(f for r in rows for f in r["flags"].split(";") if f in flags)
    return counts(dict(found))


def tier_counts(rows: list[Row]) -> dict[str, dict[str, int]]:
    """Rows per tier, per data type."""
    out: dict[str, Counter[str]] = defaultdict(Counter)
    for r in rows:
        out[r["data_type"]][r["tier"]] += 1
    return {k: dict(v) for k, v in out.items()}
