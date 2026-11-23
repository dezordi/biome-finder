"""`survey`: count public INSDC data per candidate biome and write a survey folder."""

import csv
import re
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import asdict
from pathlib import Path
from typing import Any

import yaml

from . import report
from .clients import ena, mgnify
from .clients.http import batches
from .constants import (
    BIN_LAYER,
    CONTIG_LAYER,
    EXCLUDED_TAXA,
    KINGDOM_TAXA,
    LONG_READ_PLATFORMS,
    MAG_LAYER,
    METAGENOMIC_WGS,
    MIN_READS,
)
from .log import warn
from .models import EnaAssembly, EnaBin, EnaMag, EnaRun, EnaSample, EnaSearch, MgnifyCatalogue

Config = dict[str, Any]  # a biome YAML (biomes/<name>.yaml)
Record = dict[str, Any]  # one survey record: an ENA row (models.py) plus the survey's own columns

KINGDOM_TERMS = {
    "prokaryotes": re.compile(r"bacteri|archae|prokaryot|\b16s\b", re.I),
    "eukaryotes": re.compile(
        r"eukaryot|protist|fung(?:i|al)|\balga|diatom|phytoplankton|dinoflagell|ciliat"
        r"|haptophyt|chlorophyt|oomycet|\byeast|amoeb|foraminifer|radiolar|\b18s\b",
        re.I,
    ),
}

SURVEY_COLUMNS = [
    "biome",
    "kingdom",
    "metagenome_taxon",
    "tax_id",
    "catalogue",
    "catalogue_updated",
    "catalogue_species",
    "catalogue_genomes",
    "catalogue_mags",
    "mgnify_samples",
    "mags",
    "mags_new",
    "mags_mgnify",
    "mags_to_add",
    "mag_studies",
    "mgnify_bins",
    "mgnify_bins_new",
    "mgnify_bins_to_add",
    "contigs",
    "contigs_mgnify",
    "contigs_kingdom",
    "contigs_to_bin",
    "contig_studies",
    "runs",
    "runs_1m",
    "runs_1m_kingdom",
    "runs_1m_new",
    "runs_1m_kingdom_new",
    "runs_to_assemble",
    "runs_1m_paired",
    "runs_1m_long",
    "gbases_1m",
    "gbases_1m_kingdom",
    "gbases_to_assemble",
    "read_studies_1m",
]
RECORD_COLUMNS = [
    "data_type",
    "accession",
    "study_accession",
    "from_mgnify",
    "sample_accession",
    "metagenome_taxon",
    "first_public",
    "is_new",
    "read_count",
    "base_count",
    "library_layout",
    "instrument_platform",
    "center_name",
    "assembly_software",
    "kingdom_evidence",
    "kingdom_match",
    "in_catalogue",
    "route",
    "scope_taxon",
    "tax_id",
    "scientific_name",
    "source_samples",
    "host",
    "country",
    "lat",
    "lon",
    "collection_date",
]
STUDY_COLUMNS = [
    "study_accession",
    "from_mgnify",
    "metagenome_taxa",
    "samples",
    "first_public",
    "last_public",
    "mags",
    "mags_new",
    "mags_mgnify",
    "mags_to_add",
    "mgnify_bins",
    "mgnify_bins_new",
    "mgnify_bins_to_add",
    "contigs",
    "contigs_mgnify",
    "contigs_kingdom",
    "contigs_to_bin",
    "runs",
    "runs_1m",
    "runs_1m_kingdom",
    "runs_1m_kingdom_new",
    "runs_to_assemble",
    "runs_1m_new",
    "runs_1m_paired",
    "runs_1m_long",
    "gbases_1m",
    "gbases_1m_kingdom",
    "gbases_to_assemble",
]


def summarize(
    mags: Sequence[Record], contigs: Sequence[Record], runs: Sequence[Record], cutoff: str, bins: Sequence[Record] = ()
) -> dict[str, int | float | str]:
    """Counts for one set of records. `cutoff` (YYYY-MM-DD or "") splits off records public after it."""

    def new(rows: Iterable[Record]) -> int | str:
        return sum(r["first_public"] >= cutoff for r in rows) if cutoff else ""

    def gbases(rows: Iterable[Record]) -> float:
        return round(sum(r["base_count"] or 0 for r in rows) / 1e9, 1)

    big = [r for r in runs if (r["read_count"] or 0) >= MIN_READS]
    long_ = [r for r in big if r["instrument_platform"] in LONG_READ_PLATFORMS]
    to_assemble = [r for r in big if r.get("route") == "assemble"]
    return {
        "mags": len(mags),
        "mags_new": new(mags),
        "mags_mgnify": sum(r.get("from_mgnify", False) for r in mags),
        "mags_to_add": sum(not r.get("in_catalogue", False) for r in mags),
        "mag_studies": len({r["study_accession"] for r in mags}),
        "mgnify_bins": len(bins),
        "mgnify_bins_new": new(bins),
        "mgnify_bins_to_add": sum(not r.get("in_catalogue", False) for r in bins),
        "contigs": len(contigs),
        "contigs_mgnify": sum(r.get("from_mgnify", False) for r in contigs),
        "contigs_kingdom": sum(r.get("kingdom_match", True) for r in contigs),
        "contigs_to_bin": sum(r.get("route") == "bin" for r in contigs),
        "contig_studies": len({r["study_accession"] for r in contigs}),
        "runs": len(runs),
        "runs_1m": len(big),
        "runs_1m_kingdom": sum(r.get("kingdom_match", True) for r in big),
        "runs_1m_new": new(big),
        "runs_1m_kingdom_new": new([r for r in big if r.get("kingdom_match", True)]),
        "runs_to_assemble": len(to_assemble),
        "gbases_to_assemble": gbases(to_assemble),
        "runs_1m_paired": sum(r["library_layout"] == "PAIRED" for r in big)
        - sum(r["library_layout"] == "PAIRED" for r in long_),
        "runs_1m_long": len(long_),
        "gbases_1m": gbases(big),
        "gbases_1m_kingdom": gbases([r for r in big if r.get("kingdom_match", True)]),
        "read_studies_1m": len({r["study_accession"] for r in big}),
    }


def host_filter(cfg: Config) -> str:
    """ENA query clause for host-associated biomes. Hosts are recorded inconsistently, so it
    matches the host tax ID *or* any keyword in the free-text `host` field."""
    terms = [f"host_tax_id={cfg['host_tax_id']}"] if cfg.get("host_tax_id") else []
    terms += [f'host="*{k}*"' for k in cfg.get("host_keywords") or []]
    return f" AND ({' OR '.join(terms)})" if terms else ""


def host_matches(sample: EnaSample, cfg: Config) -> bool:
    """Python twin of host_filter() for one ENA sample (case-insensitive, like ENA)."""
    if cfg.get("host_tax_id") and sample.host_tax_id == str(cfg["host_tax_id"]):
        return True
    return any(k.lower() in sample.host.lower() for k in cfg.get("host_keywords") or [])


def by_id(samples: Iterable[EnaSample]) -> dict[str, EnaSample]:
    """Samples keyed by both their primary and secondary IDs."""
    return {i: s for s in samples for i in s.ids}


def host_mags(mags: list[Record], cfg: Config) -> list[Record]:
    """MAGs of a host-associated biome, kept only if every source sample (`derived from`)
    matches the host filter: ENA MAG samples don't carry their source sample's host, and a
    co-assembled MAG counts only if all its sources are in the biome. Fills `source_samples`
    and `host` on every MAG record."""
    origin = ena.derived_from(m["sample_accession"] for m in mags).data
    hosts = by_id(ena.sample_hosts({s for src in origin.values() for s in src}).data)
    kept, unresolved = [], 0
    for m in mags:
        src = origin.get(m["sample_accession"], [])
        rows = [hosts.get(s) for s in src]
        m["source_samples"] = ";".join(src)
        m["host"] = ";".join(sorted({r.host for r in rows if r and r.host}))
        unresolved += not src or None in rows
        if src and all(r and host_matches(r, cfg) for r in rows):
            kept.append(m)
    if unresolved:
        warn(f"{cfg['name']}: {unresolved} of {len(mags)} MAGs have no resolvable source sample")
    return kept


def taxa_filter(cfg: Config) -> str:
    """ENA clause removing out-of-scope genome taxa (defaults plus `exclude_taxa_extra`)."""
    return "".join(f" AND NOT tax_tree({t})" for t in EXCLUDED_TAXA + (cfg.get("exclude_taxa_extra") or []))


def tag(
    rows: Iterable[EnaMag | EnaAssembly | EnaBin | EnaRun],
    data_type: str,
    acc_field: str,
    taxon: str,
    cutoff: str,
    **extra: Any,
) -> list[Record]:
    """ENA rows as survey records (see RECORD_COLUMNS)."""
    return [
        {
            **asdict(r),
            **extra,
            "data_type": data_type,
            "accession": getattr(r, acc_field),
            "metagenome_taxon": taxon,
            "is_new": (r.first_public >= cutoff) if cutoff else "",
        }
        for r in rows
    ]


def mgnify_bins(cfg: Config, scope: dict[int, str], records: list[Record], cutoff: str) -> list[Record]:
    """Genomes MGnify binned from this biome's data. MGnify deposits them as `binned metagenome`
    analyses in its own ENA studies, so they never show up as MAGs. A bin counts if every
    source sample (`derived from`) is a sample of this biome's assemblies or runs, which carries
    the metagenome taxon and host rules over; it takes the taxon of its first source sample."""
    taxon_of = {r["sample_accession"]: r["metagenome_taxon"] for r in records}
    studies = [r["study_accession"] for r in records if r["from_mgnify"]]
    bins = []
    for batch in batches(studies, 50):
        in_studies = " OR ".join(f'study_accession="{a}"' for a in batch)
        for t, name in scope.items():
            query = f"{BIN_LAYER} AND tax_tree({t}){taxa_filter(cfg)} AND ({in_studies})"
            bins += tag(
                ena.search(EnaSearch("analysis", query), EnaBin).data,
                "bin",
                "analysis_accession",
                "",
                cutoff,
                scope_taxon=name,
                from_mgnify=True,
            )
    origin = ena.derived_from(b["sample_accession"] for b in bins).data
    primary = {
        i: s.sample_accession
        for i, s in by_id(ena.sample_hosts({s for src in origin.values() for s in src}).data).items()
    }
    kept = []
    for b in bins:
        src = [primary.get(s, s) for s in origin.get(b["sample_accession"], [])]
        if src and all(s in taxon_of for s in src):
            b["source_samples"] = ";".join(src)
            b["metagenome_taxon"] = taxon_of[src[0]]
            kept.append(b)
    return kept


def study_kingdom(text: str, kingdom: str) -> str:
    """`study_text` if the study's title or description names the config's kingdom
    (`KINGDOM_TERMS`), `other_kingdom` if it names only the other one, else `none`. Assemblies
    and runs are whole communities with no organism, so this and the genomes already recovered
    from their sample are the survey's only kingdom signal for them."""
    other = next(k for k in KINGDOM_TERMS if k != kingdom)
    if KINGDOM_TERMS[kingdom].search(text):
        return "study_text"
    return "other_kingdom" if KINGDOM_TERMS[other].search(text) else "none"


def genome_samples(genomes: list[Record]) -> set[str]:
    """Source samples (`derived from`) of MAGs and bins, under both their primary and
    secondary IDs. Fills `source_samples` on genomes that lack it."""
    origin = ena.derived_from(g["sample_accession"] for g in genomes if not g.get("source_samples")).data
    for g in genomes:
        if not g.get("source_samples"):
            g["source_samples"] = ";".join(origin.get(g["sample_accession"], []))
    src = {s for g in genomes for s in g["source_samples"].split(";") if s}
    return src | {i for s in ena.sample_hosts(src).data for i in s.ids}


def kingdom_evidence(cfg: Config, genomes: list[Record], contigs: list[Record], runs: list[Record]) -> None:
    """Label each assembly and run with its evidence for the config's kingdom, strongest first:
    `genome` (its sample already yielded a MAG or bin of that kingdom), `study_text`,
    `other_kingdom`, `none`. `kingdom_match`: eukaryotes need `genome` or `study_text`;
    prokaryotes, being everywhere, count unless the study is about eukaryotes only.

    Then `route`, the step each counted record would enter at (a sample enters at its most
    processed data): assemblies whose sample has no genome yet are to `bin`; runs ≥ 1 M reads
    whose sample has neither a genome nor an assembly are to `assemble`; the rest are `covered`."""
    texts = {s.study_accession: s.text for s in ena.studies({r["study_accession"] for r in contigs + runs}).data}
    from_genomes = genome_samples(genomes)
    for r in contigs + runs:
        r["kingdom_evidence"] = (
            "genome"
            if r["sample_accession"] in from_genomes
            else study_kingdom(texts.get(r["study_accession"], ""), cfg["kingdom"])
        )
        r["kingdom_match"] = (
            r["kingdom_evidence"] in ("genome", "study_text")
            if cfg["kingdom"] == "eukaryotes"
            else r["kingdom_evidence"] != "other_kingdom"
        )
    assembled = {r["sample_accession"] for r in contigs if r["kingdom_match"]}
    for r in contigs:
        r["route"] = "" if not r["kingdom_match"] else "covered" if r["sample_accession"] in from_genomes else "bin"
    for r in runs:
        covered = r["sample_accession"] in from_genomes or r["sample_accession"] in assembled
        r["route"] = (
            "" if not r["kingdom_match"] or (r["read_count"] or 0) < MIN_READS else "covered" if covered else "assemble"
        )


def catalogue_membership(cat: MgnifyCatalogue | None, genomes: list[Record]) -> int | str:
    """Flag MAGs and bins already in the MGnify catalogue `cat` (matched by accession or by the
    genome's own sample). Returns the catalogue's MAG count ("" without a catalogue)."""
    listed = mgnify.catalogue_genomes(cat).data if cat else []
    ids = {i for g in listed for i in (g.genome_accession, g.sample_accession) if i}
    for g in genomes:
        g["in_catalogue"] = g["accession"] in ids or g["sample_accession"] in ids
    return sum(g.genome_type == "MAG" for g in listed) if cat else ""


def under(lineage: str, parent: str) -> bool:
    """True if `lineage` is `parent` or one of its descendants."""
    return lineage == parent or lineage.startswith(parent + ":")


def in_scope(lineage: str, include: list[str], exclude: list[str]) -> bool:
    """Biome scope rule: whole subtrees of `mgnify_lineages`, minus whole subtrees of
    `mgnify_lineages_exclude`."""
    return any(under(lineage, p) for p in include) and not any(under(lineage, p) for p in exclude)


def mgnify_samples(cfg: Config) -> int:
    """MGnify-analysed samples in the biome's lineage scope."""
    include, exclude = cfg.get("mgnify_lineages") or [], cfg.get("mgnify_lineages_exclude") or []
    nodes: dict[str, int] = {}
    for lineage in include + exclude:
        subtree = mgnify.biome_subtree(lineage)
        if not subtree.ok:
            warn(f"{cfg['name']}: MGnify lineage {lineage!r} does not exist")
        nodes.update({b.lineage: b.samples_count for b in subtree.data})
    return sum(n for lin, n in nodes.items() if in_scope(lin, include, exclude))


def survey_biome(cfg: Config) -> tuple[list[dict[str, Any]], list[Record], str]:
    """One survey row per metagenome taxon of a biome config, plus a TOTAL row; the records
    behind them; and the cutoff date. MAGs are searched once per kingdom taxon, so each is
    tagged with its group; a record is `from_mgnify` if it sits in one of MGnify's ENA studies."""
    cat = mgnify.catalogue(cfg["mgnify_catalogue"]).data if cfg.get("mgnify_catalogue") else None
    if cfg.get("mgnify_catalogue") and not cat:
        raise ValueError(f"{cfg['name']}: MGnify catalogue {cfg['mgnify_catalogue']!r} not found")
    if cat and cat.catalogue_type != cfg["kingdom"]:
        raise ValueError(
            f"{cfg['name']}: {cfg['mgnify_catalogue']} is a {cat.catalogue_type} "
            f"catalogue, but the config's kingdom is {cfg['kingdom']}"
        )
    cutoff = cat.updated_at[:10] if cat else ""
    base = {
        "biome": cfg["name"],
        "kingdom": cfg["kingdom"],
        "catalogue": cfg.get("mgnify_catalogue", ""),
        "catalogue_updated": cutoff,
        "catalogue_species": cat.genome_count if cat else "",
        "catalogue_genomes": cat.unclustered_genome_count if cat else "",
        "mgnify_samples": mgnify_samples(cfg),
    }
    run_host = host_filter(cfg)

    scope = {}
    for t in KINGDOM_TAXA[cfg["kingdom"]]:
        found = ena.taxon(t).data
        scope[t] = found.scientific_name if found else str(t)
    mgnify_studies = ena.emg_studies().data
    per_taxon: dict[str, tuple[str, list[Record], list[Record], list[Record]]] = {}
    all_mags: list[Record] = []
    all_contigs: list[Record] = []
    all_runs: list[Record] = []
    for taxon in cfg["metagenome_taxa"]:
        found = ena.taxon_by_name(taxon).data
        if not found:
            warn(f"{cfg['name']}: unknown ENA taxon {taxon!r}, skipped")
            continue
        taxid = found.tax_id
        mags = [
            m
            for t, name in scope.items()
            for m in tag(
                ena.search(
                    EnaSearch(
                        "wgs_set", f'{MAG_LAYER} AND metagenome_source="{taxon}" AND tax_tree({t}){taxa_filter(cfg)}'
                    ),
                    EnaMag,
                ).data,
                "mag",
                "accession",
                taxon,
                cutoff,
                scope_taxon=name,
            )
        ]
        if run_host:
            mags = host_mags(mags, cfg)
        contigs = tag(
            ena.search(EnaSearch("analysis", f"{CONTIG_LAYER} AND tax_eq({taxid}){run_host}"), EnaAssembly).data,
            "contigs",
            "analysis_accession",
            taxon,
            cutoff,
        )
        runs = tag(
            ena.search(EnaSearch("read_run", f"tax_eq({taxid}) AND {METAGENOMIC_WGS}{run_host}"), EnaRun).data,
            "reads",
            "run_accession",
            taxon,
            cutoff,
        )
        for r in mags + contigs + runs:
            r["from_mgnify"] = r["study_accession"] in mgnify_studies
        per_taxon[taxon] = (taxid, mags, contigs, runs)
        all_mags += mags
        all_contigs += contigs
        all_runs += runs
    bins = mgnify_bins(cfg, scope, all_contigs + all_runs, cutoff)
    kingdom_evidence(cfg, all_mags + bins, all_contigs, all_runs)
    base["catalogue_mags"] = catalogue_membership(cat, all_mags + bins)
    rows = [
        {
            **base,
            "metagenome_taxon": taxon,
            "tax_id": taxid,
            **summarize(mags, contigs, runs, cutoff, [b for b in bins if b["metagenome_taxon"] == taxon]),
        }
        for taxon, (taxid, mags, contigs, runs) in per_taxon.items()
    ]
    rows.append(
        {**base, "metagenome_taxon": "TOTAL", "tax_id": "", **summarize(all_mags, all_contigs, all_runs, cutoff, bins)}
    )
    return rows, all_mags + bins + all_contigs + all_runs, cutoff


def studies(records: list[Record], cutoff: str) -> list[dict[str, Any]]:
    """One row per study: its counts, samples and publication date range."""
    by_study = defaultdict(list)
    for r in records:
        by_study[r["study_accession"]].append(r)
    rows = []
    for study, recs in by_study.items():

        def of(data_type: str) -> list[Record]:
            return [r for r in recs if r["data_type"] == data_type]

        dates = [r["first_public"] for r in recs]
        rows.append(
            {
                "study_accession": study,
                "from_mgnify": recs[0].get("from_mgnify", False),
                "metagenome_taxa": ";".join(sorted({r["metagenome_taxon"] for r in recs})),
                "samples": len({r["sample_accession"] for r in recs if r["data_type"] in ("contigs", "reads")}),
                "first_public": min(dates),
                "last_public": max(dates),
                **summarize(of("mag"), of("contigs"), of("reads"), cutoff, of("bin")),
            }
        )
    return sorted(rows, key=lambda r: -len(by_study[r["study_accession"]]))


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run(biome_paths: Iterable[str], out_dir: str | Path) -> None:
    """Survey every biome config into `out_dir`: survey.tsv, biomes/<name>/{records,studies}.tsv
    and the dashboard report.html."""
    out = Path(out_dir)
    rows = []
    for path in biome_paths:
        with open(path) as f:
            cfg = yaml.safe_load(f)
        if cfg.get("kingdom") not in KINGDOM_TAXA:
            raise ValueError(f"{path}: kingdom must be one of {', '.join(KINGDOM_TAXA)}")
        biome_rows, records, cutoff = survey_biome(cfg)
        rows += biome_rows
        biome_dir = out / "biomes" / cfg["name"]
        biome_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(biome_dir / "records.tsv", RECORD_COLUMNS, records)
        write_tsv(biome_dir / "studies.tsv", STUDY_COLUMNS, studies(records, cutoff))
    write_tsv(out / "survey.tsv", SURVEY_COLUMNS, rows)
    report.write(out)
