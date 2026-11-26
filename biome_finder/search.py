"""`search`: build samples.tsv for one biome, starting from its survey records."""

import csv
from pathlib import Path
from typing import Any

import yaml

from . import score
from .clients import ena, ncbi
from .constants import EXCLUDED_TAXA, KINGDOM_TAXA, LONG_READ_PLATFORMS
from .models import EnaSampleMetadata, NcbiTaxdump, NcbiTaxon
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
        records = [r for r in csv.DictReader(f, delimiter="\t") if r["kingdom_match"] != "False"]
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


def to_row(r: dict[str, str]) -> Row:
    """A survey record as a samples.tsv row. A genome's biome is its source samples' biome, so
    `sample_accession` holds those (`;`-joined for co-assemblies) when ENA links them."""
    reads = r["data_type"] == "reads"
    return {
        "record_accession": r["accession"],
        "data_type": r["data_type"],
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


def sample_evidence(rows: list[Row], metadata: list[EnaSampleMetadata]) -> None:
    """Fill `reported_biome` from each row's samples (a genome's source samples), add the hard
    filters their text matches to `flags` (score.HARD_FILTERS), and put rows carrying an
    excluding flag in tier X. A genome is flagged if any of its source samples is."""
    by_id = {i: s for s in metadata for i in s.ids}
    for r in rows:
        samples = [by_id[a] for a in r["sample_accession"].split(";") if a in by_id]
        r["reported_biome"] = reported_biome(samples)
        hard = dict.fromkeys(f for s in samples for f in score.hard_flags(s))
        flags = [f for f in r["flags"].split(";") if f] + list(hard)
        r["flags"] = ";".join(flags)
        if score.EXCLUDING_FLAGS & set(flags):
            r["tier"] = "X"


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
    rows = [to_row(r) for r in survey_records(Path(survey_dir), cfg["name"]) if r["first_public"] >= since]
    release = f"{taxdump[:7]}-01" if taxdump else ncbi.latest_taxdump().data
    dump = ncbi.taxdump(release, {r["reported_tax_id"] for r in rows}).data
    if not dump:
        raise ValueError(f"NCBI taxdump release {release!r} not found in the archive")
    resolve_taxonomy(rows, cfg, dump)
    samples = {a for r in rows for a in r["sample_accession"].split(";") if a}
    sample_evidence(rows, ena.sample_metadata(samples).data)
    out = Path(out)
    old: list[Row] = []
    if out.exists():
        with open(out, newline="") as f:
            old = list(csv.DictReader(f, delimiter="\t"))
    write_tsv(out, SAMPLES_COLUMNS, merge(old, rows))
