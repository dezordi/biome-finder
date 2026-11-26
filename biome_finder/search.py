"""`search`: build samples.tsv for one biome, starting from its survey records."""

import csv
from pathlib import Path
from typing import Any

import yaml

from .constants import LONG_READ_PLATFORMS
from .survey import write_tsv

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
    "tax_id",
    "tax_name",
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


def survey_records(survey_dir: Path, name: str) -> list[dict[str, str]]:
    """The biome's survey records (survey/biomes/<name>/records.tsv), minus assemblies and runs
    whose sample has no evidence for the config's kingdom (`kingdom_match` false)."""
    path = survey_dir / "biomes" / name / "records.tsv"
    if not path.exists():
        raise ValueError(f"{path} not found: run `biome-finder survey` on this biome first")
    with open(path, newline="") as f:
        return [r for r in csv.DictReader(f, delimiter="\t") if r["kingdom_match"] != "False"]


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


def merge(old: list[Row], new: list[Row]) -> list[Row]:
    """New rows with the human and `prepare` columns of the old row of the same record. Old rows
    the new search no longer finds are kept unchanged, so no curation is lost."""
    by_acc = {r["record_accession"]: r for r in old}
    merged = []
    for r in new:
        prev = by_acc.pop(r["record_accession"], {})
        merged.append({**r, **{k: prev[k] for k in KEPT_COLUMNS if prev.get(k)}})
    return merged + list(by_acc.values())


def run(biome_path: str, survey_dir: str | Path, out: str | Path, since: str = "") -> None:
    """Write `out` (samples.tsv) for one biome config, merged into `out` if it exists."""
    with open(biome_path) as f:
        cfg = yaml.safe_load(f)
    rows = [to_row(r) for r in survey_records(Path(survey_dir), cfg["name"]) if r["first_public"] >= since]
    out = Path(out)
    old: list[Row] = []
    if out.exists():
        with open(out, newline="") as f:
            old = list(csv.DictReader(f, delimiter="\t"))
    write_tsv(out, SAMPLES_COLUMNS, merge(old, rows))
