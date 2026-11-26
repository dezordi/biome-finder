import csv

import pytest

from biome_finder import search
from biome_finder.clients import ncbi
from biome_finder.models import NCBI_RANKS, NcbiTaxdump, NcbiTaxon, Response


def taxon(tax_id, name, rank, ancestors, **ranks):
    fields = {"class_" if r == "class" else r: ranks.get(r) for r in NCBI_RANKS}
    return NcbiTaxon(tax_id, name, rank, **fields, ancestors=frozenset(ancestors))


DUMP = NcbiTaxdump(
    "2026-11-01",
    merged={"99": "1234"},
    deleted=frozenset({"77"}),
    taxa={
        "1234": taxon(
            "1234", "Nitrospira", "genus", {"131567", "2"}, domain="Bacteria", phylum="Nitrospirota", genus="Nitrospira"
        ),
        "5555": taxon("5555", "Nitzschia sp.", "species", {"131567", "2759"}, domain="Eukaryota"),
        "9999": taxon("9999", "phage X", "species", {"10239"}, domain="Duplodnaviria"),
        "1": taxon("1", "freshwater metagenome", "species", {"408169"}, species="freshwater metagenome"),
    },
)


@pytest.fixture(autouse=True)
def taxdump(monkeypatch):
    monkeypatch.setattr(ncbi, "latest_taxdump", lambda: Response(200, "2026-11-01"))
    monkeypatch.setattr(ncbi, "taxdump", lambda release, ids: Response(200, DUMP))


FIELDS = ["data_type", "accession", "study_accession", "sample_accession", "source_samples", "first_public"]
FIELDS += ["kingdom_match", "in_catalogue", "tax_id", "scientific_name", "read_count", "library_layout"]
FIELDS += ["instrument_platform", "metagenome_taxon"]


def record(acc, data_type="reads", kingdom_match="True", **kw):
    base = dict.fromkeys(FIELDS, "")
    base.update(data_type=data_type, accession=acc, sample_accession="S" + acc, first_public="2024-01-01")
    base.update(kingdom_match=kingdom_match, read_count="5000000", library_layout="PAIRED")
    base.update(instrument_platform="ILLUMINA", metagenome_taxon="freshwater metagenome", **kw)
    return base


def write_survey(tmp_path, records):
    d = tmp_path / "survey" / "biomes" / "lake"
    d.mkdir(parents=True)
    with open(d / "records.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(records)
    (tmp_path / "lake.yaml").write_text("name: lake\nkingdom: prokaryotes\n")
    (tmp_path / "survey" / "survey.tsv").write_text(
        "biome\tmetagenome_taxon\ttax_id\nlake\tfreshwater metagenome\t1\nlake\tTOTAL\t\nriver\tfreshwater metagenome\t2\n"
    )


def read(path):
    with open(path, newline="") as f:
        return {r["record_accession"]: r for r in csv.DictReader(f, delimiter="\t")}


def test_search_keeps_kingdom_records_and_curation_on_rerun(tmp_path):
    write_survey(
        tmp_path,
        [
            record("R1"),
            record("R2", kingdom_match="False"),
            record("R3", read_count="1000", library_layout="SINGLE"),
            record("M1", data_type="mag", kingdom_match="", source_samples="SA;SB", in_catalogue="True", tax_id="1234"),
        ],
    )
    out = tmp_path / "samples.tsv"
    search.run(str(tmp_path / "lake.yaml"), tmp_path / "survey", out)
    rows = read(out)
    assert set(rows) == {"R1", "R3", "M1"}
    assert list(next(iter(rows.values()))) == search.SAMPLES_COLUMNS
    assert rows["R3"]["flags"] == "single_end;low_reads"
    assert rows["R1"]["reported_tax_id"] == "1" and rows["R1"]["tax_lineage"].endswith("s__freshwater metagenome")
    assert rows["M1"]["sample_accession"] == "SA;SB" and rows["M1"]["flags"] == "in_existing_catalogue"

    curated = {**rows["R1"], "decision": "skip", "curation_note": "control"}
    gone = {**rows["R3"], "record_accession": "R9"}
    search.write_tsv(out, search.SAMPLES_COLUMNS, [curated, gone])
    search.run(str(tmp_path / "lake.yaml"), tmp_path / "survey", out)
    rows = read(out)
    assert rows["R1"]["decision"] == "skip" and rows["R1"]["curation_note"] == "control"
    assert set(rows) == {"R1", "R3", "M1", "R9"}


def test_search_since_filters_by_first_public(tmp_path):
    write_survey(tmp_path, [record("R1"), record("R2", first_public="2026-11-01")])
    out = tmp_path / "samples.tsv"
    search.run(str(tmp_path / "lake.yaml"), tmp_path / "survey", out, since="2026-01-01")
    assert set(read(out)) == {"R2"}


def test_taxonomy_is_resolved_and_genomes_out_of_scope_are_flagged():
    rows = [
        {"record_accession": a, "data_type": t, "reported_tax_id": tax, "flags": ""}
        for a, t, tax in [
            ("M_merged", "mag", "99"),
            ("M_euk", "mag", "5555"),
            ("M_virus", "mag", "9999"),
            ("M_deleted", "mag", "77"),
            ("M_unknown", "mag", "42"),
            ("R1", "reads", "1"),
        ]
    ]
    search.resolve_taxonomy(rows, {"kingdom": "prokaryotes"}, DUMP)
    by = {r["record_accession"]: r for r in rows}
    assert by["M_merged"]["domain"] == "prokaryotes"  # through the merged ID
    assert by["M_merged"]["flags"] == "" and by["M_merged"]["taxdump_date"] == "2026-11-01"
    assert by["M_merged"]["tax_lineage"] == "d__Bacteria;k__;p__Nitrospirota;c__;o__;f__;g__Nitrospira;s__"
    assert by["M_euk"]["flags"] == "out_of_scope_taxon"  # the other kingdom
    assert by["M_virus"]["flags"] == "out_of_scope_taxon"
    assert by["M_deleted"]["flags"] == "taxid_deleted" and by["M_deleted"]["tax_lineage"] == ""  # a curator decides
    assert by["M_unknown"]["flags"] == "taxid_unknown"
    assert by["R1"]["flags"] == "" and by["R1"]["domain"] == ""  # no scope for runs
    assert by["R1"]["tax_lineage"] == "d__;k__;p__;c__;o__;f__;g__;s__freshwater metagenome"
