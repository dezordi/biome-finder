import csv

import pytest

from biome_finder import score, search
from biome_finder.clients import ena, ncbi
from biome_finder.models import NCBI_RANKS, EnaSampleMetadata, NcbiTaxdump, NcbiTaxon, Response


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
    monkeypatch.setattr(ena, "sample_metadata", lambda ids: Response(200, []))


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


def meta(acc, title="", description="", **kw):
    fields = dict.fromkeys(EnaSampleMetadata.__dataclass_fields__, "")
    return EnaSampleMetadata(
        **{**fields, "sample_accession": acc, "sample_title": title, "sample_description": description, **kw}
    )


@pytest.mark.parametrize(
    "text, flags",
    [  # real freshwater-water-prok sample texts (2026-11)
        ("blank.extraction2 urban biome research facility", ["control"]),
        ("GSC:MIxS MIMS:6.0 negative control home", ["control"]),
        ("Overlying water from a 30-day control microcosm without apple snails.", []),  # an experiment's arm
        ("mock sample for coassembly ErkenSummer", []),  # ENA's virtual co-assembly sample
        ("mock community DNA standard", ["mock"]),
        ("Synthetic microbial communities from Ohio, USA", ["synthetic"]),
        ("shotgun sequenced culture enrichment of surface water microbiome", ["enrichment"]),
        ("Enriched microbiome from surface water small river", ["enrichment"]),
        ("Walnut Creek spiked with Pseudomonas putida & Escherichia coli", ["spiked"]),
        ("Lake Erie surface water", []),
    ],
)
def test_hard_filters(text, flags):
    assert score.hard_flags(meta("S1", text)) == flags


def test_sample_evidence_fills_reported_biome_and_sends_hard_filters_to_tier_x():
    rows = [
        {"record_accession": "M1", "sample_accession": "SA;SB", "flags": ""},
        {"record_accession": "R1", "sample_accession": "SA", "flags": "low_reads"},
        {"record_accession": "R2", "sample_accession": "SC", "flags": ""},
    ]
    metadata = [
        meta(
            "SA",
            "lake water",
            broad_scale_environmental_context="freshwater lake biome",
            scientific_name="freshwater metagenome",
        ),
        meta("SB", "lake water", "spiked with E. coli", broad_scale_environmental_context="river biome"),
        meta("SC", secondary_sample_accession="ERS9", environmental_medium="water"),
    ]
    search.sample_evidence(rows, metadata)
    m1, r1, r2 = rows
    assert m1["reported_biome"] == (
        "env_broad_scale=freshwater lake biome | river biome; organism=freshwater metagenome"
    )
    assert m1["flags"] == "spiked" and m1["tier"] == "X"  # one source sample is enough
    assert r1["flags"] == "low_reads" and "tier" not in r1
    assert r2["reported_biome"] == "env_medium=water"
