import csv

import pytest

from biome_finder import score, search
from biome_finder.clients import ena, ncbi, ols
from biome_finder.models import NCBI_RANKS, EnaSampleMetadata, EnaStudy, NcbiTaxdump, NcbiTaxon, OlsTerm, Response


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
    monkeypatch.setattr(ena, "studies", lambda ids: Response(200, []))
    monkeypatch.setattr(ols, "with_descendants", lambda term: Response(200, []))


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


def test_sample_evidence_fills_reported_biome_and_sends_hard_filters_to_tier_d():
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
    assert m1["flags"] == "spiked" and m1["tier"] == "D"  # one source sample is enough
    assert r1["flags"] == "low_reads" and "tier" not in r1
    assert r2["reported_biome"] == "env_medium=water"


LAKE = score.Terms.of(
    [
        OlsTerm("ENVO:00000021", "freshwater lake", ("FreshwaterLake",)),
        OlsTerm("ENVO:01000252", "freshwater lake biome", ()),
    ]
)
SEDIMENT = score.Terms.of([OlsTerm("ENVO:00002007", "sediment", ()), OlsTerm("ENVO:00000546", "lake sediment", ())])
CFG = {
    "catalogue_biome": "root:Environmental:Aquatic:Freshwater",
    "metagenome_taxa": ["freshwater metagenome"],
    "keywords_include": ["freshwater", "lake water"],
    "keywords_exclude": ["marine", "mock", "sediment"],
}


@pytest.mark.parametrize(
    "value, match",
    [
        ("freshwater lake biome [ENVO:01000252]", True),
        ("freshwater lake biome[ENVO_01000252]", True),  # as submitters write it
        ("Freshwater Lake", True),  # the name alone
        ("ENVO:1000252", True),  # a dropped leading zero
        ("saline lake", False),
        ("freshwater lake water", False),  # names match exactly, never as substrings
    ],
)
def test_envo_terms_match_by_id_or_exact_name(value, match):
    assert bool(LAKE.match(value)) is match


def test_an_included_label_beats_an_excluded_synonym():
    pond = score.Terms.of([OlsTerm("ENVO:00000055", "saline evaporation pond", ("lake",))])
    lake = score.Terms.of([OlsTerm("ENVO:00000020", "lake", ())])
    items = score.sample_items(meta("S1", local_environmental_context="lake"), lake, pond, None, None, set())
    assert [e.claim for e in items] == ["in"]
    items = score.sample_items(meta("S1", local_environmental_context="salt pond"), lake, pond, None, None, set())
    assert items == []  # neither a label nor a synonym of the included terms


def test_combine_and_tiers():
    envo = score.Evidence("sample_metadata", "env_local_scale", "in", 1.0, score.W_ENVO, "")
    name = score.Evidence("metagenome_name", "organism", "in", 1.0, score.W_METAGENOME, "")
    marine = score.Evidence("sample_metadata", "text", "out", 1.0, score.W_SAMPLE_EXCLUDE, "marine")
    assert score.combine([envo, name]) == (86, False) and score.tier(86, False) == "A"  # 1 − 0.2·0.7
    assert score.combine([name]) == (30, False) and score.tier(30, False) == "C"  # a default name alone
    assert score.combine([name], siblings=1.0) == (40, False)
    value, conflict = score.combine([envo, name, marine])
    assert conflict and score.tier(value, conflict) == "C"


def run_biome(metadata, rows, texts=None):
    search.biome_evidence(rows, CFG, metadata, texts or {}, LAKE, SEDIMENT)
    return {r["record_accession"]: r for r in rows}


def row(acc, sample, study="P1", data_type="reads"):
    return {
        "record_accession": acc,
        "data_type": data_type,
        "sample_accession": sample,
        "study_accession": study,
        "flags": "",
    }


def test_biome_evidence_scores_samples_and_genomes_take_their_worst_source_sample():
    metadata = [
        meta(
            "S1",
            "Metagenome or environmental sample from freshwater metagenome",
            local_environmental_context="freshwater lake",
            scientific_name="freshwater metagenome",
        ),
        meta(
            "S2",
            "Metagenome or environmental sample from freshwater metagenome",
            scientific_name="freshwater metagenome",
        ),
        meta("S3", "mock sample for coassembly lake water", scientific_name="freshwater metagenome"),
        meta("S4", "lake sediment core", environmental_medium="lake sediment [ENVO:00000546]"),
        meta(
            "S5", "marine water", local_environmental_context="freshwater lake", scientific_name="freshwater metagenome"
        ),
    ]
    rows = [
        row("R1", "S1"),
        row("R2", "S2"),
        row("R3", "S3"),
        row("R4", "S4"),
        row("R5", "S5"),
        row("M1", "S1;S5", study="P2", data_type="mag"),
    ]
    by = run_biome(metadata, rows, {"P1": "Freshwater lakes of Sweden", "P2": "freshwater metagenome MAGs"})
    assert by["R1"]["tier"] == "A" and by["R1"]["curated_biome"] == CFG["catalogue_biome"]
    assert "sample_metadata:env_local_scale:in:freshwater lake" in by["R1"]["evidence"]
    assert "sample_metadata:text" not in by["R2"]["evidence"]  # a default title only repeats the taxon
    assert "out:mock" not in by["R3"]["evidence"]  # ENA's co-assembly sample, no mock community
    assert by["R4"]["tier"] == "D" and "excluded_envo" in by["R4"]["flags"].split(";")
    assert by["R5"]["tier"] == "C" and by["R5"]["needs_curation"]  # a conflict
    assert by["R5"]["flags"] == "off_biome_keyword"  # the reason, visible
    assert by["M1"]["tier"] == "C" and by["M1"]["needs_curation"]  # S1 is fine, but a co-assembly needs S5 too
    assert "project_description:P2" not in by["M1"]["evidence"]  # the study only names the taxon


LINEAGE_CFG = {
    "mgnify_lineages": ["root:Environmental:Aquatic:Freshwater:Lake", "root:Environmental:Aquatic:Freshwater:Lentic"],
    "mgnify_lineages_exclude": ["root:Environmental:Aquatic:Freshwater:Lentic:Sediment"],
    "lineage_equivalents": {
        "root:Environmental:Aquatic:Freshwater:Reservoir": "root:Environmental:Aquatic:Freshwater:Lake"
    },
}


@pytest.mark.parametrize(
    "lineage, claim",
    [
        ("root:Environmental:Aquatic:Freshwater:Lake", ("in", 1.0)),
        ("root:Environmental:Aquatic:Freshwater:Lentic:Water", ("in", 1.0)),  # a sub-lineage
        ("root:Environmental:Aquatic:Freshwater:Lentic:Sediment", ("out", 1.0)),  # an excluded subtree
        ("root:Environmental:Aquatic:Freshwater", ("in", 0.8)),  # an ancestor: 4 of the target's 5 levels
        ("root:Environmental:Aquatic:Freshwater:Wetlands:Bog", ("out", 1.0)),
        ("root:Environmental:Aquatic:Freshwater:Reservoir", ("in", 1.0)),  # through lineage_equivalents
    ],
)
def test_lineage_claim(lineage, claim):
    assert score.lineage_claim(lineage, LINEAGE_CFG) == claim


def test_mixed_mgnify_studies_say_nothing():
    assert score.third_party_items("mgnify:MGYS1", "root:Mixed", LINEAGE_CFG) == []


def test_virtual_co_assembly_samples_are_replaced_by_their_members(monkeypatch):
    virtual = meta("SAMEA9", "mock sample for coassembly Erken sample of samples: ERS1,ERS2")
    members = [meta("SAMEA1", secondary_sample_accession="ERS1"), meta("SAMEA2", secondary_sample_accession="ERS2")]
    monkeypatch.setattr(
        ena, "sample_metadata", lambda ids: Response(200, [m for m in members if set(m.ids) & set(ids)])
    )
    rows = [row("M1", "SAMEA9", data_type="mag"), row("R1", "SAMEA5")]
    assert search.expand_co_assemblies(rows, [virtual]) == members
    assert rows[0]["sample_accession"] == "SAMEA1;SAMEA2" and rows[1]["sample_accession"] == "SAMEA5"


@pytest.mark.parametrize(
    "rec, evidence",
    [
        (
            {"data_type": "mag", "assembly_type": "Metagenome-Assembled Genome (MAG)"},
            "ena.assembly_type=metagenome-assembled genome (mag)",
        ),
        (
            {"data_type": "contigs", "assembly_type": "", "scientific_name": "lake water metagenome"},
            "ena.assembly_type=(empty); ena.scientific_name=lake water metagenome",
        ),
        ({"data_type": "reads"}, "ena.library_source=METAGENOMIC; ena.library_strategy=WGS or OTHER"),
    ],
)
def test_data_type_evidence(rec, evidence):
    assert search.data_type_evidence(rec) == evidence


GOLD_CFG = {
    "gold_ecosystems": ["Environmental > Aquatic > Freshwater"],
    "gold_ecosystems_exclude": ["Environmental > Aquatic > Freshwater > Lake > Sediment"],
}


@pytest.mark.parametrize(
    "path, claim",
    [
        ("Environmental > Aquatic > Freshwater > Lake > Unclassified", ("in", 1.0)),
        (
            "Environmental > Aquatic > Freshwater > Unclassified > Unclassified",
            ("in", 1.0),
        ),  # Unclassified says nothing
        ("Environmental > Aquatic > Freshwater > Lake > Sediment", ("out", 1.0)),
        ("Environmental > Aquatic > Unclassified > Unclassified > Unclassified", ("in", 0.67)),  # an ancestor
        ("Environmental > Aquatic > Marine > Coastal > Unclassified", ("out", 1.0)),
    ],
)
def test_gold_claim(path, claim):
    got = score.gold_claim(path, GOLD_CFG)
    assert (got[0], round(got[1], 2)) == claim


def test_gold_export_maps_biosamples_to_ecosystems(tmp_path, monkeypatch):
    import openpyxl
    from biome_finder.clients import gold, http

    wb = openpyxl.Workbook()
    levels = ["ECOSYSTEM", "ECOSYSTEM CATEGORY", "ECOSYSTEM TYPE", "ECOSYSTEM SUBTYPE", "SPECIFIC ECOSYSTEM"]
    for kind, gid, path in [
        ("Biosample", "Gb1", ["Environmental", "Aquatic", "Freshwater", "Lake", None]),
        ("Organism", "Go1", ["Environmental", "Aquatic", "Marine", "Coastal", "Unclassified"]),
    ]:
        ws = wb.create_sheet(kind)
        ws.append([f"{kind.upper()} GOLD ID", *(f"{kind.upper()} {lv}" for lv in levels)])
        ws.append([gid, *path])
    ws = wb.create_sheet("Sequencing Project")
    ws.append(["NCBI BIOSAMPLE ACCESSION", "BIOSAMPLE GOLD ID", "ORGANISM GOLD ID"])
    ws.append(["SAMN1", "Gb1", None])
    ws.append(["SAMN2", None, "Go1"])  # an isolate
    ws.append([None, "Gb1", None])  # no BioSample: skipped
    wb.save(tmp_path / "gold.xlsx")
    monkeypatch.setattr(gold, "CACHE_DIR", tmp_path / "cache")
    for _ in range(2):  # the second call reads the cache
        found = {g.biosample_accession: g for g in gold.biosamples(tmp_path / "gold.xlsx").data}
        assert found["SAMN1"].ecosystem_path == "Environmental > Aquatic > Freshwater > Lake > "
        assert found["SAMN2"].gold_id == "Go1" and len(found) == 2
