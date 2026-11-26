import csv

import pytest

from biome_finder import constants

from biome_finder import score, search
from biome_finder.clients import ena, ncbi, ols
from biome_finder.models import (
    NCBI_RANKS,
    from_row,
    EnaSampleMetadata,
    NcbiTaxdump,
    NcbiTaxon,
    OlsTerm,
    Response,
    SampleRow,
    SurveyRecord,
)


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
    monkeypatch.setattr(ena, "run_libraries", lambda ids: Response(200, []))
    monkeypatch.setattr(ena, "genomes", lambda query: Response(200, []))
    monkeypatch.setattr(ncbi, "assembly_summary", lambda accessions: Response(200, {}))
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
    assert list(next(iter(rows.values()))) == constants.SAMPLES_COLUMNS
    assert rows["R3"]["flags"] == "single_end;low_reads;sample_low_reads"
    assert rows["R1"]["reported_tax_id"] == "1" and rows["R1"]["tax_lineage"].endswith("s__freshwater metagenome")
    assert rows["M1"]["sample_accession"] == "SA;SB" and rows["M1"]["flags"] == "in_existing_catalogue"

    curated = {**rows["R1"], "decision": "skip", "curation_note": "control"}
    gone = {**rows["R3"], "record_accession": "R9"}
    search.write_tsv(out, constants.SAMPLES_COLUMNS, [SampleRow.from_tsv(curated), SampleRow.from_tsv(gone)])
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
        SampleRow(record_accession=a, data_type=t, reported_tax_id=tax)
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
    by = {r.record_accession: r for r in rows}
    assert by["M_merged"].domain == "prokaryotes"  # through the merged ID
    assert by["M_merged"].flags == [] and by["M_merged"].taxdump_date == "2026-11-01"
    assert by["M_merged"].tax_lineage == "d__Bacteria;k__;p__Nitrospirota;c__;o__;f__;g__Nitrospira;s__"
    assert by["M_euk"].flags == ["out_of_scope_taxon"]  # the other kingdom
    assert by["M_virus"].flags == ["out_of_scope_taxon"]
    assert by["M_deleted"].flags == ["taxid_deleted"] and by["M_deleted"].tax_lineage == ""  # a curator decides
    assert by["M_unknown"].flags == ["taxid_unknown"]
    assert by["R1"].flags == [] and by["R1"].domain == ""  # no scope for runs
    assert by["R1"].tax_lineage == "d__;k__;p__;c__;o__;f__;g__;s__freshwater metagenome"


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
        ("Overlying water from a 30-day control microcosm without apple snails.", ["manipulated"]),  # for curation
        ("mock sample for coassembly ErkenSummer", []),  # ENA's virtual co-assembly sample
        ("mock community DNA standard", ["mock"]),
        ("Synthetic microbial communities from Ohio, USA", ["synthetic"]),
        ("shotgun sequenced culture enrichment of surface water microbiome", ["enrichment"]),
        ("Enriched microbiome from surface water small river", ["enrichment"]),
        ("Walnut Creek spiked with Pseudomonas putida & Escherichia coli", ["spiked"]),
        ("Lake Erie surface water", []),
        ("Sorted cell/s from McNutts Creek water", ["single_cell"]),
        ("Lake Biwa virion metagenome (2016-2017), October hypolimnion", ["virome"]),
        ("Epilimnion water (VLP virome) shore of Caldonazzo Lake", ["virome"]),
        ("Freshwater microbial communities from Lake Mendota - viral genome information", []),
        ("non-axenic culture inoculated with surface water from freshwater lake", ["culture"]),
        ("non-axenic enrichment culture of lake water", ["enrichment"]),  # an enrichment is a culture
        ("lake water from an aquaculture region", []),
        ("Metagenome sample from freshwater metagenome Culture independent microbiome Potomac River", []),
    ],
)
def test_hard_filters(text, flags):
    assert score.hard_flags(meta("S1", text)) == flags


def test_sample_evidence_fills_reported_biome_and_sends_hard_filters_to_tier_d():
    rows = [
        SampleRow("M1", "mag", sample_accession=["SA", "SB"]),
        SampleRow("R1", "reads", sample_accession=["SA"], flags=["low_reads"]),
        SampleRow("R2", "reads", sample_accession=["SC"]),
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
    assert m1.reported_biome == "env_broad_scale=freshwater lake biome | river biome; organism=freshwater metagenome"
    assert m1.flags == ["spiked"] and m1.tier == "D"  # one source sample is enough
    assert r1.flags == ["low_reads"] and r1.tier == ""
    assert r2.reported_biome == "env_medium=water"


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
    envo = score.Evidence("sample_metadata", "env_local_scale", "in", 1.0, constants.W_ENVO, "")
    name = score.Evidence("metagenome_name", "organism", "in", 1.0, constants.W_METAGENOME, "")
    marine = score.Evidence("sample_metadata", "text", "out", 1.0, constants.W_SAMPLE_EXCLUDE, "marine")
    assert score.combine([envo, name]) == (86, False) and score.tier(86, False) == "A"  # 1 − 0.2·0.7
    assert score.combine([name]) == (30, False) and score.tier(30, False) == "C"  # a default name alone
    assert score.combine([name], siblings=1.0) == (40, False)
    value, conflict = score.combine([envo, name, marine])
    assert conflict and score.tier(value, conflict) == "C"


def run_biome(metadata, rows, texts=None, gold=None):
    search.biome_evidence(rows, {**CFG, **GOLD_CFG}, metadata, texts or {}, LAKE, SEDIMENT, gold_samples=gold)
    return {r.record_accession: r for r in rows}


def row(acc, sample, study="P1", data_type="reads"):
    return SampleRow(acc, data_type, sample_accession=sample.split(";"), study_accession=study)


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
    assert by["R1"].tier == "A" and by["R1"].curated_biome == CFG["catalogue_biome"]
    assert "sample_metadata:env_local_scale:in:freshwater lake" in by["R1"].evidence
    assert not any(e.startswith("sample_metadata:text") for e in by["R2"].evidence)  # a default title repeats the taxon
    assert not any("out:mock" in e for e in by["R3"].evidence)  # ENA's co-assembly sample, no mock community
    assert by["R4"].tier == "D" and "excluded_envo" in by["R4"].flags
    assert by["R5"].tier == "C" and by["R5"].needs_curation  # a conflict
    assert by["R5"].flags == ["off_biome_keyword"]  # the reason, visible
    assert by["M1"].tier == "C" and by["M1"].needs_curation  # S1 is fine, but a co-assembly needs S5 too
    assert not any(e.startswith("project_description:P2") for e in by["M1"].evidence)  # it only names the taxon


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
    assert rows[0].sample_accession == ["SAMEA1", "SAMEA2"] and rows[1].sample_accession == ["SAMEA5"]


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
    assert search.data_type_evidence(from_row(SurveyRecord, rec)) == evidence


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
    from biome_finder.clients import gold

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


def test_enrichment_goes_to_curation_but_spiked_samples_are_excluded():
    metadata = [
        meta("S1", "culture enrichment of lake water", local_environmental_context="freshwater lake"),
        meta("S2", "lake water spiked with E. coli", local_environmental_context="freshwater lake"),
    ]
    rows = [row("R1", "S1"), row("R2", "S2")]
    search.sample_evidence(rows, metadata)
    by = run_biome(metadata, rows)
    assert by["R1"].flags == ["enrichment"] and by["R1"].tier == "C" and by["R1"].needs_curation
    assert by["R2"].flags == ["spiked"] and by["R2"].tier == "D"


def test_a_sample_whose_own_evidence_points_only_elsewhere_is_off_biome():
    metadata = [
        meta("S1", "river sediment core", scientific_name="freshwater sediment metagenome"),  # SAMN07376294-like
        meta(
            "S2",
            "river sediment core",
            local_environmental_context="freshwater lake",
            scientific_name="freshwater metagenome",
        ),
    ]
    rows = [row("M1", "S1", data_type="mag"), row("M2", "S2", data_type="mag")]
    from biome_finder.models import GoldBiosample

    gold = {"S1": GoldBiosample("S1", "Gb1", "Environmental > Aquatic > Freshwater > Lake > Sediment")}
    by = run_biome(metadata, rows, {"P1": "Freshwater rivers"}, gold)  # the study's own text doesn't rescue S1
    assert by["M1"].tier == "D" and "off_biome" in by["M1"].flags
    alone = run_biome(metadata[:1], [row("M3", "S1", data_type="mag")], {"P1": "Freshwater rivers"})
    assert alone["M3"].tier == "C"  # without GOLD only free text says sediment: a curator decides
    assert by["M2"].tier == "C" and "off_biome" not in by["M2"].flags  # mixed: a curator decides


def test_run_flags_catch_amplicon_libraries_and_samples_without_a_usable_run():
    from biome_finder.models import EnaRunLibrary

    rows = [
        SampleRow(a, "reads", sample_accession=[s], read_count=n)
        for a, s, n in [("R1", "S1", 5_000_000), ("R2", "S2", 90_000), ("R3", "S3", 150_000), ("R4", "S3", 4_000_000)]
    ]
    libs = [
        EnaRunLibrary("R1", "OTHER", "PCR"),  # ERR1669025-like: amplicon under strategy OTHER
        EnaRunLibrary("R2", "OTHER", "other"),
        EnaRunLibrary("R3", "WGS", "RANDOM PCR"),  # PCR during library prep is not amplicon
    ]
    search.run_flags(rows, libs)
    flags = {r.record_accession: r.flags for r in rows}
    assert flags == {"R1": ["amplicon_only"], "R2": ["sample_low_reads"], "R3": [], "R4": []}  # R4's sample: a big run
    assert rows[2].data_type_evidence.endswith("ena.library_selection=RANDOM PCR")


def test_a_generic_keyword_does_not_rescue_a_sample_gold_files_as_sediment():
    sed = score.Evidence("third_party_biome", "gold:Gb1", "out", 1.0, 0.9, "… > Lake > Sediment")
    word = score.Evidence("sample_metadata", "text", "in", 1.0, constants.W_SAMPLE_KEYWORD, "freshwater")
    envo = score.Evidence("sample_metadata", "env_medium", "in", 1.0, constants.W_ENVO, "lake water")
    assert score.off_biome([sed, word])  # SRR7514475: alias …_Sed, isolation_source "freshwater"
    assert not score.off_biome([sed, envo])  # GOLD against an ENVO term: a conflict for a curator
    assert not score.off_biome([word])


def test_a_virome_study_flags_its_rows():
    rows = [row("R1", "S1", study="PRJNA672207"), row("R2", "S2", study="P2")]
    search.study_flags(
        rows, {"PRJNA672207": "To investigate the freshwater virome from the Yangtze River", "P2": "lake water"}
    )
    assert [r.flags for r in rows] == [["virome"], []]


def genome(acc, sample, package="", strain="", env="", isolate="", culture=""):
    from biome_finder.models import EnaGenome

    return EnaGenome(
        acc, "", "contig", env, package, strain, culture, isolate, sample, "PRJ1", "2", "Bacteria", "2024-01-01", ""
    )


@pytest.mark.parametrize(
    "g, data_type",
    [
        (genome("G1", "S1", package="MIMAG.water.6.0"), "mag"),
        (genome("G2", "S2", package="MISAG.water.6.0"), "sag"),
        (genome("G3", "S3", env="true"), "mag"),
        (genome("G4", "S4", package="MIGS.ba.water.6.0", strain="LW1"), "isolate"),
        (
            genome("G5", "S5", package="Microbe, viral or environmental", isolate="UBA10892"),
            "unknown",
        ),  # an NCBI MAG, untold
    ],
)
def test_genome_type(g, data_type):
    assert search.genome_type(g)[0] == data_type


def test_genome_rows_skip_survey_genomes_and_flag_type_strains(monkeypatch):
    monkeypatch.setattr(ena, "derived_from", lambda ids: Response(200, {"S1": ["ERS9"]}))
    from biome_finder.models import EnaSample

    monkeypatch.setattr(ena, "sample_hosts", lambda ids: Response(200, [EnaSample("SAMEA9", "ERS9", "", "")]))
    found = [
        genome("G1", "S1", package="MIMAG.water.6.0"),
        genome("G2", "SURVEY", package="MIMAG.water.6.0"),  # the survey's MAG, as a GCA
        genome("G3", "S3", strain="DSM 123, type strain"),
    ]
    rows = {r.record_accession: r for r in search.genome_rows(found, {"SURVEY"})}
    assert set(rows) == {"G1", "G3"}
    assert rows["G1"].sample_accession == ["SAMEA9"]  # a MAG's source sample, as the primary ID
    assert rows["G3"].data_type == "isolate" and rows["G3"].flags == ["type_strain"]


def test_genome_query_splits_phrases_into_words():
    q = search.genome_query({"kingdom": "prokaryotes", "keywords_include": ["lake water"]})
    assert '(isolation_source="*lake*" AND isolation_source="*water*")' in q
    assert 'NOT assembly_type="metagenome-assembled genome (mag)"' in q


def test_a_keyword_alone_never_makes_a_sample_off_biome():
    # SAMEA114327776: "cell culturing; microcosm", only weak evidence for freshwater
    micro = score.Evidence("sample_metadata", "text", "out", 1.0, constants.W_SAMPLE_EXCLUDE, "microcosm")
    word = score.Evidence("sample_metadata", "text", "in", 1.0, constants.W_SAMPLE_KEYWORD, "freshwater")
    assert not score.off_biome([micro, word])


def test_a_mesocosm_study_sends_its_rows_to_curation():
    rows = [row("M1", "S1", study="PRJEB56877")]
    search.study_flags(rows, {"PRJEB56877": "metagenomic and 16S rRNA data from a mesocosm experiment at Jiricka Pond"})
    assert rows[0].flags == ["manipulated"]


def test_sample_row_fields_are_the_table_columns_and_survive_a_round_trip():
    from dataclasses import fields

    columns = [f.name for f in fields(SampleRow) if f.metadata.get("column", f.name) is not None]
    assert columns == constants.SAMPLES_COLUMNS
    r = SampleRow("M1", "mag", sample_accession=["S1", "S2"], flags=["a", "b"], evidence=["x:y", "z"], score=80)
    r.needs_curation = False
    assert SampleRow.from_tsv(r.to_tsv()) == SampleRow(**{**r.__dict__, "read_count": 0})
    assert r.to_tsv()["sample_accession"] == "S1;S2" and r.to_tsv()["evidence"] == "x:y; z"
