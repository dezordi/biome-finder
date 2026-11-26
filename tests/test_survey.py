from biome_finder.clients import ena, mgnify
from biome_finder.models import (
    EnaLinkedWgs,
    EnaBin,
    EnaMag,
    EnaRun,
    EnaSample,
    EnaStudy,
    MgnifyCatalogue,
    MgnifyCatalogueGenome,
    Response,
    from_row,
)
from biome_finder.survey import summarize


def run(acc, study, reads, layout="PAIRED", platform="ILLUMINA", date="2024-01-01"):
    return {
        "run_accession": acc,
        "study_accession": study,
        "read_count": reads,
        "base_count": reads * 300,
        "first_public": date,
        "library_layout": layout,
        "instrument_platform": platform,
    }


def sample(primary, secondary="", host="", host_tax_id=""):
    return EnaSample(primary, secondary, host, host_tax_id)


def ok(data):
    return lambda *args, **kwargs: Response(200, data)


def test_from_row_converts_numbers_and_empties():
    r = from_row(EnaRun, {"run_accession": "R1", "read_count": "5", "base_count": "", "lat": "10.5;10.6"})
    assert r.read_count == 5 and r.base_count is None and r.lat == "10.5;10.6" and r.lon == ""
    assert r.country == ""  # missing column
    g = from_row(
        MgnifyCatalogueGenome,
        {"Genome": "MGYG1", "Genome_type": "MAG", "Genome_accession": "ERZ1", "Sample_accession": "SAMEA1"},
    )
    assert g.genome_type == "MAG" and g.sample_accession == "SAMEA1"


def test_summarize_counts_and_cutoff():
    runs = [
        run("R1", "S1", 5_000_000),  # big, paired, old
        run("R2", "S1", 100_000),  # below 1M: ignored
        run("R3", "S2", 2_000_000, "SINGLE", date="2026-01-01"),  # big, single, new
        run("R4", "S3", 1_000_000, "SINGLE", "OXFORD_NANOPORE", "2026-02-01"),  # big, long, new
        {**run("R5", "S3", 0), "read_count": None, "base_count": None},  # missing counts
    ]
    mags = [
        {"study_accession": "M1", "first_public": "2025-01-01"},
        {"study_accession": "M1", "first_public": "2026-01-01"},
    ]
    s = summarize(mags, [], runs, "2025-07-16")
    assert s["runs"] == 5 and s["runs_1m"] == 3
    assert s["runs_1m_new"] == 2
    assert s["runs_1m_paired"] == 1 and s["runs_1m_long"] == 1
    assert s["read_studies_1m"] == 3
    assert s["gbases_1m"] == round(8_000_000 * 300 / 1e9, 1)
    assert s["mags"] == 2 and s["mags_new"] == 1 and s["mag_studies"] == 1


def test_summarize_without_catalogue_leaves_new_empty():
    s = summarize([], [], [run("R1", "S1", 2_000_000)], "")
    assert s["runs_1m_new"] == "" and s["mags_new"] == ""


def test_host_mags_needs_every_source_sample_on_the_host(monkeypatch):
    from biome_finder import survey

    origin = {"M_ok": ["SAMEA1"], "M_coassembly": ["SAMEA1", "ERS2"], "M_keyword": ["SAMEA3"], "M_none": []}
    hosts = [
        sample("SAMEA1", host_tax_id="4565"),
        sample("SAMEA2", "ERS2", "Zea mays", "4577"),
        sample("SAMEA3", host="winter WHEAT"),
    ]
    monkeypatch.setattr(ena, "derived_from", ok(origin))
    monkeypatch.setattr(ena, "sample_hosts", ok(hosts))
    mags = [{"sample_accession": a} for a in origin]
    kept = survey.host_mags(mags, {"name": "wheat", "host_tax_id": 4565, "host_keywords": ["wheat"]})
    assert [m["sample_accession"] for m in kept] == ["M_ok", "M_keyword"]
    assert mags[1]["source_samples"] == "SAMEA1;ERS2" and mags[1]["host"] == "Zea mays"


def test_parse_tsv_rejects_truncated_page():
    from biome_finder.clients.ena import parse_tsv

    ok = "accession\tfirst_public\nA1\t2024-01-01\n"
    assert parse_tsv(ok, ["accession", "first_public"]) == [{"accession": "A1", "first_public": "2024-01-01"}]
    assert parse_tsv(ok + "Premature EOF\n", ["accession", "first_public"]) is None


def test_host_filter():
    from biome_finder.survey import host_filter

    assert host_filter({"host": "none"}) == ""
    assert (
        host_filter({"host_tax_id": 4565, "host_keywords": ["wheat", "Triticum"]})
        == ' AND (host_tax_id=4565 OR host="*wheat*" OR host="*Triticum*")'
    )


def test_taxa_filter_defaults_plus_extra():
    from biome_finder.survey import taxa_filter

    clause = taxa_filter({"exclude_taxa_extra": [2763]})
    for t in (10239, 33208, 58023, 2870, 2763):
        assert f"NOT tax_tree({t})" in clause


def test_studies_count_records_and_samples_per_route():
    from biome_finder.survey import studies, tag

    runs = tag(
        [
            from_row(EnaRun, run("R1", "S1", 5_000_000)),
            from_row(EnaRun, run("R2", "S1", 2_000_000)),
            from_row(EnaRun, run("R3", "S2", 3_000_000)),
            from_row(EnaRun, run("R4", "S4", 4_000_000)),
        ],
        "reads",
        "run_accession",
        "soil metagenome",
        "",
    )
    for r, s, route in zip(runs, ["SA1", "SA1", "SA2", "SA1"], ["assemble", "assemble", "covered", "assemble"]):
        r["sample_accession"], r["route"] = s, route
    mags = tag(
        [
            from_row(EnaMag, {"accession": "M1", "study_accession": "S1", "first_public": "2024-01-01"}),
            from_row(EnaMag, {"accession": "M2", "study_accession": "S1", "first_public": "2024-01-01"}),
        ],
        "mag",
        "accession",
        "soil metagenome",
        "",
    )
    mags[1]["in_catalogue"] = True
    contig = {
        "data_type": "contigs",
        "study_accession": "S3",
        "sample_accession": "SA3",
        "metagenome_taxon": "soil metagenome",
    }
    contig |= {"first_public": "2023-01-01", "route": "covered", "covered_by": "S1"}
    rows = {r["study_accession"]: r for r in studies(runs + mags + [contig])}
    assert rows["S1"]["runs_to_assemble"] == 2 and rows["S1"]["samples_to_assemble"] == 1
    assert (
        rows["S4"]["runs_to_assemble"] == 1 and rows["S4"]["samples_to_assemble"] == 0
    )  # SA1 counts in S1, the larger
    assert rows["S1"]["mags_to_add"] == 1  # M2 is already in the catalogue
    assert "S2" not in rows  # its only run is covered, so it has nothing on a route
    assert rows["S3"]["assemblies_with_mags"] == 1 and rows["S3"]["mags_studies"] == "S1"


def test_mgnify_bins_need_every_source_sample_in_the_biome(monkeypatch):
    from biome_finder import survey

    records = [
        {
            "sample_accession": "SAMEA1",
            "metagenome_taxon": "marine metagenome",
            "study_accession": "PRJEB_EMG",
            "from_mgnify": True,
        },
        {
            "sample_accession": "SAMEA2",
            "metagenome_taxon": "seawater metagenome",
            "study_accession": "PRJEB_X",
            "from_mgnify": False,
        },
    ]
    bins = [
        from_row(
            EnaBin,
            {
                "analysis_accession": a,
                "study_accession": "PRJEB_EMG",
                "sample_accession": s,
                "first_public": "2025-01-01",
            },
        )
        for a, s in [("B_in", "SB1"), ("B_mixed", "SB2"), ("B_out", "SB3")]
    ]
    monkeypatch.setattr(ena, "search", ok(bins))
    monkeypatch.setattr(ena, "derived_from", ok({"SB1": ["ERS1"], "SB2": ["SAMEA1", "SAMEA9"], "SB3": []}))
    monkeypatch.setattr(ena, "sample_hosts", ok([sample("SAMEA1", "ERS1")]))
    kept = survey.mgnify_bins({"name": "t"}, {2759: "Eukaryota"}, records, "")
    assert [b["accession"] for b in kept] == ["B_in"]
    assert kept[0]["metagenome_taxon"] == "marine metagenome" and kept[0]["from_mgnify"] is True


def test_study_kingdom_reads_title_and_description():
    from biome_finder.survey import study_kingdom

    assert study_kingdom("Size-fractionated protist and diatom communities", "eukaryotes") == "study_text"
    assert study_kingdom("Bacterial and archaeal 16S diversity", "eukaryotes") == "other_kingdom"
    assert study_kingdom("Bacteria and fungi in soil", "eukaryotes") == "study_text"  # both named
    assert study_kingdom("Microbial community of lake water", "eukaryotes") == "none"
    assert study_kingdom("Bacterial communities", "prokaryotes") == "study_text"
    assert study_kingdom("algal bloom, not 'falgarone'", "eukaryotes") == "study_text"


def test_kingdom_match_rules(monkeypatch):
    from biome_finder import survey

    texts = [
        EnaStudy("S_euk", "picoeukaryote diversity", ""),
        EnaStudy("S_prok", "bacterial 16S", ""),
        EnaStudy("S_none", "microbial community", ""),
    ]
    monkeypatch.setattr(ena, "studies", ok(texts))
    monkeypatch.setattr(ena, "derived_from", ok({}))
    monkeypatch.setattr(ena, "sample_hosts", ok([]))

    def data():
        return [
            {"study_accession": s, "sample_accession": smp, "read_count": 5_000_000}
            for s, smp in [("S_euk", "A"), ("S_prok", "B"), ("S_none", "C"), ("S_prok", "SAMEA_G")]
        ]

    euk = data()
    survey.kingdom_evidence({"kingdom": "eukaryotes"}, [{"study_accession": "G", "source_samples": "SAMEA_G"}], [], euk)
    assert [r["kingdom_evidence"] for r in euk] == ["study_text", "other_kingdom", "none", "genome"]
    assert [r["kingdom_match"] for r in euk] == [True, False, False, True]
    prok = data()
    survey.kingdom_evidence({"kingdom": "prokaryotes"}, [], [], prok)
    assert [r["kingdom_match"] for r in prok] == [False, True, True, True]  # only euk-only studies drop


def test_kingdom_evidence_is_shared_by_a_samples_records(monkeypatch):
    from biome_finder import survey

    texts = [EnaStudy("S_euk", "Picoeukaryote and protist diversity", ""), EnaStudy("S_none", "lake water", "")]
    monkeypatch.setattr(ena, "studies", ok(texts))
    monkeypatch.setattr(ena, "derived_from", ok({}))
    monkeypatch.setattr(ena, "sample_hosts", ok([]))

    def data():
        return [{"study_accession": s, "sample_accession": "X", "read_count": 5_000_000} for s in ("S_euk", "S_none")]

    euk = data()
    survey.kingdom_evidence({"kingdom": "eukaryotes"}, [], [], euk)
    assert [r["kingdom_evidence"] for r in euk] == ["study_text", "study_text"]  # one verdict per sample
    assert euk[1]["kingdom_terms"] == "eukaryot;protist"
    prok = data()
    survey.kingdom_evidence({"kingdom": "prokaryotes"}, [], [], prok)
    assert [r["kingdom_evidence"] for r in prok] == ["none", "none"]  # the neutral study outweighs other_kingdom
    assert all(r["kingdom_match"] for r in prok)


def test_routes_each_sample_enters_at_its_most_processed_data(monkeypatch):
    from biome_finder import survey

    monkeypatch.setattr(ena, "studies", ok([]))
    monkeypatch.setattr(ena, "derived_from", ok({}))
    monkeypatch.setattr(ena, "sample_hosts", ok([sample("S_genome", "ERS1")]))
    genomes = [{"study_accession": "G", "sample_accession": "M1", "source_samples": "ERS1"}]
    contigs = [{"study_accession": "P", "sample_accession": s} for s in ("S_genome", "S_asm")]
    runs = [
        {"study_accession": "P", "sample_accession": s, "read_count": n}
        for s, n in [("S_genome", 5_000_000), ("S_asm", 5_000_000), ("S_reads", 5_000_000), ("S_small", 1000)]
    ]
    survey.kingdom_evidence({"kingdom": "prokaryotes"}, genomes, contigs, runs)
    assert [r["route"] for r in contigs] == ["covered", "bin"]
    assert [r["route"] for r in runs] == ["covered", "covered", "assemble", ""]
    assert [r["covered_by"] for r in contigs] == ["G", ""]  # the study holding the genome
    assert genomes[0]["source_samples"] == "S_genome"  # derived from's secondary ID, as the primary one


def test_catalogue_membership_matches_accession_or_sample(monkeypatch):
    from biome_finder import survey

    listed = [
        MgnifyCatalogueGenome("MGYG1", "MAG", "ERZ1", "SAMEA9"),
        MgnifyCatalogueGenome("MGYG2", "Isolate", "GCA_1", "SAMEA8"),
    ]
    monkeypatch.setattr(mgnify, "catalogue_genomes", ok(listed))
    cat = MgnifyCatalogue("c", "prokaryotes", "2024-01-01", 1, 2, "https://ftp/c")
    genomes = [
        {"accession": "ERZ1", "sample_accession": "X"},
        {"accession": "CAB01", "sample_accession": "SAMEA9"},
        {"accession": "ERZ2", "sample_accession": "Y"},
    ]
    assert survey.catalogue_membership(cat, genomes) == 1
    assert [g["in_catalogue"] for g in genomes] == [True, True, False]
    assert survey.catalogue_membership(None, genomes) == "" and not any(g["in_catalogue"] for g in genomes)


def test_ena_search_rejects_unknown_result_and_unfiltered_search():
    import pytest
    from biome_finder.models import EnaSearch

    EnaSearch("read_run", "tax_eq(408172)")
    EnaSearch("sample", accessions=("SAMEA1",))
    with pytest.raises(ValueError, match="unknown ENA result"):
        EnaSearch("read_runs", "tax_eq(408172)")
    with pytest.raises(ValueError, match="needs a query or accessions"):
        EnaSearch("read_run")


def test_linked_assemblies_cover_samples_survey_missed(monkeypatch):
    from biome_finder import survey

    def wgs(acc, sample, kind, organism):
        return EnaLinkedWgs(acc, "P2", sample, "2024-01-01", kind, "410658", organism, "", "")

    found = [
        wgs("NCBI_ASM", "S_asm", "", "soil metagenome"),  # an NCBI metagenome assembly as a WGS set
        wgs("MAG_X", "S_mag", "Metagenome-Assembled Genome (MAG)", "Nitrospira sp."),
        wgs("MAG_Y", "S_bin", "Metagenome-Assembled Genome (MAG)", "Nitrospira sp."),
        wgs("MITO", "S_free", "binned metagenome", "freshwater metagenome"),  # kingdom unknown: covers nothing
    ]
    seen = {}

    def fake(samples, clause=""):
        seen["samples"], seen["clause"] = set(samples), clause
        return Response(200, found)

    monkeypatch.setattr(ena, "sample_assemblies", fake)
    evidence = {"kingdom_evidence": "none", "kingdom_terms": "", "kingdom_match": True, "covered_by": ""}
    runs = [
        {"accession": f"R_{s}", "sample_accession": s, "metagenome_taxon": "soil metagenome", "route": "assemble"}
        | evidence
        for s in ("S_asm", "S_mag", "S_free")
    ]
    contigs = [{"accession": "ERZ1", "sample_accession": "S_bin", "route": "bin"} | evidence]
    new = survey.linked_assemblies({"kingdom": "prokaryotes"}, contigs, runs, set(), "")
    assert seen["samples"] == {"S_asm", "S_mag", "S_free", "S_bin"}
    assert seen["clause"].startswith(" AND (tax_tree(2) OR tax_tree(2157) OR tax_tree(408169))")
    assert [(r["accession"], r["route"], r["data_type"]) for r in new] == [("NCBI_ASM", "bin", "contigs")]
    assert [r["route"] for r in runs] == ["covered", "covered", "assemble"]
    assert contigs[0]["route"] == "covered" and contigs[0]["covered_by"] == "P2"  # a genome covers it
