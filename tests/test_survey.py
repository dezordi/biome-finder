from biome_finder.survey import summarize


def run(acc, study, reads, layout="PAIRED", platform="ILLUMINA", date="2024-01-01"):
    return {"run_accession": acc, "study_accession": study, "read_count": str(reads),
            "base_count": str(reads * 300), "first_public": date,
            "library_layout": layout, "instrument_platform": platform}


def test_summarize_counts_and_cutoff():
    runs = [
        run("R1", "S1", 5_000_000),                                   # big, paired, old
        run("R2", "S1", 100_000),                                     # below 1M: ignored
        run("R3", "S2", 2_000_000, "SINGLE", date="2026-01-01"),      # big, single, new
        run("R4", "S3", 1_000_000, "SINGLE", "OXFORD_NANOPORE", "2026-02-01"),  # big, long, new
        {**run("R5", "S3", 0), "read_count": "", "base_count": ""},  # missing counts
    ]
    mags = [{"study_accession": "M1", "first_public": "2025-01-01"},
            {"study_accession": "M1", "first_public": "2026-01-01"}]
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
    from biome_finder import sources, survey
    origin = {"M_ok": ["SAMEA1"], "M_coassembly": ["SAMEA1", "ERS2"], "M_keyword": ["SAMEA3"],
              "M_none": []}
    hosts = {"SAMEA1": {"host": "", "host_tax_id": "4565"},
             "ERS2": {"host": "Zea mays", "host_tax_id": "4577"},
             "SAMEA3": {"host": "winter WHEAT", "host_tax_id": ""}}
    monkeypatch.setattr(sources, "derived_from", lambda samples: origin)
    monkeypatch.setattr(sources, "sample_hosts", lambda samples: hosts)
    mags = [{"sample_accession": a} for a in origin]
    kept = survey.host_mags(mags, {"name": "wheat", "host_tax_id": 4565, "host_keywords": ["wheat"]})
    assert [m["sample_accession"] for m in kept] == ["M_ok", "M_keyword"]
    assert mags[1]["source_samples"] == "SAMEA1;ERS2" and mags[1]["host"] == "Zea mays"


def test_parse_tsv_rejects_truncated_page():
    from biome_finder.sources import parse_tsv
    ok = "accession\tfirst_public\nA1\t2024-01-01\n"
    assert parse_tsv(ok, ["accession", "first_public"]) == [{"accession": "A1", "first_public": "2024-01-01"}]
    assert parse_tsv(ok + "Premature EOF\n", ["accession", "first_public"]) is None


def test_host_filter():
    from biome_finder.survey import host_filter
    assert host_filter({"host": "none"}) == ""
    assert host_filter({"host_tax_id": 4565, "host_keywords": ["wheat", "Triticum"]}) == \
        ' AND (host_tax_id=4565 OR host="*wheat*" OR host="*Triticum*")'


def test_in_scope_subtree_minus_excludes():
    from biome_finder.survey import in_scope
    soil = ["root:Environmental:Terrestrial:Soil"]
    no = ["root:Environmental:Terrestrial:Soil:Permafrost", "root:Environmental:Terrestrial:Soil:Wetlands"]
    assert in_scope("root:Environmental:Terrestrial:Soil", soil, no)
    assert in_scope("root:Environmental:Terrestrial:Soil:Forest soil", soil, no)
    assert not in_scope("root:Environmental:Terrestrial:Soil:Wetlands:Permafrost", soil, no)
    assert not in_scope("root:Environmental:Terrestrial:Soil:Permafrost", soil, no)
    assert not in_scope("root:Environmental:Terrestrial:Soilless", soil, no)  # prefix, not a child
    assert not in_scope("root:Environmental:Terrestrial", soil, no)


def test_taxa_filter_defaults_plus_extra():
    from biome_finder.survey import taxa_filter
    clause = taxa_filter({"exclude_taxa_extra": [2763]})
    for t in (10239, 33208, 58023, 2870, 2763):
        assert f"NOT tax_tree({t})" in clause


def test_studies_groups_records_and_counts_samples():
    from biome_finder.survey import studies, tag
    runs = tag([run("R1", "S1", 5_000_000), run("R2", "S1", 2_000_000), run("R3", "S2", 3_000_000)],
               "reads", "run_accession", "soil metagenome", "")
    for r, s in zip(runs, ["SA1", "SA1", "SA2"]):
        r["sample_accession"] = s
    mags = tag([{"accession": "M1", "study_accession": "S1", "first_public": "2024-01-01"}],
               "mag", "accession", "soil metagenome", "")
    rows = {r["study_accession"]: r for r in studies(runs + mags, "")}
    assert rows["S1"]["runs"] == 2 and rows["S1"]["samples"] == 1 and rows["S1"]["mags"] == 1
    assert rows["S2"]["mags"] == 0
    assert list(rows)[0] == "S1"  # largest study first


def test_mgnify_bins_need_every_source_sample_in_the_biome(monkeypatch):
    from biome_finder import sources, survey
    records = [{"sample_accession": "SAMEA1", "metagenome_taxon": "marine metagenome",
                "study_accession": "PRJEB_EMG", "from_mgnify": True},
               {"sample_accession": "SAMEA2", "metagenome_taxon": "seawater metagenome",
                "study_accession": "PRJEB_X", "from_mgnify": False}]
    bins = [{"analysis_accession": a, "study_accession": "PRJEB_EMG", "sample_accession": s,
             "first_public": "2025-01-01"} for a, s in [("B_in", "SB1"), ("B_mixed", "SB2"), ("B_out", "SB3")]]
    monkeypatch.setattr(sources, "ena_search", lambda result, query, fields: bins)
    monkeypatch.setattr(sources, "derived_from", lambda samples: {"SB1": ["ERS1"], "SB2": ["SAMEA1", "SAMEA9"], "SB3": []})
    monkeypatch.setattr(sources, "sample_hosts", lambda ids: {"ERS1": {"sample_accession": "SAMEA1"}})
    kept = survey.mgnify_bins({"name": "t"}, {2759: "Eukaryota"}, records, "")
    assert [b["accession"] for b in kept] == ["B_in"]
    assert kept[0]["metagenome_taxon"] == "marine metagenome" and kept[0]["from_mgnify"] is True


def test_study_kingdom_reads_title_and_description():
    from biome_finder.survey import study_kingdom
    assert study_kingdom("Size-fractionated protist and diatom communities", "eukaryotes") == "study_text"
    assert study_kingdom("Bacterial and archaeal 16S diversity", "eukaryotes") == "other_kingdom"
    assert study_kingdom("Bacteria and fungi in soil", "eukaryotes") == "study_text"   # both named
    assert study_kingdom("Microbial community of lake water", "eukaryotes") == "none"
    assert study_kingdom("Bacterial communities", "prokaryotes") == "study_text"
    assert study_kingdom("algal bloom, not 'falgarone'", "eukaryotes") == "study_text"


def test_kingdom_match_rules(monkeypatch):
    from biome_finder import sources, survey
    texts = {"S_euk": "picoeukaryote diversity", "S_prok": "bacterial 16S", "S_none": "microbial community"}
    monkeypatch.setattr(sources, "study_texts", lambda accs: texts)
    monkeypatch.setattr(sources, "derived_from", lambda samples: {})
    monkeypatch.setattr(sources, "sample_hosts", lambda ids: {})

    def data():
        return [{"study_accession": s, "sample_accession": smp, "read_count": "5000000"} for s, smp in
                [("S_euk", "A"), ("S_prok", "B"), ("S_none", "C"), ("S_prok", "SAMEA_G")]]
    euk = data()
    survey.kingdom_evidence({"kingdom": "eukaryotes"}, [{"source_samples": "SAMEA_G"}], [], euk)
    assert [r["kingdom_evidence"] for r in euk] == ["study_text", "other_kingdom", "none", "genome"]
    assert [r["kingdom_match"] for r in euk] == [True, False, False, True]
    prok = data()
    survey.kingdom_evidence({"kingdom": "prokaryotes"}, [], [], prok)
    assert [r["kingdom_match"] for r in prok] == [False, True, True, True]  # only euk-only studies drop


def test_routes_each_sample_enters_at_its_most_processed_data(monkeypatch):
    from biome_finder import sources, survey
    monkeypatch.setattr(sources, "study_texts", lambda accs: {})
    monkeypatch.setattr(sources, "derived_from", lambda samples: {})
    monkeypatch.setattr(sources, "sample_hosts", lambda ids: {})
    genomes = [{"sample_accession": "M1", "source_samples": "S_genome"}]
    contigs = [{"study_accession": "P", "sample_accession": s} for s in ("S_genome", "S_asm")]
    runs = [{"study_accession": "P", "sample_accession": s, "read_count": n} for s, n in
            [("S_genome", "5000000"), ("S_asm", "5000000"), ("S_reads", "5000000"), ("S_small", "1000")]]
    survey.kingdom_evidence({"kingdom": "prokaryotes"}, genomes, contigs, runs)
    assert [r["route"] for r in contigs] == ["covered", "bin"]
    assert [r["route"] for r in runs] == ["covered", "covered", "assemble", ""]


def test_catalogue_membership_matches_accession_or_sample(monkeypatch):
    from biome_finder import sources, survey
    listed = [{"Genome_type": "MAG", "Genome_accession": "ERZ1", "Sample_accession": "SAMEA9"},
              {"Genome_type": "Isolate", "Genome_accession": "GCA_1", "Sample_accession": "SAMEA8"}]
    monkeypatch.setattr(sources, "mgnify_catalogue_genomes", lambda cat: listed)
    genomes = [{"accession": "ERZ1", "sample_accession": "X"}, {"accession": "CAB01", "sample_accession": "SAMEA9"},
               {"accession": "ERZ2", "sample_accession": "Y"}]
    assert survey.catalogue_membership({"catalogue_id": "c"}, genomes) == 1
    assert [g["in_catalogue"] for g in genomes] == [True, True, False]
    assert survey.catalogue_membership(None, genomes) == "" and not any(g["in_catalogue"] for g in genomes)
