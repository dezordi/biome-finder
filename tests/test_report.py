"""Smoke test: the dashboard renders from a minimal survey folder."""

from biome_finder import report, survey


def test_report_writes_dashboard(tmp_path):
    recs = [
        {"data_type": "mag", "accession": "M1", "study_accession": "S1", "first_public": "2020-01-01",
         "scope_taxon": "Bacteria", "metagenome_taxon": "soil metagenome", "from_mgnify": "False",
         "in_catalogue": "False"},
        {"data_type": "bin", "accession": "B1", "study_accession": "S3", "first_public": "2024-01-01",
         "scope_taxon": "Bacteria", "metagenome_taxon": "soil metagenome", "from_mgnify": "True",
         "in_catalogue": "True"},
        {"data_type": "reads", "accession": "R1", "study_accession": "S1", "sample_accession": "SA1",
         "first_public": "2023-05-01", "read_count": "5000000", "base_count": "1500000000",
         "instrument_platform": "ILLUMINA", "library_layout": "PAIRED", "lat": "10.2", "lon": "-50.1",
         "metagenome_taxon": "soil metagenome", "kingdom_match": "True", "route": "assemble"},
        {"data_type": "contigs", "accession": "E1", "study_accession": "S2", "sample_accession": "SA2",
         "first_public": "2021-02-01", "metagenome_taxon": "soil metagenome", "kingdom_match": "True",
         "route": "bin", "from_mgnify": "True"},
    ]
    for biome in ("soil-prok", "lake-euk"):
        d = tmp_path / "biomes" / biome
        d.mkdir(parents=True)
        survey.write_tsv(d / "records.tsv", survey.RECORD_COLUMNS, recs)
    counts = {"mags": 1, "mags_mgnify": 0, "mgnify_bins": 1, "mags_to_add": 1, "mgnify_bins_to_add": 0,
              "contigs_kingdom": 1, "runs_1m_kingdom": 1, "gbases_1m_kingdom": 1.5, "contigs_to_bin": 1,
              "runs_to_assemble": 1, "mgnify_samples": 3}
    rows = []
    for biome, extra in [("soil-prok", {"kingdom": "prokaryotes", "catalogue": "soil-v1-0", "catalogue_updated": "2022-01-01",
                                        "catalogue_species": 100, "catalogue_genomes": 120, "catalogue_mags": 90,
                                        "runs_1m_kingdom_new": 1}),
                         ("lake-euk", {"kingdom": "eukaryotes"})]:
        for taxon in ("soil metagenome", "TOTAL"):
            rows.append({**counts, **extra, "biome": biome, "metagenome_taxon": taxon})
    survey.write_tsv(tmp_path / "survey.tsv", survey.SURVEY_COLUMNS, rows)
    report.write(tmp_path)
    html = (tmp_path / "report.html").read_text()
    assert "All catalogues" in html and "One catalogue" in html and "soil-v1-0" in html
    assert "Genomes to add" in html
