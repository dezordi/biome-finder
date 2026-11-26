# search

`search` builds `samples.tsv` for one biome: one row per public INSDC data record (genome,
metagenome assembly or read run) that could go into its catalogue. The curator decides the
rows, and `prepare` turns the accepted ones into workflow inputs.

```{note}
In development (milestone M2). This version builds the table from the survey records and
merges it with an existing `samples.tsv`. Taxonomy resolution, scoring, tiers and isolate
genomes come next; their columns are already in the table, empty.
```

## Usage

```bash
biome-finder survey --biomes biomes/<name>.yaml --out survey/
biome-finder search --biome biomes/<name>.yaml --survey survey/ --out samples.tsv
```

| Option | Description |
|---|---|
| `--biome` | The biome YAML config. |
| `--survey` | Survey folder holding `biomes/<name>/records.tsv` (default `survey/`). Run [survey](survey.md) on the biome first. |
| `--out` | The table (default `samples.tsv`). If it exists, the new search is merged into it. |
| `--since` | Only records public on or after this date (`YYYY-MM-DD`), for catalogue updates. |

## Records

`search` starts from the biome's survey records: MAGs, MGnify bins, metagenome assemblies and
read runs. Assemblies and runs whose sample has no evidence for the config's kingdom
(`kingdom_match` false in the survey) are left out. Lower-level records of a sample that
already has more processed data stay in the table; `prepare` marks them `not_selected`.

A genome's biome is the biome of the samples it was assembled from, so its `sample_accession`
holds its source samples (`derived from`, `;`-joined for co-assemblies) as primary IDs
(SAMEA/SAMN/SAMD), the same IDs its sample's assemblies and runs carry.

## Re-running

`search` merges by `record_accession`. It rewrites the automatic columns and keeps
`decision`, `curation_note`, `destination`, `local_path` and `prepare_status` from the
existing row. Rows the new search no longer finds are kept unchanged, so no curation is lost.

## Columns

The column order is fixed (PROPOSAL.md §3.6); columns are added, never renamed. Filled in
this version:

| Column | Content |
|---|---|
| `record_accession` | MAG (ENA WGS set), ERZ assembly or bin, or run accession |
| `data_type` | `mag`, `bin`, `contigs` or `reads`, from the survey search that found the record |
| `source_db` | `ENA` |
| `sample_accession` | The record's sample; for genomes, their source samples |
| `study_accession` | The record's ENA study (PRJ*) |
| `reported_tax_id`, `reported_tax_name` | As shown in the ENA record |
| `flags` | `in_existing_catalogue` (already in the MGnify catalogue), and for runs `long_read`, `single_end`, `low_reads` (fewer than 200,000 reads) |
| `read_stats` | Runs: read count, layout, platform |
| `first_public` | Date the record became public |
| `curation_note`, `decision` | For the curator: `decision` is `keep` or `skip` and overrides the tier |
| `destination`, `local_path`, `prepare_status` | Written by `prepare` |
