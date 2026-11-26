# search

`search` builds `samples.tsv` for one biome: one row per public INSDC data record (genome,
metagenome assembly or read run) that could go into its catalogue. The curator decides the
rows, and `prepare` turns the accepted ones into workflow inputs.

```{note}
In development (milestone M2). This version builds the table from the survey records and
merges it with an existing `samples.tsv`, with the taxonomy resolved against NCBI, each
sample's reported biome, the hard filters (tier D) and the rule-based biome score (curation
actions 2–5, tiers A–C). Third-party biomes (MGnify, GOLD), coordinates and isolate genomes
come next; their columns are already in the table, empty.
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
| `-v` / `--verbose` | Before the command (`biome-finder -v search …`): also log every HTTP request, cached or not. |
| `-v` / `--verbose` | Before the command (`biome-finder -v search …`): also log every HTTP request, cached or not. |
| `--taxdump` | NCBI taxonomy release (`YYYY-MM-DD`): the monthly `new_taxdump_<YYYY-MM>-01.zip` of that month. Default: the latest. |

Progress goes to stderr, one timestamped line per step: records read, taxonomy flags, samples
found in ENA, hard-filter counts, ENVO terms expanded, tiers per data type, and the merge.

Progress goes to stderr, one timestamped line per step: records read, taxonomy flags, samples
found in ENA, hard-filter counts, ENVO terms expanded, tiers per data type, and the merge.

## Records

`search` starts from the biome's survey records: MAGs, MGnify bins, metagenome assemblies and
read runs. Assemblies and runs whose sample has no evidence for the config's kingdom
(`kingdom_match` false in the survey) are left out. Lower-level records of a sample that
already has more processed data stay in the table; `prepare` marks them `not_selected`.

A genome's biome is the biome of the samples it was assembled from, so its `sample_accession`
holds its source samples (`derived from`, `;`-joined for co-assemblies) as primary IDs
(SAMEA/SAMN/SAMD), the same IDs its sample's assemblies and runs carry.

## Taxonomy

Each row keeps the taxon as ENA shows it (`reported_tax_id`, `reported_tax_name`) and its lineage
in one NCBI taxonomy release (`tax_lineage`, `domain`, `taxdump_date`), from the
[taxdump archive](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/). The release is
downloaded once (about 150 MB) into the cache (`~/.cache/biome-finder/files/`); only the rows of
the table's tax IDs are read from it.

- A merged tax ID (`merged.dmp`) is followed to the current one. ENA keeps the reported tax ID
  and name current, so the resolved ones aren't stored: names differ only by NCBI renames of the
  same taxon (e.g. `Binatia bacterium` → `Binatia bacterium (SeqCode)`).
- A genome whose resolved lineage is outside the config's kingdom, or under an excluded taxon
  (Viruses, Metazoa, Tracheophyta, Phaeophyceae, plus `exclude_taxa_extra`), is flagged
  `out_of_scope_taxon`.
- A tax ID that doesn't resolve is flagged `taxid_deleted` (in `delnodes.dmp`) or
  `taxid_unknown` (not in the release, e.g. newer than it), and is left to the curator: ENA can
  still use a tax ID NCBI deleted. In freshwater-water-prok, 13 MAGs have one (e.g. 3691919,
  *uncultured Methylopumilus sp.*, deleted in NCBI 2026-11-01 and current in ENA).
- Assemblies and runs carry their metagenome taxon (e.g. `freshwater metagenome`), so the
  organism scope doesn't apply to them, and their `domain` is empty. Survey finds them by that
  taxon (`tax_eq`), so a record without its own taxon in the survey takes the one it was
  found by, with its tax ID from `survey.tsv`.

`tax_lineage` holds the eight NCBI ranks (`rankedlineage.dmp`) with a prefix each, the taxon's
own rank included, and an empty name where the lineage has no taxon at that rank:

```text
d__Bacteria;k__Pseudomonadati;p__Pseudomonadota;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae;g__;s__uncultured Neisseriaceae bacterium
d__;k__;p__;c__;o__;f__;g__;s__freshwater metagenome
```

The prefixes are `d` (domain; realm for viruses), `k` (kingdom), `p`, `c`, `o`, `f`, `g`, `s`.
NCBI gives metagenomes the species rank only.

## Sample metadata and hard filters

`search` reads each sample's fields from ENA (title, description, isolation source, the MIxS
`env_broad_scale` / `env_local_scale` / `env_medium`, organism, host, body site, coordinates).
`reported_biome` lists the biome fields as found, `key=value; …`; for a co-assembled genome,
values that differ between its source samples are `|`-joined:

```text
env_broad_scale=freshwater habitats; env_local_scale=lake; env_medium=water; organism=freshwater metagenome
```

The hard filters (curation action 0) match the sample's title, description, isolation source
and env fields. A row matching one, or flagged `out_of_scope_taxon`, goes to tier `D`; a genome
matches if any of its source samples does.

| Flag | Matches | Not matched |
|---|---|---|
| `control` | technical controls: negative, blank, extraction, kit, reagent or no-template controls, NTC | an experiment's untreated arm ("control microcosm"), which is a real sample: the biome evidence decides it |
| `mock` | mock communities | "mock sample for coassembly", ENA's virtual sample for a co-assembly |
| `synthetic` | synthetic communities | |
| `enrichment` | enrichments, enriched microbiomes or cultures | |
| `spiked` | spiked samples, spike-ins | |

The patterns were checked against freshwater-water-prok (2026-11): 961 runs, 114 MAGs and one
assembly went to tier D, mostly enrichments (e.g. PRJNA1177841, "culture enrichment of surface
water microbiome") and JGI's spiked Walnut Creek samples. Amplicon data and organelle genomes
are not filtered yet.

## Biome score and tiers

Each sample is scored on its own metadata plus its record's study. Evidence items are
`action:source:claim:quote`, where the claim is `in` (the config's biome), `out` (an excluded
or other biome) or `excluded_envo`:

| Action | Source | Claim | Weight |
|---|---|---|---|
| `sample_metadata` | an `env_*` field naming an `envo_include` term or one under it, by ENVO ID or exact name (label or synonym) | in | 0.8 |
| `sample_metadata` | an `env_*` field naming an `envo_exclude` term or one under it | excluded_envo | tier D |
| `sample_metadata` | a `keywords_include` word in the sample's text (title, description, isolation source, env fields) | in | 0.4 |
| `sample_metadata` | a `keywords_exclude` word in the sample's text | out | 0.6 (a conflict) |
| `project_description` | a `keywords_include` word in the study's title or description | in | 0.6 |
| `project_description` | a `keywords_exclude` word there | out | 0.3: study texts name other habitats in passing |
| `metagenome_name` | the sample's organism is one of `metagenome_taxa` | in | 0.3: often a default |
| `metagenome_name` | the organism is another metagenome | out | 0.3 |

Keywords are matched at a word start (`estuar` matches "estuarine") and skip the metagenome
taxon names, which default titles only repeat ("Metagenome or environmental sample from
freshwater metagenome"), and ENA's "mock sample for coassembly".

- **Score** (0–100): agreeing items combine as `1 − Π(1 − w)`, opposing ones scale that down
  the same way. Action 3 (sibling samples) adds up to ±10: the share of the study's samples
  whose own metadata points in, minus the share pointing out.
- **Tier**: `A` ≥ 80, `B` 60–79, `C` below 60 or with a conflict (an opposing item of weight
  ≥ 0.6), `D` with an excluding flag. `needs_curation` is tier C or a conflict;
  `curated_biome` is the config's `catalogue_biome` for tiers A and B.
- **Genomes** take their worst source sample, since a co-assembled genome is accepted only if
  all its source samples are.

The weights are PROPOSAL.md's starting values, not yet calibrated against a curated
freshwater set. In freshwater-water-prok (2026-11): MAGs 89% A, 9% C; assemblies 74% A, 16% C;
runs 37% A, 33% C; 14% of all rows in tier C.

ENVO shares synonyms between terms ("lake" is a synonym of pond, reservoir and saline
evaporation pond), so a value matches by ENVO ID first, then by exact label, then by exact
synonym; the strongest match decides between `envo_include` and `envo_exclude`, and a tie
excludes.

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
| `tax_lineage` | The eight NCBI ranks, `d__…;k__…;p__…;c__…;o__…;f__…;g__…;s__…`, empty after the prefix where there is no taxon at that rank |
| `domain` | `prokaryotes` or `eukaryotes` (the config's `kingdom` values), from the resolved lineage; empty for metagenomes |
| `taxdump_date` | The release used (`YYYY-MM-01`) |
| `flags` | `in_existing_catalogue` (already in the MGnify catalogue), the hard filters, `excluded_envo`, `out_of_scope_taxon`, `taxid_deleted`, `taxid_unknown`, and for runs `long_read`, `single_end`, `low_reads` (fewer than 200,000 reads) |
| `reported_biome` | The samples' biome fields as found in ENA |
| `curated_biome` | The config's `catalogue_biome`, for tiers A and B |
| `score` | The biome score, 0–100 (the worst source sample for a genome) |
| `tier` | `A`, `B`, `C`, or `D` for rows with an excluding flag |
| `needs_curation` | Tier C or a conflict |
| `evidence` | The evidence items, `action:source:claim:quote`, `; `-joined |
| `read_stats` | Runs: read count, layout, platform |
| `first_public` | Date the record became public |
| `curation_note`, `decision` | For the curator: `decision` is `keep` or `skip` and overrides the tier |
| `destination`, `local_path`, `prepare_status` | Written by `prepare` |
