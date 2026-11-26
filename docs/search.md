# search

`search` builds `samples.tsv` for one biome: one row per public INSDC data record (genome,
metagenome assembly or read run) that could go into its catalogue. The curator decides the
rows, and `prepare` turns the accepted ones into workflow inputs.

```{note}
In development (milestone M2). This version builds the table from the survey records and
merges it with an existing `samples.tsv`, with the taxonomy resolved against NCBI, each
sample's reported biome, the hard filters (tier D) and the rule-based biome score (curation
actions 1–5, tiers A–C). Coordinates and isolate genomes come next; their columns are already in the table, empty.
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
| `--gold` | The GOLD data export (`goldData.xlsx`), for GOLD ecosystem evidence; see [GOLD ecosystems](#gold-ecosystems). Without it, `search` runs without GOLD. |
| `-v` / `--verbose` | Before the command (`biome-finder -v search …`): also log every HTTP request, cached or not. |
| `--taxdump` | NCBI taxonomy release (`YYYY-MM-DD`): the monthly `new_taxdump_<YYYY-MM>-01.zip` of that month. Default: the latest. |

Progress goes to stderr, one timestamped line per step: records read, taxonomy flags, samples
found in ENA, hard-filter counts, ENVO terms expanded, tiers per data type, and the merge.

## Records

`search` starts from the biome's survey records: MAGs, MGnify bins, metagenome assemblies and
read runs. Assemblies and runs whose sample has no evidence for the config's kingdom
(`kingdom_match` false in the survey) are left out. Lower-level records of a sample that
already has more processed data stay in the table; `prepare` marks them `not_selected`.

A genome's biome is the biome of the samples it was assembled from, so its `sample_accession`
holds its source samples (`derived from`, `;`-joined for co-assemblies) as primary IDs
(SAMEA/SAMN/SAMD), the same IDs its sample's assemblies and runs carry. A genome assembled
from ENA's virtual co-assembly sample ("mock sample for coassembly … of samples: ERS…,ERS…")
gets the samples that sample lists instead: in freshwater-water-prok, 41 virtual samples stand
for 272 samples, behind about 4,100 MAGs.

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
and env fields. A row matching one (except `enrichment`, see below), or flagged `out_of_scope_taxon`, goes to
tier `D`; a genome
matches if any of its source samples does.

| Flag | Matches | Not matched |
|---|---|---|
| `control` | technical controls: negative, blank, extraction, kit, reagent or no-template controls, NTC | an experiment's untreated arm ("control microcosm"), which is a real sample: the biome evidence decides it |
| `mock` | mock communities | "mock sample for coassembly", ENA's virtual sample for a co-assembly |
| `synthetic` | synthetic communities | |
| `enrichment` | enrichments, enriched microbiomes or cultures: **tier C, not D**. An enrichment culture is fine if no organisms were added while culturing, which only a curator can tell | |
| `spiked` | spiked samples, spike-ins | |
| `virome` | the viral fraction: virome, virion, virus-like particles (VLP), viral fraction or metagenome. Also matched in the study's title and description, since a virome's sample often has only a default title (e.g. the Yangtze River virome, PRJNA672207); a study mixing viromes and cellular metagenomes would lose the latter too | "viral genome information" |
| `single_cell` | JGI's "Sorted cell/s from …" samples: single-cell (SAG) data, excluded by default | |
| `culture` | a culture ("non-axenic culture", "cell culture"): **tier C**, like an enrichment, which also stands for it | "culture-independent" |

Runs get two more flags from their ENA library and their sample, both tier D: `amplicon_only`
when the library strategy is AMPLICON or its selection is PCR (amplicon data that the survey's
WGS/OTHER query lets through: e.g. the Earth Microbiome Project runs, OTHER/PCR with under
150k reads), and `sample_low_reads` when none of the sample's runs reaches 200,000 reads, so
its reads can't be routed. The library also goes to `data_type_evidence`.

The patterns were checked against freshwater-water-prok (2026-11): 380 rows are spiked (mostly
JGI's Walnut Creek samples, "spiked with Pseudomonas putida & Escherichia coli"), 45 technical
controls and 5 synthetic communities; 646 rows are enrichments, sent to curation (e.g.
PRJNA1177841, "culture enrichment of surface water microbiome"). Amplicon data and organelle genomes
are not filtered yet.

## Biome score and tiers

Each sample is scored on its own metadata plus its record's study. Evidence items are
`action:source:claim:quote`, where the claim is `in` (the config's biome), `out` (an excluded
or other biome) or `excluded_envo`:

| Action | Source | Claim | Weight |
|---|---|---|---|
| `third_party_biome` | the biome lineage MGnify's curators gave the study (API v2 `/studies/`), under `mgnify_lineages` and none of `mgnify_lineages_exclude` | in | 0.9 |
| `third_party_biome` | an ancestor of a target lineage (e.g. `…:Freshwater` for `…:Freshwater:Lake`) | in, specificity = its depth over the target's (0.8) | 0.9 |
| `third_party_biome` | any other lineage (`…:Wetlands:Bog`, `…:Freshwater:Sediment`) | out | 0.9 (a conflict) |
| `third_party_biome` | the sample's GOLD ecosystem (`--gold`), against `gold_ecosystems` minus `gold_ecosystems_exclude`, the same way; trailing `Unclassified` levels are dropped first | in / out | 0.9 |
| `sample_metadata` | an `env_*` field naming an `envo_include` term or one under it, by ENVO ID or exact name (label or synonym) | in | 0.8 |
| `sample_metadata` | an `env_*` field naming an `envo_exclude` term or one under it | excluded_envo | tier D |
| `sample_metadata` | a `keywords_include` word in the sample's text (title, description, isolation source, env fields) | in | 0.4 |
| `sample_metadata` | a `keywords_exclude` word in the sample's text | out | 0.6 (a conflict), flagged `off_biome_keyword` |
| `project_description` | a `keywords_include` word in the study's title or description | in | 0.6 |
| `project_description` | a `keywords_exclude` word there | out | 0.3: study texts name other habitats in passing |
| `metagenome_name` | the sample's organism is one of `metagenome_taxa` | in | 0.3: often a default |
| `metagenome_name` | the organism is another metagenome | out | 0.3 |

A lineage is first mapped through `lineage_equivalents`. MGnify's `root:Mixed`, for studies
spanning several biomes, says nothing about a sample and is ignored.

Keywords are matched at a word start (`estuar` matches "estuarine") and skip the metagenome
taxon names, which default titles only repeat ("Metagenome or environmental sample from
freshwater metagenome"), and ENA's "mock sample for coassembly".

- **Score** (0–100): agreeing items combine as `1 − Π(1 − w)`, opposing ones scale that down
  the same way. Action 3 (sibling samples) adds up to ±10: the share of the study's samples
  whose own metadata points in, minus the share pointing out.
- **Tier**: `A` ≥ 80, `B` 60–79, `C` below 60 or with a conflict (an opposing item of weight
  ≥ 0.6), `D` with an excluding flag. `needs_curation` is tier C (tier D is excluded, not curated);
  `curated_biome` is the config's `catalogue_biome` for tiers A and B.
- **Off-biome samples** go to tier `D` (`off_biome`): the sample's own evidence (its metadata,
  organism and GOLD ecosystem, not its study's) points to another biome at conflict weight
  (≥ 0.6), and none of it points to the config's at that weight. A study's text or MGnify
  biome doesn't rescue a sample, since a study can mix biomes, and neither does a generic
  keyword or the default metagenome name (e.g. GOLD `… > Lake > Sediment` against "freshwater"
  in the isolation source). A strong item both ways, such as an ENVO term against GOLD, is a
  conflict: tier C. In freshwater-water-prok this excludes 1,346 MAGs whose source samples are salt or
  soda lakes, ocean, wastewater or sediment (e.g. SAMN07376294: GOLD `… > River > Sediment`,
  organism `freshwater sediment metagenome`).
- **Genomes** take their worst source sample, since a co-assembled genome is accepted only if
  all its source samples are.

The weights are PROPOSAL.md's starting values, not yet calibrated against a curated
freshwater set. In freshwater-water-prok (2026-11): MAGs 88% A, 9% C; assemblies 75% A, 17% C;
runs 31% A, 39% C; 16% of all rows in tier C. 35 of its 4,219 studies have an MGnify biome.

ENVO shares synonyms between terms ("lake" is a synonym of pond, reservoir and saline
evaporation pond), so a value matches by ENVO ID first, then by exact label, then by exact
synonym; the strongest match decides between `envo_include` and `envo_exclude`, and a tie
excludes.

## GOLD ecosystems

JGI's GOLD gives many BioSamples a curated five-level ecosystem path (e.g. `Environmental >
Aquatic > Freshwater > Lake > Unclassified`). GOLD publishes it only as an Excel export, on a
page that checks for a browser, so the curator downloads it once:

1. Open <https://gold.jgi.doe.gov/downloads> in a browser (no login needed; the first visit
   may show a short "Just a moment…" check).
2. Download **Public Studies/Biosamples/SPs/APs/Organisms Excel** (about 240 MB,
   `goldData.xlsx`). The other three files (SRA-based biome APs, GOLD CVs, the ecosystem
   paths) aren't needed.
3. Pass it to `search`: `--gold goldData.xlsx`. The file is git-ignored in the project folder.

The first run reads the export (about three minutes) and caches its BioSample → ecosystem table
in `~/.cache/biome-finder/files/gold-<key>.tsv`; later runs read the cache, until the file
changes. A sequencing project links an NCBI BioSample to a GOLD Biosample (metagenomes) or
Organism (isolates), which holds the path. Download a new export now and then: the page shows
when it was last generated.

In freshwater-water-prok (export of 2026-11-25), 10,406 of its 25,229 samples have a GOLD
ecosystem, mostly `… > Freshwater > Lake` (5,909) and `… > River` (1,893). With it, MAGs and
assemblies are 91% tier A, and runs 45% A, 34% C (from 31% A, 39% C without GOLD).

## Benchmark

`benchmark/freshwater-water-prok_benchmark150.tsv` (local, git-ignored, like the curated files
below) holds 150 rows drawn from the freshwater
`samples.tsv` (stratified by tier and data type, at most four per study; seed 20261126), with
the tool's predictions apart in `…_benchmark150_tool.tsv`, so curators work blind. They were
curated with the curate-sample skill (`curated_freshwater-water-prok.tsv`) and, 40 of them, by
hand (`curated_freshwater-water-prok-filipe.tsv`). On 2026-11-26, against the skill's 150:
tier A precision 47/47, tier A+B 64/67, tier D 64/64 skips, tier C 19/150 (13%), no kept record
in tier D, and data type 150/150 (SAGs as `single_cell`). The three A/B rows a curator skipped
lack the evidence in their own metadata: a sample from a mesocosm study that doesn't say so, a
salmon aquaculture environment with no sample description, and a co-culture whose sample names
only the cyanobacterium and the inoculum. The 40 curated by hand agree with the skill's decisions on all 40 (after review), with
tier A 12/12 and tier D 20/20.

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
| `data_type_evidence` | What decided it: `ena.assembly_type=…` for genomes and assemblies (an untyped assembly found by its sample also shows its metagenome organism), the library for runs |
| `source_db` | `ENA` |
| `sample_accession` | The record's sample; for genomes, their source samples |
| `study_accession` | The record's ENA study (PRJ*) |
| `reported_tax_id`, `reported_tax_name` | As shown in the ENA record |
| `tax_lineage` | The eight NCBI ranks, `d__…;k__…;p__…;c__…;o__…;f__…;g__…;s__…`, empty after the prefix where there is no taxon at that rank |
| `domain` | `prokaryotes` or `eukaryotes` (the config's `kingdom` values), from the resolved lineage; empty for metagenomes |
| `taxdump_date` | The release used (`YYYY-MM-01`) |
| `flags` | `in_existing_catalogue` (already in the MGnify catalogue), the hard filters, `excluded_envo`, `off_biome` (tier D: the sample's own evidence points only to another biome), `amplicon_only` and `sample_low_reads` (runs, tier D), `single_cell` (tier D), `culture` (tier C), `off_biome_keyword` (the sample's own text has a `keywords_exclude` word: a conflict, tier C unless another flag excludes the row), `out_of_scope_taxon`, `taxid_deleted`, `taxid_unknown`, and for runs `long_read`, `single_end`, `low_reads` (fewer than 200,000 reads) |
| `reported_biome` | The samples' biome fields as found in ENA |
| `third_party_biome` | The study's MGnify biome, `mgnify:<MGYS>:<lineage>`, and the samples' GOLD ecosystems, `gold:<Gb or Go ID>:<path>`, `; `-joined |
| `curated_biome` | The config's `catalogue_biome`, for tiers A and B |
| `score` | The biome score, 0–100 (the worst source sample for a genome) |
| `tier` | `A`, `B`, `C`, or `D` for rows with an excluding flag |
| `needs_curation` | Tier C, where every conflict lands; tier D rows are excluded, not curated |
| `evidence` | The evidence items, `action:source:claim:quote`, `; `-joined |
| `read_stats` | Runs: read count, layout, platform |
| `first_public` | Date the record became public |
| `curation_note`, `decision` | For the curator: `decision` is `keep` or `skip` and overrides the tier |
| `destination`, `local_path`, `prepare_status` | Written by `prepare` |
