# biome-finder — Project Proposal

**Status:** draft v0.6 · 2026-11-19
**Process reference:** MGnify genome catalogue curation guidelines, condensed in [Appendix A](#appendix-a--curation-guidelines) (rules), [B](#appendix-b--reference-queries) (queries) and [C](#appendix-c--mgnify-catalogues-snapshot-2026-11-19) (catalogue snapshot)

---

## 1. What biome-finder is

`biome-finder` organizes the **input data** for genome catalogue workflows. Given a target
biome, it produces one table, `samples.tsv`, and then the workflow inputs built from it:

1. **`samples.tsv`**: every public INSDC record (genome, metagenome assembly or read run)
   that could go into a catalogue for that biome. Each row carries the sample, project,
   publication, taxonomy, the biome as reported in the metadata, a pre-curated biome, a
   confidence score, and whether a human needs to check it. The same table later records
   what `prepare` did with each row (destination, local path, status). There is no second
   manifest.
2. **Workflow inputs**, routed by the most processed public data each sample has:

| `data_type` | What it is | Goes to |
|---|---|---|
| `mag`, `bin`, `isolate` | Single-organism genomes, at any assembly level (contig, scaffold, chromosome), used as they are | FASTA download → [genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) |
| `contigs` | Metagenome assemblies that were never binned (e.g. MGnify ERZ assemblies from metaSPAdes / MEGAHIT) | samplesheet → [genomes-generation](https://github.com/EBI-Metagenomics/genomes-generation) |
| `reads` | Raw metagenomic runs only | samplesheet → [miassembler](https://github.com/EBI-Metagenomics/miassembler) → genomes-generation |

The downstream workflows do the rest: assembly, binning, QC, dereplication, clustering
and annotation.

**Biome-agnostic by design.** All biome knowledge (lineages, ontology terms,
include/exclude rules, boundary decisions) lives in one YAML config per biome. The code
has no biome-specific logic. The first runs will target environmental biomes, but the same
tool must work for any environmental or host-associated biome.

**Organism scope:** Bacteria, Archaea and microbial Eukaryotes (protists, algae, fungi).
Excluded by default:

- Viruses and phages (10239)
- Animals (Metazoa, 33208)
- Vascular plants (Tracheophyta, 58023)
- Brown algae (Phaeophyceae, 2870). More generally, **all macroalgae** (multicellular,
  macroscopic seaweeds). Red and green macroalgae have no single excludable taxon, so
  they are judged per genus during curation.
- Plasmids and other mobile genetic elements
- Organelle genomes

The excluded taxa are a default list that a biome config can extend.

```
                ┌────────────────────────── biome-finder ────────────────────────────┐
biome.yaml ─────┤ survey ─► search ─► samples.tsv ◄──── curator (decision, PDFs)     │
                │                         │  ▲                                       │
                │                       prepare (writes status/paths back)           │
                │   mag / bin / isolate ──┼──► catalogue/ ───────────────────────────┼──► genomes-catalogue-pipeline
                │   (eukaryotes) ─────────┼──► catalogue/eukaryotes/datascout.csv ───┼──► datascout ──┘
                │   contigs ──────────────┼──► binning/<study>/ ─────────────────────┼──► genomes-generation ──┘
                │   reads ────────────────┴──► assembly/samplesheet.csv ─────────────┼──► miassembler
                │                    link-assemblies ◄── miassembler output ─────────┤       (then → binning/)
                └────────────────────────────────────────────────────────────────────┘
```

---

## 2. Step 0 — Choosing a biome: `survey`

```bash
biome-finder survey --biomes biomes/*.yaml --out survey/
```

The output is a folder: `survey.tsv` (one row per biome × metagenome taxon, plus TOTAL),
`biomes/<name>/records.tsv` and `studies.tsv` (every counted record, and per-study
counts) and `report.html` (an interactive dashboard: all catalogues, or one at a time).

For each biome config, it reports:

- the **existing MGnify catalogue** for that biome, if any: version, last update, species
  and genome counts (read from the MGnify API, not hard-coded) → candidate for **update**
  vs **new**
- the number of records of each `data_type` that match the biome rules (MAGs, bins,
  isolates, contigs, read runs with ≥ 1 M reads plus their total bases), split into
  records public before and after the last catalogue update
- distinct studies and linked publications

The survey only counts. The user writes the biome decision into the config, including
whether to split a biome or keep it whole ("one marine catalogue vs. marine water +
marine sediment", "rumen vs. whole gut").

**First intended targets (illustrative):** environmental biomes with no catalogue, such as
freshwater, freshwater sediment, permafrost and crop rhizospheres (wheat, rice), and
updates of existing ones (soil v1.0, barley rhizosphere v2.0). §9 lists the data sources
found for these. Nothing in the tool depends on this choice.

### 2.1 Biome config (`biomes/<name>.yaml`)

The YAML file is read with `yaml.safe_load` (PyYAML). It is both the tool's input and the
written record of the biome decision.

```yaml
name: freshwater-lake-sediment
kingdom: prokaryotes                  # prokaryotes | eukaryotes: MGnify builds one catalogue per kingdom (pipeline --kingdom)
mgnify_catalogue: ""                 # set for updates, e.g. soil-v1-0
catalogue_biome: "root:Environmental:Aquatic:Freshwater:Sediment"   # → genomes-catalogue-pipeline --biome (must exist in MGnify)
mgnify_lineages:                     # whole subtrees count (each lineage plus all its descendants)
  - "root:Environmental:Aquatic:Freshwater:Sediment"
mgnify_lineages_exclude: []         # subtrees removed from the scope, e.g. a :Sediment child of a water lineage
envo_include:                        # term + all its ENVO subclasses; illustrative: verify all IDs in OLS
  - ENVO:00000546                    # lake sediment
  - ENVO:01000252                    # freshwater lake biome (+ env_medium = sediment)
envo_exclude:                        # subtrees removed (term + subclasses), checked on every env_* field
  - ENVO:03000033                    # marine sediment
gold_ecosystems:                     # GOLD 5-level paths (also used by awesome-mag)
  - "Environmental > Aquatic > Freshwater > Lake > Sediment"
gold_ecosystems_exclude: []          # GOLD paths (and everything under them) removed
metagenome_taxa: ["freshwater sediment metagenome"]
keywords_include: ["lake sediment", "lacustrine sediment"]
keywords_exclude: ["marine", "estuar", "wastewater", "bioreactor", "enrichment",
                   "mesocosm", "microcosm", "spiked", "mock", "synthetic", "control"]
host: none                           # environmental: any host field is a conflict
# host: "Bos taurus"
host_tax_id:                         # NCBI tax ID of the host (e.g. 9913)
host_keywords: []                    # free-text host names (e.g. ["cow", "Bos"]); ENA hosts are recorded inconsistently,
                                     # so samples match host_tax_id OR any keyword in the `host` field
# body_site_include: ["rumen"]
contaminant_reference: ""            # host genome for miassembler decontamination (host-associated)
data_types: [mag, isolate, contigs, reads]   # add bin only in the cases of Appendix A.3; sag only if explicitly wanted
exclude_taxa_extra: []               # added to the default exclusions (Viruses, Metazoa, Tracheophyta, Phaeophyceae)
lineage_equivalents: {}              # evidence lineage → lineage it counts as; e.g. for a marine water biome:
                                     # {"root:Environmental:Aquatic:Marine:Coastal": "root:Environmental:Aquatic:Marine:Oceanic"}
ncbi_isolate_query: ""               # NCBI BioProject query for isolates (A.4), e.g. 'bacteria[Organism] AND "lake sediment"'
ncbi_wgs_terms: []                   # host / isolation_source terms for the NCBI WGS search (B.4)
min_study_yield: 0.3                 # mixed-biome read studies: skip if accepted runs / all runs is below this (A.5)
indicpiper_habitats: ["freshwater sediment"]  # IndicPiper habitat(s) that count as this biome (§3.9)
notes: "Surface and core sediments. Saline/soda lakes excluded."
```

---

## 3. Step 1 — The table: `search`

```bash
biome-finder search --biome biomes/<name>.yaml --survey survey/ --out samples.tsv \
    [--gold goldData.xlsx] [--since YYYY-MM-DD] [--taxdump YYYY-MM-DD] [--papers-dir papers/] [--llm] [--composition [all]]
```

**Re-running is safe.** If `samples.tsv` already exists, `search` merges by
`record_accession`. It refreshes the automatic columns and keeps the human ones
(`decision`, human-written `curation_note`). The normal loop is: search → curate → add
PDFs → search again.

**Input.** `search` starts from the biome's survey records
(`survey/biomes/<name>/records.tsv`: MAGs, MGnify bins, assemblies and runs with evidence for
the config's kingdom) and adds isolate genomes, which `survey` doesn't search.

### 3.1 Data sources (programmatic)

Everything below is public and needs no credentials, except the optional NCBI API key and
the LLM key.

| Source | What we take | Access |
|---|---|---|
| ENA Portal API | MAG layer (`wgs_set`), bins and primary metagenome assemblies (`analysis`), read runs, sample and study fields, `assembly_type` | `/ena/portal/api/search` (queries in Appendix B) |
| ENA Browser / XML | Full sample attributes, study description, `derived from` links, publications | `/ena/browser/api/xml/<acc>` |
| ENA cross-references | Publications describing or citing a record; assemblies/MAGs derived from a read study | `/ena/xref/rest/` |
| NCBI | WGS sets, BioSample attributes and package (MIMAG/MISAG/MIGS), BioProject, `assembly_summary_genbank.txt` (GCA by BioSample, `excluded_from_refseq`) | E-utilities + FTP |
| NCBI Taxonomy | Monthly `new_taxdump_YYYY-MM-01.zip` from `pub/taxonomy/taxdump_archive/` (`merged.dmp`, `delnodes.dmp`, lineages) | FTP |
| MGnify API | Existing catalogues and their genomes (update dedup); studies with a curated MGnify biome lineage; MGnify-generated ERZ assemblies | `/metagenomics/api/v2/genomes/catalogues/` for catalogues (v1 lists outdated versions); `/metagenomics/api/v2/studies/` for the studies' curated biomes (action 1); `/metagenomics/api/v1/biomes/` to browse the lineage tree |
| Metalog (EMBL) | Manually curated sample metadata, where available | Downloads |
| JGI GOLD | Curated 5-level ecosystem classification of BioSamples | Public data export ("Public Studies/Biosamples/SPs/APs/Organisms Excel"), downloaded by the curator from a browser (the page is behind a browser check) and passed as `--gold` |
| Europe PMC | Paper metadata, open-access full text (Methods, Data availability), accession text-mining | REST API |
| SPIRE + mag-provenance | Flags SPIRE bins and re-assemblies of the same source samples | [mag-provenance](https://github.com/EBI-Metagenomics/mag-provenance/tree/dev) |
| Sandpiper (SingleM) | GTDB genus-level community profile per public run, for IndicPiper (§3.9) | `sandpiper.qut.edu.au/api/condensed_csv_with_extras/<run>?taxonomy_type=gtdb` |
| MGnify analyses | SSU/LSU OTU taxonomy files of runs/assemblies MGnify already analysed, for Trapiche (§3.9) | MGnify API downloads |

### 3.2 Seed collections (from [awesome-mag](https://github.com/inspirewind/awesome-mag))

These are used to find studies and papers. A seed is used only through its INSDC
accessions. The routes below come from awesome-mag's access notes and must be confirmed
per study.

| Route | Resources (scope) | Use |
|---|---|---|
| **Genomes in INSDC** | [MGnify Genomes](https://www.ebi.ac.uk/metagenomics/browse/genomes/) (multi-biome, ENA) · [OceanDNA](https://www.nature.com/articles/s41597-022-01392-5) (marine, DDBJ PRJDB11811) · [Glacier-fed Streams MAGs](https://www.nature.com/articles/s41564-024-01874-9) (stream sediment, NCBI) · [GROWdb](https://www.nature.com/articles/s41586-024-08240-z) (rivers, NCBI) · [RGMGC](https://rummeta.njau.edu.cn/) (ruminant GIT, NCBI Assembly) · [SPIRE](https://spire.embl.de/) (multi-biome, ENA MAG layer, actually bins) · [Human Gut Archaeome](https://www.nature.com/articles/s41564-021-01020-9) (check) | `mag`/`bin` rows directly |
| **Reads or assemblies in INSDC, genomes elsewhere** | [LakePulse](https://www.nature.com/articles/s41564-023-01435-6) (lakes, ENA reads) · [RUG2](https://www.nature.com/articles/s41587-019-0202-3) (rumen, ENA reads + assemblies) · [Meta320](https://link.springer.com/article/10.1186/s40168-024-01806-z) (sheep/goat, SRA) · [GEM](https://genome.jgi.doe.gov/portal/GEMs/GEMs.home.html) (multi-biome, IMG) · [SMAG](https://microbma.github.io/project/SMAG.html) (soil, Zenodo) · [Bin Chicken RBG](https://zenodo.org/records/15220963) (multi-biome, Zenodo) · [cFMD](https://github.com/SegataLab/cFMD) (food, Zenodo) · [Microbiome Datahub](https://mdatahub.org/) (multi-biome, NIG) · [mOTUs DB](https://motus-db.org/) · [OMDB](https://omdb.microbiomics.io/about) (ocean) · [HRGM](https://www.decodebiome.org/HRGM/) / [HROM](https://www.decodebiome.org/HROM/) / [MRGM](https://www.decodebiome.org/MRGM/) (human gut/oral, mouse gut) | Study lists → `contigs`/`reads` rows (regenerate MAGs) |
| **Not INSDC** (CNGB/CNSA, NGDC GSA/GWH, NMDC-China) | [GOMC](https://db.cngb.org/maya/datasets/MDB0000002) · [TPMC](https://ngdc.cncb.ac.cn/tpmc) · [QXLSG](https://www.sciencedirect.com/science/article/pii/S2211124725012549) (lake sediment) · [PIGC](https://db.cngb.org/data_resources/project/CNP0000824) · [Caprinae gut](https://journals.asm.org/doi/10.1128/spectrum.02211-22) · [UHSG](https://doi.org/10.1002/advs.202300050) · [ICRGGC](https://nmdc.cn/icrggc/) | Excluded unless the same data is mirrored in INSDC |
| **Aggregator portals** | [MAGdb](https://magdb.nanhulab.ac.cn/) (offline) · [gcMeta](https://gcmeta.wdcm.org/) · [Anammox catalog](https://figshare.com/articles/dataset/A_comprehensive_catalog_encompassing_genes_and_genomes_reveals_the_core_community_and_functional_diversity_in_anammox_microbiota/25476583) | Only for finding source studies |
| **Out of scope (viral)** | GPD · MGV · UHGV · MetaVR | — |

### 3.3 Row unit: the data record. Curation unit: the sample

Each row of `samples.tsv` is one **data record**: a genome, a metagenome assembly or a
read run, together with its **sample** (SAMEA/SAMN/SAMD).

- The biome score and tier are computed **per sample** and repeated on each of its rows.
  A biome belongs to a sample, not to a project.
- `decision` is **per row**. A curator can reject one genome, for example a SAG, without
  rejecting the sample.
- A co-assembled genome lists all its source samples. It is accepted only if **all** of
  them are accepted.

### 3.4 Data type: MAG, isolate or something else

The data type is not always explicit, so it gets its own evidence and confidence. In
`marine_euks.tsv`:

- 689 of the 757 isolates have an empty `Assembly Type`
- 135 records ended up as `Unknown`
- 6 were `NEEDS CONFIRMATION`; the paper showed they were SAGs

Evidence, in order:

| Field / source | Indicates |
|---|---|
| ENA `assembly_type` | `Metagenome-Assembled Genome (MAG)` → mag · `binned metagenome` → bin · `primary metagenome` → contigs · `clone or isolate` → isolate (still needs pure-culture or tissue evidence). Often empty for isolates and SAGs |
| ENA assembly `DESCRIPTION` | "cell sorted by flow cytometer", "amplified" (MDA/WGA) → sag · "binned", "metagenome" → mag/bin · "strain", "culture" → isolate candidate |
| BioSample `mixed sample` / `member sample` | `in silico mixture` → pooled source samples (common for Tara SAGs/MAGs) |
| NCBI `excluded_from_refseq` | `derived from metagenome` → mag/bin · `derived from single cell` → sag |
| BioSample package | `MIMAG.*` → mag · `MISAG.*` → sag · `MIGS.ba/eu` → candidate isolate |
| BioSample attributes | `metagenome-source`, `derived_from` → mag. Source sample organism under "metagenomes" (taxid 408169) → metagenome-derived |
| Strain / culture collection IDs (ATCC, NIES, CCAP, DSM…) | Supports isolate (it never proves the isolation source) |
| Paper Methods (LLM) | "pure culture", "axenic", "single-cell sorted", "binned", "co-culture" |

Rules (Appendix A.4):

- Don't trust the word "isolate".
- `isolate` is a single-organism genome from either:
  - a **pure culture**, with evidence at **sample** level or in the paper; or
  - **tissue or whole individuals of a wild-collected organism**, not a metagenome
    (BioSample `tissue`/`dev_stage`, or a paper saying "wild-collected individuals").
    The organism must still be in scope; macroalgae are out.

  Otherwise the record becomes `mag` (if metagenome-derived) or `unknown`
  (→ `needs_curation`).
- SAGs (`sag`) are excluded by default.
- `contigs` means an unbinned **metagenome** assembly. A single-organism genome at contig or
  scaffold level is still `mag`/`isolate`. The `assembly_level` column records that level
  separately.

The `Category` column of `marine_euks.tsv` (Isolate / MAG / Unknown / NEEDS CONFIRMATION)
is the benchmark for data-type assignment.

### 3.5 Taxonomy

Each row keeps the taxonomy **as reported** and its lineage **in a given NCBI Taxonomy
release**. ENA keeps the reported tax ID and name current (in freshwater-water-prok, 2026-11,
no tax ID resolves to another and names differ only by NCBI renames of the same taxon), so
no resolved tax ID or name is stored:

- `--taxdump YYYY-MM-DD` picks `new_taxdump_<YYYY-MM>-01.zip` from the NCBI archive; the
  default is the latest release.
- Merged tax IDs are followed to the current one (`merged.dmp`). Deleted ones are flagged `taxid_deleted`
  (`delnodes.dmp`), and IDs missing from the release `taxid_unknown`. Neither is put out of
  scope: ENA can still use a tax ID NCBI deleted, so the curator decides.
- `domain` (prokaryote / eukaryote) and the out-of-scope taxa filter use the **resolved**
  lineage.
- For `contigs`/`reads`, the taxon is the sample's metagenome taxon (e.g. "lake sediment
  metagenome"). That taxon is also weak biome evidence (action 5).

### 3.6 `samples.tsv` columns

| Column | Content |
|---|---|
| `record_accession` | Genome (ENA set / GCA), ERZ assembly, or run accession |
| `data_type` | `mag` · `bin` · `isolate` · `sag` · `contigs` · `reads` · `unknown` |
| `data_type_evidence` | Fields that decided the type (§3.4), e.g. `ena.assembly_type=clone or isolate; biosample.package=MIGS.eu; paper:"axenic culture"` |
| `assembly_level` | contig / scaffold / chromosome / complete (genomes only) |
| `source_db` | ENA / NCBI / DDBJ |
| `sample_accession` | BioSample. For co-assemblies: `;`-separated source samples |
| `related_accessions` | Runs behind a `contigs` record; ERZ/MAGs derived from a `reads` record |
| `study_accession`, `secondary_study_accession`, `bioproject` | PRJ* / ERP·SRP·DRP |
| `reported_tax_id`, `reported_tax_name` | As shown in the database record |
| `tax_lineage`, `domain` | Resolved against NCBI Taxonomy at `taxdump_date`. `tax_lineage` holds the eight NCBI ranks, `d__…;k__…;p__…;c__…;o__…;f__…;g__…;s__…`, empty after the prefix where the lineage has no taxon at that rank |
| `taxdump_date` | The taxdump release used |
| `publication` | DOI / PMID (`;`-separated) and the link type (`describes` vs `cites`) |
| `paper_access` | `open_access` · `pdf_provided` · `needs_pdf` · `not_needed` · `no_paper` |
| `paper_pdf` | Expected local file name when a PDF is needed or provided (§3.11) |
| `reported_biome` | Raw metadata as found: `env_broad_scale`, `env_local_scale`, `env_medium`, `isolation_source`, organism, host, body site |
| `third_party_biome` | MGnify lineage, GOLD ecosystem or Metalog classification, if any |
| `curated_biome` | Biome lineage proposed by the scorer |
| `score` | 0–100 confidence that the **sample** belongs to the target biome |
| `tier` | `A` auto-accept · `B` accept + spot-check · `C` needs curation · `D` excluded |
| `needs_curation` | bool: tier C (where any conflict and `data_type = unknown` land); tier D rows are excluded, not curated |
| `evidence` | Ordered `action:source:claim` items, e.g. `sample_metadata:ENVO:00000546:lake sediment; paper_methods:PMC123:"cores were taken from…"` |
| `flags` | `control`, `synthetic`, `mock`, `enrichment`, `spiked`, `amplicon_only`, `single_end`, `low_reads`, `long_read`, `organelle`, `out_of_scope_taxon`, `off_biome_keyword`, `taxid_deleted`, `taxid_unknown`, `no_assembly`, `duplicate`, `spire_bin`, `same_study_mag_and_bin`, `in_existing_catalogue`, `host_conflict`, `composition_conflict` |
| `read_stats` | runs, layout, platform, min/max/mean reads (for `reads` and `contigs`) |
| `first_public` | Date (for incremental catalogue updates) |
| `curation_note` | Free text, written by the LLM or a human, in the style of the "Manual Inspection" column |
| `decision` | For the human curator: `keep` / `skip`. When filled, it overrides `tier` |
| `destination` | Written by `prepare`: `catalogue/prokaryotes/ena`, `binning/<study>`, `assembly`, … |
| `local_path` | Written by `prepare` / `link-assemblies`: downloaded FASTA or produced assembly |
| `prepare_status` | `fetched` · `fetch_failed` · `qs50_failed` · `queued_assembly` · `assembly_only` · `low_yield_study` · `queued_binning` · `assembly_qc_failed` · `deferred_long_read` · `not_selected` |

`--since` restricts the search to records that became public after a date (for catalogue
updates).

### 3.7 Scoring based on curation actions

The `Manual Inspection` column of `marine_euks.tsv` records what a curator actually did to
decide each biome. A rough keyword classification of the 848 non-empty notes:

| Action (what the curator did) | Notes | How often it was conclusive |
|---|---|---|
| Read the linked **paper** (methods/title) | 560 | ~46%. Many papers "only report in the title" or "do not inform the source" |
| Read **supplementary tables/figures** (sites, lat/lon) | 83 | High |
| Read **sample** metadata (BioSample page) | 66 | High when specific |
| Checked **organelle** status → excluded | 58 | n/a (filter) |
| Checked **assembly availability** → excluded | 31 | n/a (filter) |
| Looked up **culture collection / commercial strain** page (ATCC, NIES, CCAP, …) | 37 | ~43% |
| Read **project/study description** (incl. MGnify "samples from the following biomes: root:…") | 10+ | High |
| Checked **other samples of the same project** for consistency | several | Medium |

The scorer runs these actions in order, from cheapest to most expensive, and stops as
soon as the evidence is conclusive:

| # | Action | Automated by | Base weight (to calibrate) |
|---|---|---|---|
| 0 | **Hard filters**: public in INSDC, data available, in-scope taxon, not organelle, not control/mock/synthetic/enrichment/spiked, not amplicon-only. Controls are technical ones (negative, blank, extraction, kit); an experiment's untreated arm is a real sample left to actions 1–9 | rules | gate |
| 1 | Third-party curated biome (MGnify study lineage, GOLD, Metalog) | lookup | 0.9 |
| 2 | Structured sample metadata (ENVO `env_*`, `env_biome`, `isolation_source`, host, body site). ENA `env_tax:*`/`env_geo:*` tags are derived from the organism's taxonomy, not the sample, so they don't count (a sample can carry both `env_tax:marine` and `env_tax:freshwater`) | rules + term map | 0.8 when specific, 0.4 when generic, ~0 for `env_tax`/`env_geo` |
| 3 | Consistency with sibling samples in the same project | rules | ± 0.1 modifier |
| 4 | Project/study title and description | rules → LLM | 0.6 |
| 5 | Metagenome organism name (e.g. "freshwater metagenome") | rules | 0.3 (often a default value) |
| 6 | Paper Methods / Data availability (open access, or a provided PDF) | LLM | 0.8 if explicit, 0.3 if only title or region |
| 7 | Supplementary tables (sample sites, coordinates) | LLM | 0.8 |
| 8 | Culture collection record (isolates only) | LLM + web fetch | 0.7 |
| 9 | Coordinates sanity check (e.g. a "freshwater" sample in the open ocean) | rules | conflict flag only |
| 10 | **Community composition** of the sample's reads: Trapiche lineage prediction, IndicPiper indicator-habitat fractions (§3.9) | Trapiche (subprocess) + IndicPiper database | Trapiche 0.6 · IndicPiper 0.4 |

**Combining evidence.** Each item gives `(claimed_lineage, specificity, weight)`.
`specificity` is how deep the claim reaches into the target lineage: "Aquatic" only
partially supports "Aquatic:Freshwater:Sediment". Agreeing items are combined as
`1 − Π(1 − w·s)`; opposing items scale that down the same way, `× Π(1 − w·s)`. Conflicting
items with weight ≥ 0.6 set `needs_curation`, whatever the score. Action 3 adds up to ±0.1:
the share of the study's samples whose own metadata points in, minus the share pointing out.
A keyword match skips the metagenome taxon names, which default sample titles repeat and
action 5 already counts. An exclude keyword in a study's text weighs 0.3, below the conflict
threshold: study descriptions name other habitats in passing. Tiers: **A** ≥ 80 with no conflicts, **B** 60–79, **C** < 60 or any conflict,
**D** hard-filter fail. When unsure, exclude: that's guideline A.1 ("if you are not
convinced… leave it out").

**Calibration.** Weights and thresholds are fitted on `marine_euks.tsv`:

- **Biome:** `Curated Biome`, `Keep` and `Curator decision` are the ground truth.
- **Data type:** `Category` is the ground truth.
- **Targets:** tier A ≥ 95% concordant with the manual decision, and tier C ≤ 30% of
  records.
- **New biomes:** each biome curated with the tool adds its manually decided rows to the
  benchmark, so the weights don't fit a single biome.

### 3.8 LLM layer (`--llm`)

- **Used only for actions 4, 6, 7 and 8, plus the paper part of §3.4**, and only for
  records still uncertain after the rule-based actions.
- **Input:** the biome config plus the retrieved text: study description, paper Methods
  and Data-availability sections, supplementary table rows, or provided PDFs.
- **Output:** strict JSON with `claimed_lineage`, `specificity`, `data_type`
  (isolate/mag/sag/…), `is_control`, a **verbatim quote** and its source, and a
  `curation_note` sentence.
- **Guardrail against hallucination:** if the quote is not found in the retrieved text,
  the item is discarded. LLM-only evidence can't put a sample in tier A; it reaches B at
  most.
- **Cost control:** results are cached per accession and document hash; batch processing
  is used for large runs; model choice is set in the config, not hard-coded.

### 3.9 Composition-based biome evidence (`--composition`)

When the text and metadata leave a sample unclear, the microbial community itself is
another, independent source of evidence. Two tools provide it. Both need public reads
for the sample:

- **`reads` and `contigs` rows:** the sample's own runs.
- **`mag`/`bin` rows:** the runs of the source sample (via `derived from`).
- **Isolates:** not applicable.

**[Trapiche](https://github.com/Finn-Lab/trapiche)** (Finn lab, Python)

- **What it does:** predicts a deep biome lineage, in the same `root:…` vocabulary as
  MGnify and GOLD. It combines text-based labels with a taxonomy-based model (a
  community2vec embedding of the sample's taxonomic profile).
- **Text input:** our best text evidence (actions 1–8), passed as
  `ext_text_pred_project` / `ext_text_pred_sample`. Trapiche then uses the community
  profile to refine the lineage or contradict it. With no usable text, it falls back to
  its own BERT classifier on the project/sample description.
- **Taxonomy input:** MGnify-format taxonomy files (SSU/LSU OTU TSVs) through the MGnify
  API, so it is available only for runs MGnify has already analysed.
- **Execution:** heavy dependencies (TensorFlow, Hugging Face models), so it runs as a
  **subprocess in its own environment**: NDJSON in, NDJSON out (`trapiche input.ndjson`).
  It is not a dependency of biome-finder.
- **Evidence item:** `final_selected_prediction`, with its probability as the
  specificity factor.

**[IndicPiper](https://github.com/cliffbueno/IndicPiper)** (Bueno de Mesquita, R)

- **The database:** a pre-built table of GTDB genera that indicate one of 13 habitats.
  v2.1.0 has 767 genera for soil, marine sediment, marine water, mammalian gut,
  freshwater sediment, freshwater water, human oral, wastewater, glacier, hydrothermal
  vent, phyllosphere, human vaginal and human skin. Soil and rhizosphere are combined by
  default.
- **Profile:** the run's GTDB genus profile, fetched from the **Sandpiper** API
  (SingleM, pre-computed for public SRA/ENA runs), so nothing is computed locally.
- **Method:** exact genus match, then relative abundance summed per indicator habitat.
  That's a join and a sum, so we **reimplement it in Python** (no R dependency) and
  vendor the database CSV at a pinned version.
- **Mapping:** the biome config's `indicpiper_habitats` says which habitat(s) count as
  the target biome.
- **Evidence item:** "the target habitat dominates the indicator signal". The specificity
  factor is capped at the habitat level, so it can tell freshwater sediment from
  freshwater water, but not lake sediment from river sediment.

**Rules**

- **When it runs:** by default, only for samples still in tier C or with generic-only
  evidence after actions 0–9. `--composition all` also runs it on tier A/B samples,
  as a **mislabel check**.
- **Composition can't raise a sample to tier A on its own** (same rule as the LLM).
- **Composition contradicts strong metadata** (e.g. "lake sediment" in `env_medium`, but
  IndicPiper says marine sediment or soil): `composition_conflict` → `needs_curation`.
  This is how mislabelled samples are caught.
- **No profile available** (run not in Sandpiper and not analysed by MGnify): the action
  is skipped and the sample's evidence records it.
- **Calibration** can't use `marine_euks.tsv`, which is mostly eukaryotic isolates. Use
  MGnify studies with curated biome lineages, keeping in mind that Trapiche may have been
  trained on some of them (data leakage), plus the metagenomic rows of each new biome
  curated with the tool.

### 3.10 Deduplication

- The same genome in ENA and NCBI (an ENA set vs a GCA): keep one, using
  `assembly_summary_genbank.txt` and the BioSample.
- MAGs and bins from the same study: use only one layer → `same_study_mag_and_bin`.
- SPIRE bins and re-assemblies of the same source samples: run `mag-provenance` →
  `spire_bin` / `duplicate`.
- For updates: flag genomes already in the existing MGnify catalogue →
  `in_existing_catalogue`.
- **Each sample enters at its most processed public data:** genomes > contigs > reads.
  Lower-level rows of the same sample stay in the table, with
  `prepare_status = not_selected`, so they remain traceable. Nothing public is re-binned
  or re-assembled.
- **What a sample already has is checked by sample**, the way its ENA page lists it: every
  assembly analysis and WGS set linked to the sample, of any type or submitter, whose organism is
  in the catalogue's kingdom or is a metagenome (query in B.6). Searches by taxon miss data filed
  otherwise: NCBI metagenome assemblies stored as untyped WGS sets, submitters' bins, genomes
  without `metagenome_source`. A genome of the kingdom covers the sample; an unbinned assembly
  found makes its runs `not_selected` and becomes the sample's `contigs` row. A genome whose
  organism is the metagenome itself covers nothing: its kingdom is unknown, and organelle bins
  are filed this way too (e.g. PRJEB95915, the Freshwater Mitochondrial Atlas). `survey` already
  runs this check for every sample on a route. Co-assemblies sit under a sample of their own,
  linked to the originals by `derived from`; only `search` follows those links.

### 3.11 Papers that can't be read automatically

The automated path reads only open-access full text (Europe PMC). When a record still
needs its paper after the cheaper actions, and the paper isn't open access:

1. `search` sets `paper_access = needs_pdf` and puts the DOI in `publication`.
2. At the end of the run, `search` prints how many distinct papers are needed and how many
   records each would unlock.
3. The curator saves the PDFs into one folder, **with whatever name the publisher's site
   gives them**.
4. `search --papers-dir papers/` identifies each PDF and renames it to the canonical name:
   the DOI, lowercased, with `/` replaced by `_`
   (`10.1038/s41564-023-01435-6` → `10.1038_s41564-023-01435-6.pdf`).
   - **Where it looks for the DOI**, in order:
     1. the PDF metadata (XMP `prism:doi` / `dc:identifier`, Info dictionary)
     2. a DOI regex on the text of the first two pages
     3. fallback: the paper title (from Europe PMC/Crossref metadata), matched against the
        first page
   - **Only DOIs the table needs count.** A paper's first pages often cite other DOIs, so
     a match must be one of the `needs_pdf` DOIs.
   - **If exactly one DOI matches**, the file is renamed. The row gets
     `paper_access = pdf_provided`, `paper_pdf` is set to the new file name, and actions
     6 and 7 run on the extracted text.
   - **If no DOI matches, or more than one does**, the file is left untouched and listed
     at the end of the run.
   - Records with no DOI get `no_paper`.
5. **Supplements** often carry no DOI, so they can't be identified automatically. Save
   them as `<stem>.supp<N>.<ext>` (e.g. `10.1038_s41564-023-01435-6.supp1.xlsx`) next to
   the renamed paper, and `search` will read them with it.

Only papers that can still change a decision are requested. Records already in tier A,
or excluded by hard filters, are `not_needed`.

### 3.12 Reads: tracked vs used

All metagenomic WGS runs of an accepted sample are listed. What `prepare` routes now:

- **Short reads, paired or merged,** with at least one run ≥ 200k reads (preferably
  millions) → miassembler.
- **Single-end short reads** → `single_end` flag. Not routed unless they are long or
  high-coverage (A.5).
- **Long reads (ONT/PacBio)** → `long_read` flag, `prepare_status = deferred_long_read`.
  They stay in the table for tracking but aren't routed for now; genomes-generation bins
  only short reads.

A per-study summary equivalent to `filter_ena_results.py` is printed at the end of
`search`.

---

## 4. Step 2 — Workflow inputs: `prepare` and `link-assemblies`

```bash
biome-finder prepare --samples samples.tsv --biome biomes/X.yaml --tiers A,B --out <dir>
biome-finder link-assemblies --samples samples.tsv --miassembler-outdir <path> --out <dir>
```

`prepare` keeps the rows where `decision = keep`, or, where `decision` is empty, the rows
in the selected tiers. It writes the outputs below and then updates `destination`,
`local_path` and `prepare_status` in `samples.tsv`. Before rewriting the table it keeps a
copy as `samples.tsv.bak`.

```
<dir>/
├── catalogue/                           → genomes-catalogue-pipeline (one run per kingdom)
│   ├── prokaryotes/ena/                 *.fa (uncompressed) + checkM CSV written by fetch_ena.py
│   ├── prokaryotes/ncbi/                *.fa (uncompressed)
│   ├── eukaryotes/ena/  eukaryotes/ncbi/
│   ├── eukaryotes/datascout_samplesheet.csv   → datascout (protein evidence for BRAKER)
│   ├── params.prokaryotes.json          -params-file
│   └── params.eukaryotes.json
├── binning/                             → genomes-generation (one run per study)
│   └── <study>/samplesheet.csv  assembly_software.tsv  params.json
└── assembly/                            → miassembler
    └── samplesheet.csv
```

### 4.1 Genomes → genomes-catalogue-pipeline

- **ENA MAGs/bins:** `vendor/fetch_ena.py --unzip` per project (`--bins` for bins). The
  script applies the QS50 filter and writes the completeness/contamination CSV that the
  pipeline takes as `--ena_genomes_checkm`. It works per **project**, so for mixed-biome
  projects `prepare` keeps only the files of accepted rows.
- **NCBI genomes and isolates (GCA):** `vendor/fetch_ncbi.py --unzip`, into a separate
  folder (`--ncbi_genomes`). The pipeline needs ENA and NCBI genomes in separate folders
  and uncompressed.
- **Kingdoms are split** by the resolved `domain`, because the pipeline runs once per
  `--kingdom`.
- **Eukaryotes → datascout input.** `datascout_samplesheet.csv` has the columns
  `sample_id,genome,taxid,assembly_file`: the GCA accession, `reported_tax_id`,
  and the downloaded FASTA. The protein evidence that
  [datascout](https://github.com/EBI-Metagenomics/datascout) produces is then passed to
  the catalogue pipeline as `--protein_evidence`.
- `params.<kingdom>.json` holds `biome` (from `catalogue_biome`), `kingdom`,
  `ena_genomes`, `ena_genomes_checkm` and `ncbi_genomes`. The MGnify-only values (catalogue
  name, version, MGYG range) are placeholders for the operator.
- **For updates**, `prepare` also lists the old catalogue's genomes whose rows are now
  tier D or `skip`. That list is the input for `--remove_genomes`.

### 4.2 Contigs → genomes-generation

genomes-generation bins **one assembly per run**, with no co-binning, and runs **one
study at a time**, with the study-level parameters `--metagenome` and `--biomes`.

- **Samplesheet per study:** `id,assembly_fasta,fastq_1,fastq_2,assembly_accession`. The
  assembly and read files are referenced by their ENA FTP URLs, not downloaded.
- **`assembly_software.tsv`:** run → `Tool_vVersion`, taken from the ENA assembly record.
- **`params.json`:**
  - `metagenome`: the ENA metagenome taxon shared by the accepted samples
  - `biomes`: `"environment_biome,environment_feature,environment_material"`, from the
    samples' curated metadata
  - the biome config can override both
- **Only accepted rows go in.** That's why we don't use genomes-generation's own
  study-level ENA mode, which filters only by scientific name or environment_biome.

### 4.3 Reads → miassembler, then genomes-generation

No reads are downloaded. Run files are looked up in ENA (`fastq_ftp`, `fastq_md5`,
`library_layout`), and a samplesheet is written in miassembler's schema:

```csv
study_accession,reads_accession,fastq_1,fastq_2,library_layout,library_strategy,contaminant_reference
ERP123456,ERR1,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/.../ERR1_1.fastq.gz,ftp://.../ERR1_2.fastq.gz,paired,metagenomic,
```

- `study_accession` is the **secondary** study accession (ERP/SRP/DRP), as miassembler
  expects.
- Paths must end in `.fastq.gz`/`.fq.gz`. The ENA `fastq_ftp` value gets an `ftp://`
  scheme so Nextflow can stream the file.
- `contaminant_reference` comes from the biome config (host genome for host-associated
  biomes; empty for environmental ones).
- **Mixed-biome studies (A.5):** every WGS run of a study that has at least one accepted
  `reads` row is assembled, including off-biome runs and controls. Those extra runs get
  `prepare_status = assembly_only` and are never binned.
- **Cost guard:** if accepted runs / all runs in the study is below `min_study_yield`, the
  study is skipped (`low_yield_study`), unless its rows have `decision = keep`. The
  guideline example: don't assemble 1,000 runs to get 2 from the right biome.

After miassembler runs, `link-assemblies` reads its `assembled_runs.csv`,
`qc_failed_runs.csv` and per-run output tree. It writes the same per-study
genomes-generation inputs as §4.2, **for accepted runs only**, and updates `local_path` /
`prepare_status` (`queued_binning` or `assembly_qc_failed`) in `samples.tsv`.

### 4.4 Vendored scripts

The fetch scripts are copied unchanged from genomes-catalogue-pipeline at a pinned commit
(see `vendor/README.md`). They are GPL-3.0, so we run them as subprocesses and never
import them.

---

## 5. CLI and layout

```
biome-finder survey          --biomes biomes/*.yaml --out survey/                  → survey.tsv, biomes/*/, report.html
biome-finder search          --biome biomes/X.yaml [--llm] [--composition [all]] [--since] [--taxdump] [--papers-dir]  → samples.tsv
biome-finder prepare         --samples samples.tsv --biome biomes/X.yaml --out <dir>  → catalogue/, binning/, assembly/ (+ updates samples.tsv)
biome-finder link-assemblies --samples samples.tsv --miassembler-outdir <path> --out <dir>  → binning/<study>/ (+ updates samples.tsv)
```

```
biome_finder/
├── PROPOSAL.md  CLAUDE.md
├── pyproject.toml            # console_script: biome-finder
├── biome_finder/
│   ├── cli.py                # click group: entry point only, no command logic
│   ├── survey.py  search.py  prepare.py  link_assemblies.py   # one module per command
│   ├── report.py             # survey dashboard (Panel + plotly → report.html)
│   ├── models.py             # query and response dataclasses, one section per third-party service
│   ├── clients/              # one module per service: http.py (cache, retries), ena.py, mgnify.py,
│   │                         #   later ncbi.py, europepmc.py, gold.py, taxdump.py, sandpiper.py
│   ├── score.py              # biome evidence, data-type evidence, tiers
│   └── llm.py                # optional LLM extraction (imported only with --llm)
├── biomes/*.yaml             # one config per biome = the biome decision
├── benchmark/marine_euks.tsv # calibration set (plus future curated sets)
├── vendor/                   # fetch_ena.py, fetch_ncbi.py (+ helpers), IndicPiper genus_habitat_indicators CSV; unchanged
└── tests/
```

Dependencies: Python ≥ 3.11, `click`, `requests` and `PyYAML`; pandas, Panel and plotly for the
survey dashboard; plus `xmltodict` and `retry` for the
vendored scripts. Optional extras: the LLM SDK and `pypdf` (text extraction from provided
PDFs). Trapiche is an **external tool** in its own environment, called through its CLI
(`--trapiche-bin`); it is never installed into biome-finder.

---

## 6. Milestones

| # | Deliverable | Done when |
|---|---|---|
| M0 | This proposal and `CLAUDE.md` agreed | Reviewed |
| M1 | `survey` + first biome configs | Survey run on the initial targets; the user picks 1–2 biomes |
| M2 | `search` with rule-based scoring (actions 0–5, 9), data-type evidence, taxdump resolution | `samples.tsv` produced end to end for one biome; benchmark agreement reported |
| M3 | Action labels extracted from `marine_euks.tsv`; biome and data-type weights calibrated | Tier A ≥ 95% concordant; data type matches `Category` |
| M4 | LLM layer (actions 6–8) with quote verification, cache, `needs_pdf` + `--papers-dir` loop | Tier C share falls without lowering tier A precision |
| M4b | Composition evidence (action 10): IndicPiper via Sandpiper, Trapiche via subprocess | On a curated metagenomic set: tier C share falls further, and `--composition all` flags known mislabelled samples |
| M5 | `prepare`: catalogue genomes + params + datascout sheet, genomes-generation sheets for contigs, miassembler sheet, status written back | Each output passes its pipeline's input schema validation without manual edits |
| M5b | `link-assemblies` | genomes-generation runs on miassembler output for one study |
| M6 | A second biome of a different kind (e.g. host-associated) | Runs with a new YAML and no code changes, which proves the tool is biome-agnostic |

---

## 7. Risks and open questions

- **MAGs outside INSDC:** many large collections host their MAGs on Zenodo, Dryad, CNGB or
  NGDC. Those genomes can't be used directly; their samples appear as `contigs`/`reads` if
  the data is public in INSDC (ENA/NCBI/DDBJ), and get regenerated downstream.
- **Metadata quality:** many samples carry only a default organism name ("metagenome",
  "sediment metagenome"). Tier C will be large for some biomes until the paper-reading
  step is in place.
- **Data-type ambiguity:** isolate vs MAG vs SAG often needs the paper (§3.4). Records
  with `unknown` data type always go to curation.
- **Boundaries between biomes** (freshwater vs brackish, sediment vs water, rumen vs gut):
  each one needs an explicit rule in the YAML. In `marine_euks.tsv`, "coastal vs brackish"
  was the most common reason for "Needs check".
- **Long reads** are tracked but not routed until genomes-generation supports them.
- **Composition evidence coverage and limits:**
  - Sandpiper only has public runs it has already profiled, and Trapiche's taxonomy
    pathway only works for runs MGnify has analysed. Recent or non-MGnify runs may have
    no profile.
  - IndicPiper's 13 habitats are coarse, so its specificity is capped at the habitat
    level.
  - Trapiche declares no license in its repository (to confirm with the Finn lab).
    IndicPiper is GPL-3.0; we use only its database CSV, whose license on Zenodo must be
    confirmed before vendoring.
- **Pipeline schemas drift:** the input formats were read from each pipeline's
  `assets/schema_input.json` and README on 2026-11-19. Tests validate against those
  schemas, so an upstream change shows up as a test failure.

---

## 8. Suggested Claude Code skills, plugins and tools

### Already available in this environment

| Skill / plugin | Use in this project |
|---|---|
| `anthropic-skills:deep-research` | Literature sweeps per biome; checking Data-availability sections for INSDC accessions |
| `python-bioinfo-reviewer` | Review of the CLI and scoring code |
| `claude-api` | When building `llm.py`: current model IDs, structured outputs, PDF input, batch processing, prompt caching |
| `nextflow` plugin (`nextflow-config`, `nextflow-schema`, `launch-workflow`, `debug-local-run`) | Validating `prepare` outputs against the miassembler / genomes-generation / genomes-catalogue-pipeline / datascout schemas; running and debugging those pipelines |
| `claude-in-chrome` | Helping the curator collect the `needs_pdf` papers in their own logged-in browser |
| `anthropic-skills:xlsx` | Round-tripping `samples.tsv` through a spreadsheet for curators |
| `dataviz` | Survey and benchmark figures |
| `schedule` / `loop` | Periodic re-run of `survey` to detect new data for a biome |
| `code-review`, `simplify`, `ponytail` | Keep the codebase small and correct |

### Project skills (`.claude/skills/`, created)

The skills are self-contained: each one carries its own copy of the rules it needs, and
none refers to this proposal. When a rule here changes, update the copies listed below.

| Skill | Rules copied from this proposal |
|---|---|
| `define-biome` | §2.1 (YAML template), A.2, A.3, organism scope, IndicPiper habitats (§3.9) |
| `check-paper-data` | §3.2 (routes), A.1, A.3 (SPIRE), A.5 (finding the paper, read rules), A.6, B.4 |
| `curate-sample` | §3.3, §3.4, §3.6 (columns), §3.7 (actions and weights), §3.9, §3.11 (PDF names), A.1, A.4, A.5 |

| Skill | Used in | What it does |
|---|---|---|
| `define-biome` | M1 | Writes `biomes/<name>.yaml`: scope questions to the user, MGnify lineages checked against the API, ENVO/UBERON IDs verified in OLS, ENA metagenome taxa, GOLD path, IndicPiper habitat, keywords/host/data types; then runs `survey` |
| `check-paper-data` | M1 | Given a DOI or collection: Europe PMC → Data availability → accessions → ENA/NCBI checks → one route (`genomes_in_insdc` / `contigs_or_reads_in_insdc` / `not_insdc` / `unclear`), as in §3.2 |
| `curate-sample` | M2–M3 | Curates records with hard filters → data type (§3.4, isolate rules of A.4) → biome actions 1–10 (§3.7). Appends a row to `benchmark/curated_<biome>.tsv`, with `actions` (the action numbers performed) and quoted `evidence`, building the calibration set |

### Optional MCP servers

- A literature MCP server (PubMed / Europe PMC) for interactive lookups during curation,
  if one is approved in your environment. No MCP server is needed for ENA or NCBI: their
  REST APIs are called directly.

---

## 9. Data sources found for the initial targets

These are candidate studies to check during M1. Their INSDC status must be confirmed.

- **Freshwater lakes:**
  - LakePulse ([paper](https://www.nature.com/articles/s41564-023-01435-6.pdf)): reads in ENA PRJEB29238, MAGs on [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.zkh1893fs)
  - Swedish urban lakes, PRJEB54817 ([data](https://researchdata.se/en/catalogue/dataset/doi-10-17044-scilifelab-22270225))
  - Boreal Shield lakes, PRJNA1003647 ([OmicsDI](https://www.omicsdi.org/dataset/project/PRJNA1003647))
  - FRESH-MAP: a curated set of 636 freshwater metagenomes ([data](https://researchdata.se/sv/catalogue/dataset/doi-10-17044-scilifelab-28327964))
- **Rivers:** GROWdb, NCBI PRJNA946291 ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11666465))
- **Freshwater sediment:**
  - Glacier-fed stream sediments, NCBI ([paper](https://www.nature.com/articles/s41564-024-01874-9))
  - TPMC / TPMC-S ([NGDC](https://ngdc.cncb.ac.cn/tpmc-s)) and QXLSG lake sediments: **NGDC, not INSDC**, so they're only usable if mirrored
- **Permafrost:** Siberian permafrost MAGs ([PMC8130349](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8130349/))
- **Soil:**
  - SMAG (Zenodo; [paper](https://www.nature.com/articles/s41467-023-43000-z))
  - Microflora Danica long-read MAGs ([preprint](https://www.biorxiv.org/content/10.1101/2024.12.19.629313.full.pdf)): long reads, tracked only
  - High-quality soil MAG subset ([Zenodo](https://zenodo.org/records/15054613))
- **Rhizosphere:** wheat Hi-C dataset ([PMC11873117](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11873117/)); Antarctic soil/rhizosphere MAGs ([PMC12243570](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12243570/))
- **Cross-biome:** GEM catalog ([paper](https://www.doi.org/10.1038/S41587-020-0718-6)); SPIRE; Bin Chicken RBG

---

## Appendix A — Curation guidelines

This is a condensed version of the MGnify genome catalogue curation guidelines, the
process this tool automates. The code enforces these rules (see `CLAUDE.md`).

### A.1 What a genome must satisfy

- It is in **INSDC** (ENA / NCBI / DDBJ).
- It is **public**.
- It comes from the **correct biome**, as fixed in the biome config. This is not purely
  mechanical: metadata is often incomplete, and other sources (paper, project
  description, sibling samples) may be needed to confirm the biome.
- **Always assume ENA/NCBI metadata may contain errors.**
- **If you are not convinced the data is from the correct biome, leave it out.**

### A.2 Defining the biome

- **Environmental biomes:** one general catalogue or several subdivided ones (a single
  marine catalogue, or separate marine water and marine sediment; one soil catalogue, or
  one per soil type)? `survey` gives the counts; the decision is recorded in the YAML.
- **Host-associated biomes:** decide what the biome covers (cow rumen only, or the whole
  cow gut?). Environmental material near the host (e.g. barn faeces for a farm-animal gut
  catalogue) is excluded unless there is a good reason to include it.

### A.3 MAGs vs bins

- **ENA definitions:**
  - a **bin** is any genome assembled from metagenomic data
  - a **MAG** is a bin that was dereplicated across a biome within a project, so a project
    has at most one (best-quality) MAG per species per biome
- **Default to MAGs.** That keeps catalogues small, but it limits pan-genomes: one MAG
  study gives no pan-genomes, and two give at most two genomes per species.
- **Use bins** (`data_types: [bin, …]`) only when:
  - a project requires all of its data to be in the catalogue
  - the catalogue comes from a single ENA project and needs pan-genomes
  - very little data exists for the biome and no large influx is expected
- **Never use bins and MAGs from the same study**; the MAGs are a subset of the bins.
- **SPIRE** genomes sit in ENA's MAG layer but were not dereplicated before upload, so they
  are bins by ENA's definition. SPIRE also re-assembles datasets that may already be in a
  catalogue, so check them with
  [mag-provenance](https://github.com/EBI-Metagenomics/mag-provenance/tree/dev).

### A.4 Isolates

- There is no fixed process for isolates, and most catalogues have none, because it's hard
  to confirm their biome. When suitable isolates exist, they are valuable to include.
- **Finding them:**
  - search papers and check that the genomes are in INSDC, or
  - search NCBI BioProject with a biome query (`ncbi_isolate_query`, e.g.
    `bacteria[Organism] AND "human oral cavity"`). This finds both isolates and MAGs, and
    some hits are reads only (skip those).
- **The source must be confirmed for the genome itself.** A strain "usually isolated from"
  the biome, or isolated from it in the past, doesn't count. This mostly affects type
  strains with no isolation record: exclude them unless the source can be confirmed.
- **An isolate comes from a pure culture, or from the tissue of a wild-collected
  individual.** Co-cultures sequenced as metagenomes are not isolates, and SAGs are not
  isolates either. The reliable indicators are sample-level (not genome-level) metadata
  stating a pure culture or tissue sampling, or the paper or study description saying so.
- **Marine water collected at the coast counts as ocean**
  (`root:Environmental:Aquatic:Marine:Oceanic`). If it's seawater, it's ocean, even at the
  coast. Coastal sediment, intertidal, estuarine and brackish environments keep their own
  lineages. A marine-water biome config encodes this through `lineage_equivalents`.
- **Don't trust the word "isolate".** When in doubt, classify the genome as a MAG.

### A.5 Curating studies and samples

- **Find the publication** that generated the data:
  1. ENA "Publications in which the generation of this record is described"
  2. if that's empty, "Publications citing this record". Be careful: this tab also lists
     papers that merely cite the study.
  3. search the web for the PRJ accession and the ERP/SRP/DRP accession
  4. if there's no paper, the ENA study description is often enough
- **Check in the paper:**
  - Are any samples controls?
  - Are there samples from other biomes?
  - Is it really WGS?
  - Were samples manipulated in a way we don't want (e.g. spiked with microbes)?
- **Raw-read studies:**
  - Check that MAGs or assemblies don't already exist, using ENA "Related ENA Records",
    both for the study and for one of its samples.
  - Avoid single-end short reads, which yield few MAGs. Exceptions: paired reads merged
    before upload, and high-coverage or long-read runs.
  - Skip studies with no run of reasonable size. 200k reads is the absolute minimum;
    millions are better. It's fine if only some runs pass, since the assembly pipeline
    drops bad runs.
  - Mixed-biome studies are fine, but **assemble all runs**. Only runs from the biome go
    on to MAG generation. Weigh cost against value: if 1,000 runs must be assembled to
    get 2 from the biome, leave the study out (`min_study_yield`).
  - Studies mixing amplicon and WGS are fine; the assembly pipeline ignores amplicon runs.
- **Controls and synthetic metagenomes** are assembled but never binned.
- **Genomes we didn't generate** must be traceable to the original biological sample,
  usually through the `derived from` field. This is not always straightforward.

### A.6 Where to look for data

- MAG collection lists: [awesome-mag](https://github.com/inspirewind/awesome-mag) (§3.2).
- [Metalog](https://metalog.embl.de/): curated sample metadata. It doesn't cover every
  environment or only INSDC studies, but it saves curation time.
- [MAGdb](https://magdb.nanhulab.ac.cn/): curated, QC'ed data (currently offline).
- Papers describing MAG collections (search "metagenome assembled genomes {biome}").
  Their genomes are often outside INSDC (check the Data-availability section). If the
  reads or assemblies are public but the MAGs never were, generate the MAGs ourselves.
- The ENA MAG layer and bin layer, and metagenomic read runs (Appendix B).
- NCBI (optional; most catalogues use only ENA data). NCBI-deposited genomes don't appear
  in ENA's MAG layer; they show up as GCA accessions elsewhere in ENA (B.4).

---

## Appendix B — Reference queries

`clients/` reproduces these queries programmatically. They're kept here as the
reference for tests and for manual checks.

### B.1 ENA MAG layer

```bash
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d "result=wgs_set&query=assembly_type%3D%22metagenome-assembled%20genome%20(mag)%22&fields=accession%2Cassembly_type%2Cstudy_accession%2Csample_accession%2Cmetagenome_source%2Cisolation_source%2Cfirst_public%2Cassembly_accession%2Cset_fasta_ftp&format=tsv" \
  "https://www.ebi.ac.uk/ena/portal/api/search" > results_mags.tsv
```

### B.2 ENA bin layer (only when bins are enabled)

```bash
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" \
  -d "result=analysis&query=assembly_type%3D%22binned%20metagenome%22&fields=analysis_accession%2Cassembly_type%2Cstudy_accession%2Csample_accession%2Cenvironment_biome%2Cstudy_title%2Cisolation_source%2Clocal_environmental_context%2Cenvironment_material%2Cgenerated_ftp&format=tsv" \
  "https://www.ebi.ac.uk/ena/portal/api/search" > results_bins.tsv
```

The same `analysis` query with `assembly_type="primary metagenome"` returns the unbinned
metagenome assemblies (`contigs`).

### B.3 Metagenomic read runs (WGS or OTHER)

```bash
curl -X GET "https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=ena&fields=study_accession%2Cstudy_title%2Cproject_name%2Clibrary_strategy%2Clibrary_source%2Crun_accession%2Cread_count%2Cfirst_public&includeMetagenomes=true&limit=0&result=read_run&query=library_source%3D%22METAGENOMIC%22%20AND%20(library_strategy%3D%22WGS%22%20OR%20library_strategy%3D%22OTHER%22)" \
  -H "accept: */*" > ena_request.tsv
```

**Per-study summary** (the logic of the guidelines' `filter_ena_results.py`, which
`search` reproduces; no separate script):

- drop runs with `first_public` before `--since` (runs with an unparseable date are kept
  as "Unknown")
- drop `AMPLICON` runs and anything not `METAGENOMIC`
- group by `study_accession` and report: study title, project name, number of runs,
  min/max/mean `read_count`, and study date

### B.4 NCBI

- **WGS genomes:** the guidelines use the
  [Sequence Set Browser](https://www.ncbi.nlm.nih.gov/Traces/wgs/), with Type = `wgs`
  (excludes TSA/TLS), the Host facet, and a Term such as "gut" matched against the
  isolation source. `search` does the same through E-utilities with `ncbi_wgs_terms`.
- **GCA accessions:** map the BioSample IDs to GCA accessions in
  [`assembly_summary_genbank.txt`](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt).
  Download a fresh copy every run. Map **by sample, not by project**, because projects can
  span several biomes.
- Never duplicate genomes already found in ENA.

### B.5 Downloads

These are the vendored `fetch_ena.py` (by ENA project; QS50 pre-filter; completeness and
contamination CSV) and `fetch_ncbi.py` (by GCA list). See §4.1.

### B.6 Data already linked to a sample

```bash
# 200 samples per query; the sample list and the clause are filled in by the code
curl -X POST https://www.ebi.ac.uk/ena/portal/api/search -d result=analysis \
  --data-urlencode 'query=analysis_type=SEQUENCE_ASSEMBLY AND (sample_accession="SAMEA…" OR …) AND (tax_tree(2) OR tax_tree(2157) OR tax_tree(408169)) AND NOT tax_tree(10239) …' \
  -d fields=analysis_accession,sample_accession,assembly_type,tax_id,scientific_name,study_accession -d format=tsv
# the same with result=wgs_set and fields=accession,sample_accession,assembly_type,tax_id,scientific_name,study_accession
```

The organism must be in the catalogue's kingdom (here prokaryotes: 2, 2157) or a metagenome
(408169), minus the excluded taxa, so a bacterial MAG doesn't cover a sample for a eukaryote
catalogue. An untyped WGS set whose organism is a metagenome (`soil metagenome`) is an unbinned
assembly; a genome with a metagenome as organism has no known kingdom and covers nothing.

---

## Appendix C — MGnify catalogues snapshot (2026-11-19)

This is the baseline from [MGnify Genomes](https://www.ebi.ac.uk/metagenomics/browse/genomes)
on the date this proposal was written. `survey` reads the live list from the MGnify API;
this table is for reference and tests.

| catalogue_id | version | updated_at | species (genome_count) | genomes (unclustered) |
|---|---|---|---|---|
| zebrafish-fecal-v1-0 | 1.0 | 2022-11-28 | 79 | 101 |
| human-gut-v3-0 | 3.0 | 2026-09-29 | 6,461 | 286,621 |
| human-skin-v1-0 | 1.0 | 2025-07-25 | 579 | 4,611 |
| human-vaginal-v1-0 | 1.0 | 2023-08-15 | 280 | 618 |
| maize-rhizosphere-v1-0 | 1.0 | 2025-08-29 | 336 | 431 |
| marine-v2-0 | 2.0 | 2024-05-12 | 13,223 | 50,866 |
| marine-sediment-v1-0 | 1.0 | 2025-11-19 | 6,158 | 10,126 |
| mouse-gut-v1-0 | 1.0 | 2023-12-08 | 2,847 | 112,951 |
| pig-gut-v1-0 | 1.0 | 2022-11-29 | 1,376 | 3,972 |
| cow-rumen-v1-0-1 | 1.0.1 | 2023-01-23 | 2,729 | 5,578 |
| sheep-rumen-v1-0 | 1.0 | 2024-08-16 | 2,172 | 2,904 |
| soil-v1-0 | 1.0 | 2025-07-16 | 19,472 | 20,908 |
| tomato-rhizosphere-v1-0 | 1.0 | 2025-08-29 | 579 | 666 |
| human-oral-v1-0-1 | 1.0.1 | 2024-01-23 | 452 | 1,225 |
| barley-rhizosphere-v2-0 | 2.0 | 2025-12-12 | 86 | 103 |
| honeybee-gut-v1-0-1 | 1.0.1 | 2024-01-23 | 154 | 627 |
| non-model-fish-gut-v2-0 | 2.0 | 2023-04-26 | 178 | 254 |
| chicken-gut-v1-0-1 | 1.0.1 | 2024-01-23 | 1,322 | 13,386 |
| marine-eukaryotes-v1-0 | 1.0 | 2026-10-27 | 74 | 385 |
| human-skin-eukaryotes-v1-0 | 1.0 | 2026-10-28 | 11 | 288 |
