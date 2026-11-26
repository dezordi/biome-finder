---
name: curate-sample
description: Manually curate one or more INSDC records (sample, genome, assembly or run accessions) for a biome-finder target biome — decide the data type (MAG / bin / isolate / SAG / contigs / reads / unknown) and the biome, following the ordered curation actions, then append the result to benchmark/curated_<biome>.tsv. Use when the user asks to curate, check, verify or decide on specific accessions, rows from samples.tsv (tier C / needs_curation), or to build benchmark data.
---

# curate-sample

## Core rules

- A genome must be **in INSDC** (ENA / NCBI / DDBJ), **public**, and from the **correct
  biome** as fixed in `biomes/<name>.yaml`.
- **The sample (BioSample) is the unit of biome curation, not the project.** Mixed-biome
  projects are common.
- **A co-assembled genome is accepted only if all its source samples are accepted.**
- **Always assume ENA/NCBI metadata may contain errors.** Confirming the biome often needs
  other sources: the paper, the project description, sibling samples.
- **If you are not convinced the data is from the correct biome, leave it out**, and say
  why.
- **Genomes we didn't generate must be traceable to the original biological sample**,
  usually via `derived from`. If you can't trace it, `skip`.
- **Organism scope:** Bacteria, Archaea and microbial Eukaryotes. Out of scope:
  - Viruses (10239), Metazoa (33208), Tracheophyta (58023), Phaeophyceae (2870, brown
    algae and kelps), plus the biome's `exclude_taxa_extra`
  - **macroalgae in general**: multicellular, macroscopic seaweeds. Red and green algae
    have no single taxon to exclude, so judge them per genus (e.g. *Ulva*, *Porphyra*,
    *Gracilaria* are out; unicellular or microscopic species are in)
  - plasmids and other MGEs
  - organelle genomes

## Inputs

- One or more accessions, or rows from `samples.tsv`.
- The target biome config, `biomes/<name>.yaml`. Ask for it if it isn't obvious. Its
  `kingdom`, `mgnify_lineages` minus `mgnify_lineages_exclude`, `envo_include` minus
  `envo_exclude` (all as whole subtrees; an excluded term in any `env_*` field rules the
  sample out), `keywords_*`, `host` and `notes` define what counts as the biome.

## Step 1: Resolve the record

Find the sample, the source samples (`derived from`, for genomes), the study (primary and
secondary accessions), the reported taxon, and the linked publications.

- `https://www.ebi.ac.uk/ena/browser/api/xml/<acc>`:
  - for a sample: its attributes and links
  - for an assembly: `DESCRIPTION`, `ASSEMBLY_LEVEL`, `SAMPLE_REF`, `STUDY_REF`
- `https://www.ebi.ac.uk/ena/xref/rest/json/search?accession=<acc>`: related records and
  publications (`&source=EuropePMC` for papers only)
- `https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<GCA>/dataset_report`: the
  best single source for a genome. It gives the BioSample package, every BioSample
  attribute (`culture_collection`, `strain`, `isolate`, `host`, `env_*`,
  `isol_growth_condt`, `mixed sample`), the assembly status and the assembly level.
- ENA Portal `assembly` search: `accession` is **versionless**. To query a versioned GCA,
  use `assembly_set_accession`, and quote the values:
  ```bash
  curl -s -G https://www.ebi.ac.uk/ena/portal/api/search --data-urlencode result=assembly \
    --data-urlencode 'query=assembly_set_accession="GCA_900128745.1" OR assembly_set_accession="GCA_016586335.1"' \
    --data-urlencode fields=accession,assembly_type,assembly_level,sample_accession --data-urlencode format=tsv
  ```

## Step 2: Hard filters (action 0)

Each failure is a `skip`; record the reason.

- Not public in INSDC, or no data available (e.g. no genome assembly, only reads, when a
  genome was expected).
- Out-of-scope taxon (see Core rules), or an organelle genome.
- Control, mock or synthetic metagenome. In mixed studies these are *assembled* but never
  binned.
- Spiked sample (organisms added), or otherwise manipulated in a way we don't want. An
  enrichment culture is fine if no organisms were added while culturing: check what was done
  (sample metadata, paper methods) before excluding it.
- Amplicon-only data (not WGS).
- Off-biome material, e.g. a host conflict for an environmental biome, or barn faeces for a
  farm-animal gut biome.

## Step 3: Data type

Check the evidence in this order:

| Field / source | Indicates |
|---|---|
| ENA `assembly_type` | `Metagenome-Assembled Genome (MAG)` → mag · `binned metagenome` → bin · `primary metagenome` → contigs · `clone or isolate` → isolate (still needs pure-culture or tissue evidence). Often empty for isolates and SAGs |
| ENA assembly `DESCRIPTION` | "cell sorted by flow cytometer", "amplified" (MDA/WGA) → sag · "binned", "metagenome" → mag/bin · "strain", "culture" → isolate candidate |
| BioSample `mixed sample` / `member sample` | `in silico mixture` → pooled source samples (common for Tara SAGs/MAGs) |
| NCBI `excluded_from_refseq` | `derived from metagenome` → mag/bin · `derived from single cell` → sag |
| BioSample package | `MIMAG.*` → mag · `MISAG.*` → sag · `MIGS.ba/eu` → candidate isolate |
| BioSample attributes | `metagenome-source`, `derived_from` → mag. A source-sample organism under "metagenomes" (taxid 408169) → metagenome-derived |
| Strain / culture collection IDs (ATCC, NIES, CCAP, DSM…) | Supports isolate; never proves the isolation source |
| Paper Methods | "pure culture", "axenic", "single-cell sorted", "binned", "co-culture" |

**Isolate rules**

- **Don't trust the word "isolate."**
- `isolate` means a single-organism genome sequenced from either:
  - a **pure culture**: the reliable indicators are sample-level (not genome-level)
    metadata stating a pure culture, or the paper or study description saying so; or
  - **tissue or whole individuals of a wild-collected organism**, not a metagenome
    (e.g. BioSample `tissue`, `dev_stage`, or a paper saying "wild-collected
    individuals"). The organism must still be in scope; for example, macroalgae are out.
- Co-cultures sequenced as metagenomes are not isolates. SAGs are not isolates either:
  `sag`, excluded by default.
- **The isolation source must be confirmed for this genome.** A strain "usually isolated
  from" the biome, or isolated from it in the past, doesn't count. This mostly affects
  type strains with no isolation record: `skip` unless the source is confirmed.
- **When in doubt, classify the genome as `mag`** if it is metagenome-derived. Otherwise
  use `unknown` and set `needs_curation = TRUE`.

**Other data-type rules**

- `contigs` means an unbinned **metagenome** assembly. A single-organism genome at contig or
  scaffold level is still `mag` or `isolate`.
- **MAG vs bin:** a MAG is a bin dereplicated across a biome within a project. Never accept
  bins and MAGs from the same study.
- **SPIRE** genomes in ENA's MAG layer are bins (they were never dereplicated). Flag them.
- **Reads** (`reads` records):
  - Single-end short reads yield few MAGs. Accept them only if they are merged pairs or
    high-coverage runs.
  - Long reads are tracked but not routed for now. Note this.
  - At least one run needs ≥ 200k reads (millions preferred).
  - First check that MAGs or assemblies don't already exist ("Related ENA Records" of the
    study and of a sample). If they do, curate those instead.

## Step 4: Biome actions

Run the actions cheapest first and stop once the evidence is conclusive. Record every
action you perform, even when it finds nothing. Every claim needs a **verbatim quote** and
its URL.

| # | Action | Typical weight |
|---|---|---|
| 1 | Third-party curated biome: MGnify study lineage (`https://www.ebi.ac.uk/metagenomics/api/v1/studies?accession=<ERP/PRJ>`), GOLD ecosystem, Metalog | 0.9 |
| 2 | Structured sample metadata: ENVO `env_broad_scale` / `env_local_scale` / `env_medium`, `env_biome`, `isolation_source`, host, body site. **ENA `env_tax:*` / `env_geo:*` tags don't count**: ENA derives them from the organism's taxonomy, not from the sample (one sample can carry both `env_tax:marine` and `env_tax:freshwater`) | 0.8 when specific, 0.4 when generic ("marine", "water"), ~0 for `env_tax`/`env_geo` tags |
| 3 | Consistency with sibling samples of the same project | ± 0.1 |
| 4 | Project/study title and description. This includes MGnify-style "This project includes samples from the following biomes: root:…" | 0.6 |
| 5 | Metagenome organism name (e.g. "freshwater metagenome"). It's often a default value | 0.3 |
| 6 | Paper Methods / Data availability (see below) | 0.8 if explicit, 0.3 if only the title or a region |
| 7 | Supplementary tables: sample sites, coordinates | 0.8 |
| 8 | Culture collection record (isolates only): ATCC, NIES, CCAP, DSM, Chlamydomonas collection… | 0.7 |
| 9 | Coordinates sanity check (e.g. a "freshwater" sample in the open ocean) | conflict only |
| 10 | Community composition (see below) | Trapiche 0.6 · IndicPiper 0.4 |

**Combining:**

- `specificity` is how deep a claim reaches into the target lineage. "Aquatic" only
  partially supports "Aquatic:Freshwater:Sediment".
- Agreeing items reinforce each other.
- Any conflicting item with weight ≥ 0.6 means `needs_curation = TRUE`.
- **Composition evidence alone is never enough to `keep`.** Neither is evidence that only
  an LLM inferred without a verbatim quote.

**How the benchmark curator worked** (`marine_euks.tsv`, 848 notes):

- Reading the paper was the most common action, but it was conclusive only about 46% of
  the time. Many papers "only report [the source] in the title" or "do not inform the
  source".
- Supplementary tables and specific sample metadata were usually conclusive.
- Culture-collection pages were conclusive about 43% of the time.
- The most common open question was a boundary case ("coastal vs brackish", "glacier vs
  freshwater").

**Finding and reading the paper (action 6):**

1. Look under ENA "Publications in which the generation of this record is described".
2. If that's empty, look under "Publications citing this record", but check that the paper
   actually generated the data.
3. Search the web for the PRJ accession and for the ERP/SRP/DRP accession.
4. If there's no paper, the ENA study description is often enough.

Open-access full text: `https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML`.
If the paper isn't open access, ask the user for the PDF. Its canonical name is the DOI,
lowercased, with `/` replaced by `_`, plus `.pdf`
(`10.1038/s41564-023-01435-6` → `10.1038_s41564-023-01435-6.pdf`). Supplements are
named `<stem>.supp<N>.<ext>`.

In the paper, check:

- Are any samples controls?
- Are there samples from other biomes?
- Is it really WGS?
- Were samples manipulated (e.g. spiked with microbes)?

**Community composition (action 10).** Use it only when the evidence so far is unclear, or
to check a suspected mislabel. It needs public reads: the sample's own runs, or the source
sample's runs for a MAG. It doesn't apply to isolates.

- **IndicPiper:**
  1. Fetch the run's GTDB genus profile:
     `https://sandpiper.qut.edu.au/api/condensed_csv_with_extras/<run>?taxonomy_type=gtdb`.
  2. Keep the genus rows and strip the `d__`…`g__` prefixes.
  3. Exact-match the genera against the IndicPiper database (v2.1.0, 767 genera, 13
     habitats). Repository: https://github.com/cliffbueno/IndicPiper, file
     `genus_habitat_indicators_v2.1.0.csv`, columns `Taxonomy, Habitat, …`.
  4. Sum the relative abundance per habitat.

  The evidence is whether the biome's `indicpiper_habitats` dominate the indicator signal.
  IndicPiper can separate freshwater sediment from freshwater water, but not lake sediment
  from river sediment.
- **Trapiche** (https://github.com/Finn-Lab/trapiche): only if it is installed and the run
  has MGnify taxonomy files (SSU/LSU OTU). Pass the text-based lineage as
  `ext_text_pred_sample`. Use `final_selected_prediction`.
- If composition contradicts strong metadata, flag `composition_conflict` and set
  `needs_curation = TRUE`. If no profile exists, record that and move on.

## Step 5: Decide

- `curated_biome`: the deepest lineage the evidence supports, after applying the biome
  config's `lineage_equivalents`. **Marine water collected at the coast is
  `root:Environmental:Aquatic:Marine:Oceanic`**, not `…:Marine:Coastal`. If it's
  seawater, it's ocean, even at the coast (e.g. the Gulf of Naples, the Salish Sea).
  Coastal *sediment*, intertidal, estuarine and brackish environments keep their own
  lineages.
- `decision`:
  - `keep`: the evidence clearly supports the target biome.
  - `skip`: the record failed a hard filter, or the evidence is too generic, conflicting
    or off-biome.
- `needs_curation = TRUE` when the decision depends on a judgement the user should
  confirm. Write the open question in the note.
- `curation_note`: one or two sentences in the style of the benchmark's "Manual Inspection"
  column. Examples:
  - `The paper (https://…) confirms the source as "cores were taken from Lake X"`
  - `The sample informs the source as soil (https://www.ncbi.nlm.nih.gov/biosample/SAMD…)`
  - `No detailed information about the source in ENA or NCBI, and no publication`
  - `The paper (https://…) only reports in the title that the samples are from <region>`

## Step 6: Write `benchmark/curated_<biome name>.tsv`

- Append one row per record. If the file doesn't exist, create it with the header.
- If a row with the same `record_accession` already exists, **show the old and new rows
  and ask** before replacing it.
- Write with Python's `csv` module (tab-delimited, `csv.QUOTE_MINIMAL`). Never hand-edit
  the file with sed.

The columns have a fixed order. New columns may be added only at the end:

```
record_accession  data_type  data_type_evidence  sample_accession  study_accession  publication
reported_tax_name  reported_biome  curated_biome  decision  needs_curation
actions  evidence  curation_note  curated_by  curated_on
```

| Column | Content |
|---|---|
| `record_accession` | Genome (ENA set / GCA), ERZ assembly, or run accession |
| `data_type` | `mag` · `bin` · `isolate` · `sag` · `contigs` · `reads` · `unknown` |
| `data_type_evidence` | e.g. `ena.assembly_type=clone or isolate; biosample.package=MIGS.eu; paper:"axenic culture"` |
| `sample_accession` | BioSample; for co-assemblies, the source samples separated by `;` |
| `study_accession` | PRJ* (add the secondary ERP/SRP/DRP after a `;` if known) |
| `publication` | DOI / PMID, separated by `;` |
| `reported_tax_name` | Taxon name as shown in the record |
| `reported_biome` | Raw metadata as found: `env_*`, `isolation_source`, organism, host, body site |
| `curated_biome` | Your lineage, e.g. `root:Environmental:Aquatic:Freshwater:Sediment` |
| `decision` | `keep` / `skip` |
| `needs_curation` | `TRUE` / `FALSE` |
| `actions` | The action numbers performed, in order, e.g. `0;2;4;6` |
| `evidence` | `action:source:"quote"` items separated by `;` |
| `curation_note` | Manual Inspection-style sentence |
| `curated_by` | `claude` (the user replaces it after review) |
| `curated_on` | Today, YYYY-MM-DD |

## Output to the user

- One line per record: accession → data type → decision → curated biome, plus the note.
- The records that need curation, with their open questions.
- The PDFs needed, by canonical file name.
