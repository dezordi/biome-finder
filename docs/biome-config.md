# The biome config

Each catalogue is defined by one YAML file, `biomes/<name>.yaml`. The file is both the tool's
input and the written record of the biome decision: what goes into the catalogue, what
stays out, and why. The code has no biome-specific logic, so every query `survey`, `search`
and `prepare` sends is shaped by this file.

Like MGnify, there is one config per biome **and kingdom** (`soil-prok`, `soil-euk`), and
host-associated biomes get one config per host (barley, wheat and rice rhizosphere are
separate). The [define-biome skill](#adding-or-changing-terms) writes and revises configs
with these rules.

## Which step uses each field

`survey` is implemented; `search` and `prepare` are planned, so their columns describe the
design in `PROPOSAL.md`.

| Field | survey | search | prepare |
|---|---|---|---|
| [`name`](#name) | output names | output names | output names |
| [`kingdom`](#kingdom) | MAG search, kingdom evidence, catalogue check | taxon scope | `--kingdom`, folder |
| [`mgnify_catalogue`](#mgnify-catalogue) | catalogue context, cutoff, genomes already in it | flags genomes already in it | — |
| [`catalogue_biome`](#catalogue-biome) | — | target lineage of the score | `--biome` |
| [`mgnify_lineages`](#mgnify-lineages) | — | biome evidence (action 1) | — |
| [`mgnify_lineages_exclude`](#mgnify-lineages-exclude) | — | conflict evidence | — |
| [`lineage_equivalents`](#lineage-equivalents) | — | maps evidence lineages | — |
| [`envo_include`](#envo-include-and-envo-exclude) | — | biome evidence (action 2) | — |
| [`envo_exclude`](#envo-include-and-envo-exclude) | — | rules a sample out | — |
| [`gold_ecosystems`](#gold-ecosystems) | — | biome evidence (action 1) | — |
| [`metagenome_taxa`](#metagenome-taxa) | **the ENA searches** | candidate records, weak evidence (action 5) | — |
| [`keywords_include`](#keywords-include-and-keywords-exclude) | — | text evidence (actions 2, 4) | — |
| [`keywords_exclude`](#keywords-include-and-keywords-exclude) | — | conflicts, hard filters | — |
| [`host`](#host-host-tax-id-host-keywords) | — | host conflict check | — |
| [`host_tax_id`, `host_keywords`](#host-host-tax-id-host-keywords) | host filter | host evidence | — |
| [`body_site_include`](#body-site-include) | — | body-site evidence | — |
| [`contaminant_reference`](#contaminant-reference) | — | — | miassembler samplesheet |
| [`data_types`](#data-types) | — | which record types to collect | routing |
| [`exclude_taxa_extra`](#exclude-taxa-extra) | MAG search | out-of-scope filter | — |
| [`ncbi_isolate_query`](#ncbi-isolate-query-and-ncbi-wgs-terms) | — | NCBI isolate search | — |
| [`ncbi_wgs_terms`](#ncbi-isolate-query-and-ncbi-wgs-terms) | — | NCBI WGS search | — |
| [`min_study_yield`](#min-study-yield) | — | `low_yield_study` | skips studies |
| [`indicpiper_habitats`](#indicpiper-habitats) | — | composition evidence (action 10) | — |
| [`notes`](#notes) | — | — | — |

## Field reference

The example values come from `biomes/marine-water-euk.yaml` and
`biomes/barley-rhizosphere-prok.yaml`. Fields stay in the order of the template. Each field
says where to look up new values. Every ID goes in verified, never from memory; when a value
can't be confirmed, write it with a `# TODO verify` comment.

### Identity

#### name

The config's name, `<biome>-prok` or `<biome>-euk`. It names the survey's `biome` column,
the `survey/biomes/<name>/` folder and the catalogue in the dashboard. Keep it equal to the
file name.

#### kingdom

**Required.** `prokaryotes` (Bacteria + Archaea) or `eukaryotes`, the type of MGnify
catalogue and the genomes-catalogue-pipeline `--kingdom`.

- **survey:** MAGs and MGnify bins are searched under Bacteria (2) and Archaea (2157), or
  under Eukaryota (2759). Assemblies and runs have no organism, so the kingdom decides
  which of them count: a eukaryote config needs evidence (a eukaryote genome from the same
  sample, or a project that names eukaryotes); a prokaryote config takes everything but
  eukaryote-only projects. The survey stops if `mgnify_catalogue` is of the other type.
- **search / prepare:** the taxon scope of genomes, and one catalogue folder and one
  pipeline run per kingdom.

#### notes

Free text: the boundary decisions in plain words (what is in, what is out, and why). Not
read by the code; it is the record for the next curator.

### Catalogue

#### mgnify_catalogue

The MGnify catalogue this config would update, by its API v2 id (`soil-v1-0`,
`marine-eukaryotes-v1-0`), or `""` for a new catalogue. Its `catalogue_type` must equal
`kingdom`.

**Where to look:** [MGnify genome catalogues](https://www.ebi.ac.uk/metagenomics/browse/genomes).
Use the v2 API; v1 lists outdated versions.

```bash
curl -s "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/catalogues/"
```

- **survey:** species and genome counts, the last update (the cutoff of every `*_new`
  column) and the genome list on the MGnify FTP, used to tell which genomes are already in
  the catalogue (`*_to_add`).
- **search:** flags genomes already in the catalogue.

#### catalogue_biome

Exactly one MGnify lineage, passed to genomes-catalogue-pipeline as `--biome`. It is usually
the top of the scope, e.g. `root:Environmental:Aquatic:Marine`. `search` scores each record
against it.

**Where to look:** [MGnify biomes](https://www.ebi.ac.uk/metagenomics/browse/biomes). The
lineage must exist there:

```bash
curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/root:Environmental:Aquatic:Marine"
```

### Biome scope

The scope is always **subtrees minus excluded subtrees**: a lineage or ontology term brings
in all its descendants, and an excluded one removes its own subtree. This keeps water and
sediment, or soil and permafrost, apart.

#### mgnify_lineages

MGnify biome lineages that count as the biome, each with its whole subtree.

**Where to look:** [MGnify biomes](https://www.ebi.ac.uk/metagenomics/browse/biomes). Read
the full children listing: MGnify often files the same environment under several parents
(freshwater sediment sits under `Freshwater:Sediment`, `Lentic:Sediment`, `Lotic:Sediment`
and `Pond:Sediment`), so list each of them.

```bash
curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/root:Environmental:Aquatic:Marine"   # exists? samples-count
curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/root:Environmental:Aquatic:Marine/children?page_size=100"
```

Used by `search`: a study's MGnify lineage is third-party curated biome evidence (action 1).

#### mgnify_lineages_exclude

Descendants of `mgnify_lineages` that belong to another biome; each removes its whole
subtree. Example: marine-water-euk lists `root:Environmental:Aquatic:Marine` and excludes
its 11 sediment subtrees (`…:Marine:Sediment`, `…:Marine:Coastal:Sediment`, …). In
`search`, a lineage in an excluded subtree is evidence against the biome.

**Where to look:** the children listing of each `mgnify_lineages` entry (above); exclude
every child that belongs to another biome.

#### lineage_equivalents

A map from an evidence lineage to the lineage it counts as, for boundary rules that the
lineage tree alone gets wrong. Example: marine water collected at the coast counts as ocean,
so a marine-water config maps `root:Environmental:Aquatic:Marine:Coastal` to
`root:Environmental:Aquatic:Marine:Oceanic`, while coastal sediment, intertidal, estuarine
and brackish stay separate. Used by `search` when it compares evidence with
`catalogue_biome`.

**Where to look:** both sides are lineages from
[MGnify biomes](https://www.ebi.ac.uk/metagenomics/browse/biomes), listed as for
`mgnify_lineages`. The rule itself comes from the curation guidelines (`PROPOSAL.md`,
Appendix A), not from a database.

#### envo_include and envo_exclude

ENVO terms (`ENVO:00002150`, coastal sea water), each meaning the term plus all its
subclasses.

**Where to look:** [ENVO in OLS](https://www.ebi.ac.uk/ols4/ontologies/envo). Search by name
and check the hierarchy before choosing a term:

```bash
curl -s "https://www.ebi.ac.uk/ols4/api/search?q=coastal%20sea%20water&ontology=envo&rows=10"   # obo_id + label
```

To see which terms samples of the biome really use, look at the ENA sample fields
`broad_scale_environmental_context`, `local_environmental_context` and
`environmental_medium` (or the older `environment_biome`, `environment_feature`,
`environment_material`).

- **search:** a sample whose `env_*` fields fall under `envo_include` gets structured
  metadata evidence (action 2). A sample with **any** `env_*` field under `envo_exclude` is
  ruled out, even when another field matches: broad scale "marine biome" plus medium
  "sediment" is out of a marine-water biome. Use `envo_exclude` for the same splits as
  `mgnify_lineages_exclude`.
- ENA `env_tax:*` and `env_geo:*` sample tags come from the organism's taxonomy, not from
  the sample, so they are never biome evidence.

#### gold_ecosystems

GOLD ecosystem paths, five levels from Ecosystem to Specific Ecosystem, e.g.
`Environmental > Aquatic > Freshwater > Lake > Sediment`. In `search`, a GOLD classification
of the sample is curated biome evidence (action 1); awesome-mag seed collections use the
same paths.

**Where to look:** browse the five levels of the
[GOLD ecosystem classification](https://gold.jgi.doe.gov/ecosystem_classification), also on
[BioPortal (GOLDTERMS)](https://bioportal.bioontology.org/ontologies/GOLDTERMS). Unconfirmed
paths carry `# TODO verify`.

#### keywords_include and keywords_exclude

Lowercase phrases matched in free text: sample fields (`isolation_source`, descriptions)
and the project title and description.

- `keywords_include`: specific phrases for the biome ("seawater", "coastal water",
  "rhizosphere"). Avoid words that also describe neighbours.
- `keywords_exclude`: the boundary neighbours ("estuar", "brackish", "sediment" for marine
  water; "bulk soil", "endosphere" for rhizosphere) plus manipulations that are always out
  ("enrichment", "mesocosm", "microcosm", "spiked", "mock", "synthetic", "control").
  Prefixes work: "estuar" matches estuary and estuarine.

Used by `search` only (actions 2 and 4).

**Where to look:** the study titles in `survey/biomes/<name>/studies.tsv` after a survey run,
and the reasons in the curated benchmark tables (`benchmark/`).

#### indicpiper_habitats

The IndicPiper habitats that count as this biome, chosen from its 13 (database v2.1.0):
soil (includes rhizosphere), marine sediment, marine water, mammalian gut, freshwater
sediment, freshwater water, human oral, wastewater, glacier, hydrothermal vent,
phyllosphere, human vaginal, human skin. Use `[]` if none fits. Used by `search
--composition`: IndicPiper sums the indicator genera of the run's Sandpiper profile per
habitat, and the target habitat dominating is evidence for the biome.

**Where to look:** only this fixed list; the habitats come with the pinned IndicPiper
database.

### Organisms and data

#### metagenome_taxa

ENA/NCBI organism names under **metagenomes** (taxid 408169), e.g. `marine metagenome`,
`seawater metagenome`. This is the field that matters most for `survey`.

**Where to look:** the children of metagenomes (408169) in
[ENA](https://www.ebi.ac.uk/ena/browser/view/Taxon:408169) or
[NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=408169), or
search by words:

```bash
curl -s "https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-submission/seawater%20metagenome"
```

Keep entries with `"metagenome": "true"` that fit the scope. Generic names
(`metagenome`, `sediment metagenome`, `soil metagenome` for a rhizosphere) pull in off-biome
samples; include them only on purpose, and say so in `notes`.

- **survey:** each name is resolved to a tax ID and gets its own searches and output row.
  MAGs match it through `metagenome_source`; assemblies and runs through the sample's
  taxon, exactly (`tax_eq`, no child taxa).
- **search:** the starting set of candidate records, and weak biome evidence (action 5):
  submitters often keep a default name.

#### exclude_taxa_extra

NCBI tax IDs added to the default exclusions: Viruses (10239), Metazoa (33208),
Tracheophyta (58023) and Phaeophyceae (2870). Each removes its whole tax tree. `survey`
removes them from the MAG search; `search` from every genome.

**Where to look:** [ENA](https://www.ebi.ac.uk/ena/browser/) or
[NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), or by scientific name:

```bash
curl -s "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/Rhodophyta"
```

#### data_types

The record types to collect: `mag`, `isolate`, `contigs`, `reads`, plus `bin` and `sag`
only when wanted. Default: `[mag, isolate, contigs, reads]`. Enable `bin` only when a
project needs all its data in the catalogue, the catalogue comes from one ENA project and
needs pan-genomes, or the biome has very little data (`PROPOSAL.md`, Appendix A.3). Used by
`search` (what to collect) and `prepare` (routing). `survey` always counts MAGs, MGnify
bins, assemblies and runs.

### Host

#### host, host_tax_id, host_keywords

- `host`: `none` for environmental biomes (a host field is then a conflict in `search`), or
  the host's scientific name (`Hordeum vulgare`).
- `host_tax_id`: the host's NCBI tax ID (`4513`).
- `host_keywords`: free-text host names (`["barley", "Hordeum"]`), because ENA hosts are
  recorded inconsistently.

**survey:** assemblies and runs must match `host_tax_id` or contain a keyword in the `host`
field; MAGs are traced through `derived from` to their source samples and kept only if
every source matches. **search:** host evidence and host conflicts.

**Where to look:** the host's tax ID by its scientific name, in
[ENA](https://www.ebi.ac.uk/ena/browser/) or [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy):

```bash
curl -s "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/Hordeum%20vulgare"
```

For `host_keywords`, the `host` values ENA samples really carry:

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=sample&query=host_tax_id=4513&fields=host&limit=200"
```

#### body_site_include

Optional, host-associated biomes only: UBERON terms for the body sites that count (e.g.
rumen). Used by `search`.

**Where to look:** [UBERON in OLS](https://www.ebi.ac.uk/ols4/ontologies/uberon).

```bash
curl -s "https://www.ebi.ac.uk/ols4/api/search?q=rumen&ontology=uberon&rows=10"
```

#### contaminant_reference

The host genome miassembler removes reads against, written to its samplesheet by `prepare`.
Empty for environmental biomes.

**Where to look:** the host references the
[miassembler](https://github.com/EBI-Metagenomics/miassembler) README accepts.

### NCBI searches

#### ncbi_isolate_query and ncbi_wgs_terms

- `ncbi_isolate_query`: an NCBI BioProject query for isolates, e.g.
  `bacteria[Organism] AND "lake sediment"`. It also finds MAGs, and some hits are reads
  only.
- `ncbi_wgs_terms`: host and isolation-source terms for the NCBI WGS search.

Both are used by `search`, which looks for genomes deposited in NCBI and not in ENA.

**Where to look:** test the query on
[NCBI BioProject advanced search](https://www.ncbi.nlm.nih.gov/bioproject/advanced) first.

### Thresholds

#### min_study_yield

For read studies that mix biomes: if accepted runs / all runs in a study is below this value
(default 0.3), the study is skipped (`low_yield_study`) instead of assembling all its runs
for a few useful ones, unless a curator set `keep`. Used by `search` and `prepare`.

## Adding or changing terms

1. **Find the term** with the resource under "Where to look" in the
   [field reference](#field-reference) and verify its ID. Never write an ID from memory.
2. **Check the boundaries.** A new lineage or ENVO term brings in its whole subtree; list its
   children and exclude the ones that belong to another biome.
3. **Record why** in `notes` when the change moves a boundary.
4. **Re-run the survey** on the config and compare the counts:
   `biome-finder survey biomes/<name>.yaml --out survey/`. A big jump after adding a generic
   metagenome taxon usually means off-biome samples.

In Claude Code, the `define-biome` skill runs these lookups and writes the config in the
template's field order.
