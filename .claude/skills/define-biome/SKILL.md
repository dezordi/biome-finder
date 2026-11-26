---
name: define-biome
description: Write or revise a biome config (biomes/<name>.yaml) for biome-finder — MGnify lineage, GOLD path, verified ENVO terms, metagenome taxa, IndicPiper habitat, include/exclude keywords, host rules and boundary decisions. Use when the user wants to define, create, split, merge or adjust a target biome or catalogue scope, or before running `biome-finder survey` on a new biome.
---

# define-biome

The YAML file is both biome-finder's input and the written record of the biome decision.
The tool is biome-agnostic: everything biome-specific goes in this file, never in code.
Don't invent fields. If a biome needs something new, propose a generic field to the user.

## Scoping rules

- **Environmental biomes:** one general catalogue or several subdivided ones? For example,
  a single marine catalogue, or separate marine water and marine sediment ones; one soil
  catalogue, or one per soil type. `biome-finder survey` gives the counts; the decision is
  recorded in the YAML.
- **Host-associated biomes:** decide what the biome covers (cow rumen only, or the whole cow
  gut?). Environmental material near the host (e.g. barn faeces for a farm-animal gut
  catalogue) is excluded unless there is a good reason to include it.
- **Follow MGnify's pattern; don't invent groupings.** biome-finder will feed MGnify, so
  host-associated catalogues are **one per host species** (MGnify has separate barley,
  maize and tomato rhizosphere catalogues; never merge them into "crop rhizosphere"), and
  every biome gets one config per kingdom. Check the existing catalogues for precedent.
- **MAGs vs bins:**
  - ENA definitions: a **bin** is any genome assembled from metagenomic data. A **MAG**
    is a bin dereplicated across a biome within a project, so a project has at most one
    (best-quality) MAG per species per biome.
  - **Default to MAGs.** That keeps catalogues small, but it limits pan-genomes: one MAG
    study gives none, and two give at most two genomes per species.
  - **Enable bins** only when one of these holds:
    - a project requires all of its data to be in the catalogue
    - the catalogue comes from a single ENA project and needs pan-genomes
    - very little data exists for the biome and no large influx is expected
  - Bins and MAGs from the same study are never mixed (the MAGs are a subset of the bins).
- **Organism scope:** Bacteria, Archaea and microbial Eukaryotes.
  - **Excluded by default:** Viruses (10239), Metazoa (33208), Tracheophyta (58023),
    Phaeophyceae (2870), plasmids/MGEs and organelle genomes.
  - **Macroalgae are out of scope**, including red and green seaweeds, which are judged
    per genus during curation.
  - `exclude_taxa_extra` adds to the exclusions.
- **One kingdom per config.** MGnify builds separate prokaryote (Bacteria + Archaea) and
  eukaryote catalogues, and the pipeline runs once per `--kingdom`. Every config sets
  `kingdom: prokaryotes` or `kingdom: eukaryotes`; a biome needing both gets two configs
  named `<biome>-prok` and `<biome>-euk` (e.g. `soil-prok` and `soil-euk`).
- **Marine water at the coast counts as ocean.** For marine water biomes, map
  `…:Marine:Coastal` to `…:Marine:Oceanic` in `lineage_equivalents`. Coastal sediment,
  intertidal, estuarine and brackish environments stay separate.

## YAML template (field order is fixed)

```yaml
name: freshwater-lake-sediment
kingdom: prokaryotes                  # prokaryotes | eukaryotes: MGnify builds one catalogue per kingdom (pipeline --kingdom)
mgnify_catalogue: ""                 # set for updates, e.g. soil-v1-0
catalogue_biome: "root:Environmental:Aquatic:Freshwater:Sediment"   # → genomes-catalogue-pipeline --biome (must exist in MGnify)
mgnify_lineages:                     # whole subtrees count (each lineage plus all its descendants)
  - "root:Environmental:Aquatic:Freshwater:Sediment"
mgnify_lineages_exclude: []         # subtrees removed from the scope, e.g. a :Sediment child of a water lineage
envo_include:                        # term + all its ENVO subclasses; every ID verified in OLS
  - ENVO:00000546                    # lake sediment
envo_exclude:                        # subtrees removed (term + subclasses), checked on every env_* field
  - ENVO:03000033                    # marine sediment
gold_ecosystems:                     # GOLD 5-level paths
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
data_types: [mag, isolate, contigs, reads]   # add bin only in the cases above; sag only if explicitly wanted
exclude_taxa_extra: []               # added to the default exclusions (Viruses, Metazoa, Tracheophyta, Phaeophyceae)
lineage_equivalents: {}              # evidence lineage → lineage it counts as, e.g.
                                     # {"root:Environmental:Aquatic:Marine:Coastal": "root:Environmental:Aquatic:Marine:Oceanic"}
ncbi_isolate_query: ""               # NCBI BioProject query for isolates, e.g. 'bacteria[Organism] AND "lake sediment"'
ncbi_wgs_terms: []                   # host / isolation_source terms for the NCBI WGS search
min_study_yield: 0.3                 # mixed-biome read studies: skip if accepted runs / all runs is below this
indicpiper_habitats: ["freshwater sediment"]  # IndicPiper habitat(s) that count as this biome
notes: "Surface and core sediments. Saline/soda lakes excluded."
```

`data_types` values:

- `mag`, `bin`, `isolate`: single-organism genomes at any assembly level
- `sag`: single amplified genome (not an isolate)
- `contigs`: unbinned metagenome assemblies, e.g. MGnify ERZ
- `reads`: raw metagenomic runs

## Steps

1. **Pin down the scope with the user** before any lookup, applying the scoping rules
   above. Ask only what applies:
   - **Environmental:**
     - split or keep whole (water vs sediment, lake vs river, soil types)?
     - boundaries: freshwater vs brackish/saline, surface vs core, rhizosphere vs bulk
       soil
   - **Host-associated:**
     - host taxon?
     - which body sites count?
     - is environmental material near the host excluded (default: yes)?
   - **New or update?** For an update, get the `mgnify_catalogue` id:
     ```bash
     curl -s "https://www.ebi.ac.uk/metagenomics/api/v2/genomes/catalogues/"   # v2: v1 lists outdated versions
     ```
     Its `catalogue_type` (prokaryotes/eukaryotes) must match the config's `kingdom`; the
     survey refuses a mismatch.
   - **Bins:** only if one of the three bin conditions above holds.

2. **MGnify lineages.** Check that every lineage exists and list its children, so the user
   can pick exact nodes:
   ```bash
   curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/<lineage>"                     # exists? samples-count
   curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/<lineage>/children?page_size=100"
   ```
   `catalogue_biome` is exactly one lineage; `mgnify_lineages` can list several.
   - **Scope = subtrees minus excludes.** Each `mgnify_lineages` entry brings in its whole
     subtree. Put any descendant that belongs to another biome in `mgnify_lineages_exclude`
     (its whole subtree goes), so water and sediment, or soil and permafrost, never mix.
     Example: a soil biome lists `root:Environmental:Terrestrial:Soil` and excludes
     `…:Soil:Permafrost` and `…:Soil:Wetlands`.
   - Read the whole children listing: MGnify often files the same environment under
     several parents (freshwater sediment sits under `Freshwater:Sediment`,
     `Lentic:Sediment`, `Lotic:Sediment` and `Pond:Sediment`; wetlands under both
     `Terrestrial:Soil:Wetlands` and `Aquatic:Freshwater:Wetlands`).

3. **ENVO terms.** Verify every ID in OLS. Never write an ID from memory.
   ```bash
   curl -s "https://www.ebi.ac.uk/ols4/api/search?q=<term>&ontology=envo&rows=10"   # obo_id + label
   ```
   - `envo_include` and `envo_exclude` both mean the term **plus all its ENVO subclasses**.
     A sample is excluded if any `env_*` field (broad scale, local scale, medium) falls in
     an excluded subtree, even when another field matches `envo_include` (e.g. broad scale
     "marine biome" + medium "sediment" is out of a marine-water biome). Check the
     hierarchy before choosing: `curl -s "https://www.ebi.ac.uk/ols4/api/ontologies/envo/terms/<double-encoded IRI>/hierarchicalAncestors"`.
   - Use `envo_exclude` for the same splits as `mgnify_lineages_exclude` (water vs
     sediment, soil vs permafrost/wetland/rhizosphere).

   For body sites of host-associated biomes, use UBERON (`ontology=uberon`) and add a
   `body_site_include` list.

4. **Metagenome taxa.** These are ENA/NCBI organism names under "metagenomes", taxid
   408169.
   ```bash
   curl -s "https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-submission/<words>%20metagenome"
   ```
   Keep only entries with `"metagenome": "true"` that fit the scope. Generic names
   ("metagenome", "sediment metagenome") are weak evidence and pull in off-biome samples;
   include them only if the user agrees.

5. **GOLD ecosystem path** (5 levels, e.g.
   `Environmental > Aquatic > Freshwater > Lake > Sediment`). If you can't confirm it,
   write `# TODO verify` rather than guess.

6. **IndicPiper habitat(s).** Choose only from these 13 (database v2.1.0; soil includes
   rhizosphere): soil, marine sediment, marine water, mammalian gut, freshwater sediment,
   freshwater water, human oral, wastewater, glacier, hydrothermal vent, phyllosphere,
   human vaginal, human skin. If none fits, use `[]`.

7. **Keywords, host, data types and thresholds.**
   - `keywords_include`: specific phrases.
   - `keywords_exclude`: boundary neighbours (e.g. "marine", "estuar" for a freshwater
     biome) plus manipulations (enrichment, mesocosm, spiked, mock, synthetic, control).
   - `host: none` for environmental biomes.
   - `ncbi_isolate_query`: an NCBI BioProject query. It finds isolates *and* MAGs, and some
     hits are reads only.
   - Keep the default `min_study_yield` unless the user says otherwise.

8. **Write `biomes/<name>.yaml`** in the template's field order. Put the step-1 boundary
   decisions in `notes`, in plain words; that's the record of why.

9. **Check:**
   - Parse the file: `python -c "import yaml,sys; yaml.safe_load(open(sys.argv[1]))" biomes/<name>.yaml`.
   - If the CLI exists, run `biome-finder survey --biomes biomes/<name>.yaml` and show the
     counts. Low counts for a data type, or mostly generic evidence, are reasons to revisit
     step 1.

## Output to the user

- The YAML path.
- A table: field → value → source URL, or "user decision".
- Anything left as `TODO verify`.
