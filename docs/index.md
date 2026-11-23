# biome-finder

biome-finder is a command-line tool that finds and organizes public INSDC (ENA/NCBI/DDBJ)
data for a target biome, as input for the EBI-Metagenomics genome catalogue workflows.

For a biome defined in `biomes/<name>.yaml`, it produces:

1. `samples.tsv`: one row per data record (genome, metagenome assembly or read run). Each
   row has the sample, project, publication, taxonomy (as reported and resolved), reported
   biome, pre-curated biome, confidence score/tier and `needs_curation`.
2. Workflow inputs, routed by data type:

   | Data type | Goes to |
   |---|---|
   | `mag` / `bin` / `isolate` | [genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) (plus a [datascout](https://github.com/EBI-Metagenomics/datascout) samplesheet for eukaryotes) |
   | `contigs` | [genomes-generation](https://github.com/EBI-Metagenomics/genomes-generation) |
   | `reads` | [miassembler](https://github.com/EBI-Metagenomics/miassembler), then genomes-generation |

The tool is biome-agnostic: environmental or host-associated, all biome knowledge lives in
the YAML configs. Organism scope is Bacteria, Archaea and microbial eukaryotes.

## How the tool works

biome-finder runs in steps. Each step reads more sources than the one before it, and a
curator decides between them.

```{image} img/pipeline-overview.svg
:alt: biome-finder steps (survey, search, prepare, link-assemblies) and the data each one reads
:width: 100%
```

1. **[survey](survey.md)** answers "is there enough public data for this catalogue?". It
   reads only metadata already public in ENA, plus MGnify for context (existing
   catalogues and their genomes). It reads **no papers**, runs no LLM and queries no NCBI,
   GOLD or composition data. Its counts are proxies, based on the metagenome taxon each
   submitter chose.
2. **search** builds `samples.tsv` for the biomes the curator picked, scoring each record.
   Rule-based metadata checks run first. Only records still uncertain go on to papers
   (Europe PMC open access, plus PDFs the curator adds), the optional LLM layer and the
   composition evidence (Trapiche on MGnify taxonomy profiles, IndicPiper on Sandpiper
   GTDB profiles). The curator decides the rows left in `needs_curation` and
   re-runs `search`.
3. **prepare** routes each curated sample by its most processed public data: genomes to
   genomes-catalogue-pipeline, assemblies to genomes-generation, reads to miassembler.
4. **link-assemblies** builds the genomes-generation inputs from the miassembler output.

To regenerate the figure: `python docs/img/pipeline_overview.py docs/img/pipeline-overview.svg`.

## Status

Early development. Only [survey](survey.md) is implemented; `search`, `prepare` and
`link-assemblies` are planned.

## Installation

Requires Python 3.11 or later.

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[test]"
source .venv/bin/activate
biome-finder --help
```

To build this documentation locally:

```bash
uv pip install --python .venv/bin/python -e ".[docs]"
sphinx-build -b html docs docs/_build/html
```

```{toctree}
:maxdepth: 2
:caption: Guide

biome-config
survey
```
