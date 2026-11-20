# biome-finder

biome-finder is a command-line tool that finds and organizes public INSDC (ENA/NCBI/DDBJ)
data for a target biome, as input for the EBI-Metagenomics genome catalogue workflows. It
does not build catalogues itself: assembly, binning, QC, dereplication and annotation are
done by those workflows.

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
:caption: Commands

survey
```
