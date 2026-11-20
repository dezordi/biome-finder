# biome-finder

A command-line tool that finds and organizes public INSDC (ENA/NCBI/DDBJ) data for a
target biome, as input for the EBI-Metagenomics genome catalogue workflows. It does not
build catalogues itself.

For a biome defined in `biomes/<name>.yaml`, it produces:

1. **`samples.tsv`**: one row per data record (genome, metagenome assembly or read run).
   Each row has the sample, project, publication, taxonomy (as reported and resolved),
   reported biome, pre-curated biome, confidence score/tier and `needs_curation`.
2. **Workflow inputs**, routed by data type:

   | Data type | Goes to |
   |---|---|
   | `mag` / `bin` / `isolate` | [genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) (+ [datascout](https://github.com/EBI-Metagenomics/datascout) samplesheet for eukaryotes) |
   | `contigs` (unbinned assemblies) | [genomes-generation](https://github.com/EBI-Metagenomics/genomes-generation) |
   | `reads` | [miassembler](https://github.com/EBI-Metagenomics/miassembler), then genomes-generation |

The tool is biome-agnostic: environmental or host-associated, all biome knowledge lives in
the YAML configs. Organism scope is Bacteria, Archaea and microbial eukaryotes.

## Status

Early development. Only `survey` is implemented; `search`, `prepare` and
`link-assemblies` are planned (see `PROPOSAL.md`, §6 milestones).

## Install

Requires Python ≥ 3.11.

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[test]"
source .venv/bin/activate
```

## Usage

Count the public data available for candidate biomes (MAGs, assemblies and read runs per
metagenome taxon, plus what is new since an existing MGnify catalogue):

```bash
biome-finder survey --biomes biomes/*.yaml --out survey.tsv
```

Survey counts are proxies based on the metagenome taxon the submitter declared, not curated
results. API responses are cached in `~/.cache/biome-finder/` (override with
`BIOME_FINDER_CACHE`).

## Layout

- `biome_finder/`: the package
- `biomes/`: one YAML config per biome (the current files are drafts)
- `benchmark/`: manually curated tables used to calibrate scoring
- `vendor/`: fetch scripts copied unchanged from genomes-catalogue-pipeline
- `.claude/skills/`: Claude Code skills for defining biomes, checking papers and curating
  samples
- `PROPOSAL.md`: full design and curation guidelines

## Tests

```bash
python -m pytest tests/
```

## License

GPL-3.0-or-later (see `LICENSE`). The scripts in `vendor/` come from
EBI-Metagenomics/genomes-catalogue-pipeline under GPL-3.0.
