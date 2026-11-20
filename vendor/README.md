# vendor/

Copied unchanged from
[EBI-Metagenomics/genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline/tree/master/bin)
at commit `a86c346cca2ade3886e49aca01af19a1695ae9bc` (copied 2026-11-19).

| File | Role |
|---|---|
| `fetch_ena.py` | Downloads MAGs or bins from ENA, by project accession |
| `fetch_ncbi.py` | Downloads genomes from NCBI, by GCA accession |
| `utils.py`, `get_ENA_metadata.py` | Helpers imported by the two scripts above |

The scripts are licensed under GPL-3.0 (see the file headers). `biome-finder` runs them as
subprocesses and does not import them. Don't edit these files. To pick up an upstream fix,
copy the files again and update the commit hash above.

Extra runtime dependencies: `requests`, `xmltodict`, `retry`.
