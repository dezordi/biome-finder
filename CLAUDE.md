# CLAUDE.md

## What this project is

`biome-finder` is a Python CLI that **organizes the input data for genome catalogue
workflows**. For a target biome it produces:

1. `samples.tsv`: one row per public INSDC data record (genome, metagenome assembly or read
   run) that could go into a catalogue for that biome. Each row has the sample, project,
   publication, taxonomy (as reported and resolved), reported biome, pre-curated biome,
   confidence score, and curation status. The same table also records what `prepare` did
   with each row (`destination`, `local_path`, `prepare_status`).
2. After curation, the inputs for the downstream EBI-Metagenomics workflows, routed by
   `data_type` (the most processed public data of each sample):
   - `mag` / `bin` / `isolate`: single-organism genomes at any assembly level → downloaded
     FASTAs + params for
     [genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline).
     Eukaryotes also get a samplesheet for
     [datascout](https://github.com/EBI-Metagenomics/datascout).
   - `contigs`: unbinned metagenome assemblies, e.g. MGnify ERZ → samplesheets for
     [genomes-generation](https://github.com/EBI-Metagenomics/genomes-generation)
   - `reads` → samplesheet for [miassembler](https://github.com/EBI-Metagenomics/miassembler).
     Afterwards, `link-assemblies` builds the genomes-generation inputs from its output.

Assembly, binning, QC, dereplication and annotation are done by those workflows, outside
this repo.

The tool is **biome-agnostic**: environmental or host-associated, any biome. All biome
knowledge lives in `biomes/*.yaml`. Never add biome-specific logic to the code; if a biome
needs something new, add a generic config field.

The tool will feed MGnify, so **scoping follows MGnify's conventions**: one catalogue per
kingdom (`<biome>-prok` / `<biome>-euk`) and one per host (e.g. barley, maize and tomato
rhizosphere are separate catalogues, never a merged "crop rhizosphere"). When in doubt, mirror
an existing MGnify catalogue (API v2 `/genomes/catalogues/`).

Organism scope: Bacteria, Archaea and microbial Eukaryotes. Excluded by default, using the
resolved NCBI lineage: Viruses (10239), Metazoa (33208), Tracheophyta (58023) and
Phaeophyceae (2870), plus plasmids and other MGEs, and organelle genomes. All macroalgae
are out of scope; red and green ones are judged per genus during curation. A biome config
extends the exclusions (`exclude_taxa_extra`). Like MGnify, catalogues are per kingdom:
every config sets `kingdom: prokaryotes` (Bacteria + Archaea) or `kingdom: eukaryotes`,
which is also the genomes-catalogue-pipeline `--kingdom`.

Read `PROPOSAL.md` for the design. Its Appendix A holds the curation guidelines, Appendix B
the reference ENA/NCBI queries, and Appendix C the MGnify catalogue snapshot. Don't
duplicate them here.

## Layout

- `biome_finder/`: package. `cli.py` is the click entry point only; each command's logic lives in its own module (`survey.py`, `search.py`, `prepare.py`, …). Third-party services are reached only through `clients/`, one module per service (`ena.py`, `mgnify.py`, later `ncbi.py`, `europepmc.py`, …, sharing the cached HTTP layer in `clients/http.py`), with their query and response models in `models.py`. Constants shared by commands (ENA data-layer clauses, organism scope) live in `constants.py`. Scoring lives in `score.py` and the LLM in `llm.py`
- `biomes/*.yaml`: one config per biome (lineages, GOLD/ENVO terms, host/body site,
  include/exclude keywords, data types, notes). The file is the biome decision. Read it with
  `yaml.safe_load`.
- `benchmark/`: manually curated tables used to calibrate scoring. `marine_euks.tsv` is the
  first one. Its `Curated Biome` / `Keep` columns are
  the biome ground truth, and `Category` is the data-type ground truth. The
  `curate-sample` skill appends to `benchmark/curated_<biome>.tsv`: fixed columns, one row
  per record, and `actions` lists the curation actions performed.
- `.claude/skills/`: `define-biome`, `check-paper-data`, `curate-sample`. They are
  self-contained, each with its own copy of the rules it needs. When a rule in
  PROPOSAL.md changes, also update the skills that copy it (mapping in PROPOSAL.md §8).
- `docs/`: Sphinx docs (Read the Docs theme), one page per command. When a command's
  behaviour or output columns change, update its page. Check with
  `sphinx-build -W -b html docs docs/_build/html`.
- `vendor/`: `fetch_ena.py`, `fetch_ncbi.py` and helpers copied **unchanged** from
  EBI-Metagenomics/genomes-catalogue-pipeline (pinned commit in `vendor/README.md`, GPL-3.0).
  Run them as subprocesses; don't import, edit or reformat them. Also holds the IndicPiper
  `genus_habitat_indicators_v<version>.csv` database, pinned and unchanged.

## Commands

```bash
# system python is 3.9 (no PyYAML); use the project venv (Python 3.12 via uv)
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e ".[test]"
# then run tools as .venv/bin/biome-finder / .venv/bin/python -m pytest
biome-finder survey  --biomes biomes/*.yaml --out survey/   # survey.tsv, biomes/*/, report.html
biome-finder search  --biome biomes/<name>.yaml --survey survey/ --out samples.tsv \
                     [--gold goldData.xlsx] [--llm] [--composition [all]] [--since YYYY-MM-DD] [--taxdump YYYY-MM-DD] \
                     [--papers-dir papers/]
biome-finder prepare --samples samples.tsv --biome biomes/<name>.yaml --tiers A,B --out <dir>
biome-finder link-assemblies --samples samples.tsv --miassembler-outdir <path> --out <dir>
python -m pytest tests/
```

## Curation rules the code must enforce

These come from PROPOSAL.md Appendix A. Never relax them silently.

- Data must be **public in INSDC (ENA/NCBI/DDBJ)**. MAGs hosted only on Zenodo, Dryad, CNGB,
  NGDC or a lab portal don't count as genomes. Their samples appear as `contigs`/`reads` if
  that data is public in INSDC.
- Rows are data records, but the **sample** (BioSample) is the unit of biome curation: score
  and tier are per sample, `decision` is per row. A co-assembled genome is accepted only if
  all its source samples are accepted.
- Exclude: controls, mock and synthetic metagenomes, spiked samples, organelle genomes,
  out-of-scope taxa, amplicon-only data, and off-biome material (host conflict, barn faeces
  for a gut biome, and so on; a sample whose own evidence points only to another biome is
  `off_biome`, whatever its study says). An enrichment culture is not excluded automatically: it is
  fine if no organisms were added while culturing, so it goes to curation.
- Mixed-biome read studies: assemble **all** WGS runs of the study, controls and off-biome
  runs included (`assembly_only`), but bin only the accepted runs. Skip the study
  (`low_yield_study`) if accepted runs / all runs < `min_study_yield`, unless a curator set
  `keep`.
- Genomes we didn't generate must be traceable to the original biological sample (via
  `derived from`). Map NCBI genomes to GCA **by BioSample, not by project**.
- Never mix MAGs and bins from the same study. Default to MAGs; the YAML enables bins.
- **Data type needs its own evidence** (ENA `assembly_type` and assembly `DESCRIPTION`,
  NCBI `excluded_from_refseq`, BioSample package, `metagenome-source`/`derived_from`, paper).
  Don't trust the word "isolate". `isolate` requires pure-culture evidence at **sample**
  level or in the paper, or tissue / whole individuals of a wild-collected organism.
- ENA `env_tax:*` / `env_geo:*` sample tags are derived from the organism's taxonomy, so
  they are never biome evidence.
- Biome scope is **subtrees minus excluded subtrees**: `mgnify_lineages` minus
  `mgnify_lineages_exclude`, `envo_include` minus `envo_exclude` (ENVO terms include their
  subclasses). An excluded term in any `env_*` field rules the sample out.
- Marine water at the coast counts as ocean (`…:Marine:Oceanic`). Encode it per biome with
  `lineage_equivalents`, never in code.
  SAGs are not isolates and are excluded by default. Exclude type strains with no
  confirmed source. If the type stays unclear, use `unknown` → `needs_curation`.
- `contigs` means an unbinned metagenome assembly. A single genome at contig or scaffold
  level is still `mag`/`isolate`.
- Reads: list every WGS metagenomic run of an accepted sample. Route only short reads,
  paired or merged, with at least one run ≥ 200k reads. Keep long reads in the table
  (`long_read`, `deferred_long_read`), but don't route them for now.
- Each sample enters at its most processed public data: genomes > contigs > reads.
  Lower-level rows stay in the table as `not_selected`. Never re-assemble or re-bin
  something that is already public. Check this by sample (every assembly or WGS set ENA links
  to it, of any type; PROPOSAL.md §3.10, B.6), not only through the searches by taxon. Only a
  genome of the catalogue's kingdom covers a sample.
- Deduplicate ENA sets against NCBI GCAs by BioSample. Flag SPIRE bins and re-assemblies
  (mag-provenance). For updates, flag genomes already in the existing catalogue.
- **When in doubt, leave it out.** A wrong-biome genome in a catalogue is worse than a
  missing one.

## Scoring conventions

- Evidence items are `(action, source, claimed_lineage, specificity, weight, quote)`. The
  actions mirror what curators did in the benchmark's `Manual Inspection` column: sample
  metadata, project description, sibling samples, paper methods, supplementary tables,
  culture collection page.
- Run the cheap actions first and stop once the evidence is conclusive. The LLM runs only on
  records still uncertain after the rule-based steps.
- An LLM evidence item must include a verbatim quote that is present in the retrieved text,
  or it is dropped. LLM-only evidence can reach tier B at most, never A.
- Any conflict with weight ≥ 0.6 sets `needs_curation`, whatever the score.
- Composition evidence (action 10: Trapiche + IndicPiper) is independent of the metadata.
  By default it runs only for samples still uncertain after actions 0–9; `--composition all`
  also uses it as a mislabel check. Like LLM evidence, it can never reach tier A on its own.
  A disagreement with strong metadata sets `composition_conflict` → `needs_curation`.
- IndicPiper: reimplement the genus exact-match + per-habitat sum in Python on Sandpiper
  GTDB profiles. Don't add R. Trapiche: call its CLI as a subprocess (NDJSON in/out) in its
  own environment, passing our text evidence as `ext_text_pred_*`. Never import it, because
  it pulls in TensorFlow.
- A filled `decision` column (human curation) always overrides `tier`.
- After changing any weight, threshold or rule, re-run the benchmark and report the tier A
  precision, the tier C share and the data-type agreement.

## Papers

- Request a paper only when it can still change a decision; records already in tier A or
  excluded are `not_needed`.
- Non-open-access papers get `paper_access = needs_pdf`. The curator drops the PDFs into
  `--papers-dir` with any name. `search` finds each PDF's DOI (XMP/Info metadata → regex
  on the first two pages → title match) and renames the file to the canonical name: the
  DOI, lowercased, with `/` replaced by `_`, plus `.pdf`
  (`10.1038/s41564-023-01435-6` → `10.1038_s41564-023-01435-6.pdf`).
- A match counts only if the DOI is one of the table's `needs_pdf` DOIs, because the first
  pages cite other papers. With zero or several matches, leave the file untouched and
  report it; never guess.
- Supplements can't be identified automatically; the curator names them
  `<stem>.supp<N>.<ext>`.
- `--papers-dir` picks up the provided files (`pdf_provided`). Text is extracted with
  `pypdf`, so LLM quotes can be checked against it.

## Taxonomy

- Keep `reported_tax_id` / `reported_tax_name` exactly as shown in the record.
- Resolve `tax_lineage` and `domain` with the NCBI taxdump archive release chosen by
  `--taxdump YYYY-MM-DD` (`new_taxdump_<YYYY-MM>-01.zip`; the default is the latest). Follow
  `merged.dmp`, and flag `delnodes.dmp` hits as `taxid_deleted`. ENA keeps the reported tax ID
  and name current, so no resolved tax ID or name is stored.
- `domain` and the out-of-scope filter always use the resolved lineage; the datascout `taxid`
  is `reported_tax_id`.

## Re-running and editing `samples.tsv`

- `search` on an existing `samples.tsv` merges by `record_accession`. It refreshes the
  automatic columns and never touches `decision` or human-written `curation_note`.
- `prepare` and `link-assemblies` only write `destination`, `local_path` and
  `prepare_status`, and they save `samples.tsv.bak` before rewriting the table.

## External services

- ENA Portal/Browser/xref, NCBI E-utilities, FTP and taxdump archive, MGnify API, Europe PMC,
  GOLD exports (downloaded by the curator: the page is behind a browser check, never bypass it), Sandpiper API (GTDB run profiles), MGnify analysis taxonomy files. Cache responses on disk (default `~/.cache/biome-finder/`, keyed by URL or
  accession).
- Respect rate limits: NCBI allows 3 requests/s without an API key and 10 with one (via the
  `NCBI_API_KEY` env var). Use backoff on HTTP 429 and 5xx.
- `prepare` downloads genomes only. For contigs and reads it writes ENA FTP URLs and md5s,
  and never downloads those files.
- Don't put personal emails or names in User-Agent headers or API params. Secrets
  (`ANTHROPIC_API_KEY`, `NCBI_API_KEY`) come only from env vars and are never written to files.

## Downstream input contracts

These were read from each pipeline's `assets/schema_input.json` and README on 2026-11-19.
Validate outputs against those schemas in tests; if upstream changes, fix the writer, not
the test.

- **miassembler** `--samplesheet` (CSV): `study_accession` (secondary: ERP/SRP/DRP),
  `reads_accession`, `fastq_1`, `fastq_2`, `library_layout` (`single`/`paired`),
  `library_strategy` (`metagenomic`…), and optionally `contaminant_reference`. FASTQ paths
  must end in `.fastq.gz`/`.fq.gz`.
- **genomes-generation** `--input` (CSV, one study per run, one assembly per run, no
  co-binning, short reads only): `id` (run), `assembly_fasta` (`.fa.gz`/`.fasta.gz`),
  `fastq_1`, `fastq_2`, `assembly_accession` (ERZ). Also `--assembly_software_file`
  (run → `Tool_vVersion`), `--metagenome` (ENA metagenome taxon name) and `--biomes`
  (`"environment_biome,environment_feature,environment_material"`).
- **genomes-catalogue-pipeline**: uncompressed genomes, with ENA and NCBI in separate
  folders (`--ena_genomes`, `--ncbi_genomes`); `--ena_genomes_checkm` is the CSV written by
  `fetch_ena.py`; `--biome` is an MGnify lineage; one run per `--kingdom`
  (prokaryotes/eukaryotes); eukaryotes also take `--protein_evidence` from datascout. Params
  go in a `-params-file` JSON. MGnify-only values (catalogue name, version, MGYG range) are
  placeholders.
- **datascout** `--samplesheet`: `sample_id`, `genome` (GCA), `taxid` (resolved), optional
  `assembly_file`.

## Code style

- Keep it small: stdlib first (`csv`, `json`), plus `click`, `requests`, `PyYAML` and `openpyxl`
  (only to read the GOLD export). The
  optional extras are the LLM SDK and `pypdf`. Ask before adding any other dependency.
- Reports: `report.py` only, with pandas, Panel and plotly. The survey report is one
  self-contained dashboard (`report.html`): Panel components (tabs, cards, `Number`
  indicators, `Tabulator` tables) holding one plotly express figure per card, saved with
  `save(embed=True, max_opts=<catalogues>)` so every catalogue view is precomputed and the
  file works without a server. Panel templates can't be embedded, so don't use them. Use
  the libraries' own features; never write HTML by hand. Command logic stays stdlib; no
  pandas outside reports.
- Third-party calls are typed. Every service gets a section in `models.py`: a query model
  where the query has structure (e.g. `EnaSearch`; a plain identifier otherwise) and frozen
  dataclasses for what comes back, with only the fields biome-finder uses (for ENA they are
  also the fields requested). Clients return `Response[T]`: the HTTP status (404 = not found;
  other errors raise) plus the data. Parse rows with `from_row`, and keep a source's value
  as `str` when its format varies (ENA `lat` can be `;`-joined). Command modules never call
  `clients/http.py` or build URLs themselves.
- Type hints on every function. Check with
  `uvx mypy --ignore-missing-imports --python-executable .venv/bin/python biome_finder`.
- Format with black (`uvx black biome_finder tests docs`, line length 120 in `pyproject.toml`).
  Keep comments few: explain *why* in the function's docstring, not in loose comments.
- `samples.tsv` has a fixed column order (PROPOSAL.md §3.6). Columns are added, never
  renamed, because curators keep spreadsheets of past runs.
- Every non-trivial rule in `score.py` has a test, ideally using a real accession from the
  benchmark.
- Code, comments, docs and commit messages are in English.

## Git

- Commit or push only when asked. Run outputs (`survey/`, `samples.tsv`, `papers/`) are
  git-ignored.
