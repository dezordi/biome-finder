---
name: check-paper-data
description: Given a DOI, PMID, paper or MAG-collection name (e.g. an awesome-mag entry), find its data-availability statement, list every accession, check each in ENA/NCBI and classify how biome-finder can use it — genomes in INSDC, only contigs/reads in INSDC (regenerate MAGs), or not in INSDC (unusable). Use when triaging seed collections, a paper the user mentions, or "is this dataset usable / where are the MAGs deposited".
---

# check-paper-data

## Rules

- **Usable data must be public in INSDC: ENA, NCBI or DDBJ only.** CNGB/CNSA, NGDC
  (GSA/GWH), NMDC-China, Zenodo, Figshare, Dryad, JGI IMG and lab portals are not INSDC.
- MAG collection papers often keep their genomes outside INSDC. If the reads or
  assemblies are public in INSDC but the MAGs never were, the dataset is still usable:
  the MAGs get regenerated (`contigs` rows go to genomes-generation; `reads` rows go to
  miassembler, then genomes-generation).
- NCBI-deposited genomes don't show up in ENA's MAG layer. They appear as GCA accessions
  elsewhere in ENA, so check NCBI as well.
- **SPIRE** genomes sit in ENA's MAG layer, but they were never dereplicated: by ENA's
  definition they are **bins**. SPIRE also re-assembles datasets that may already be in
  a catalogue, so flag them for
  [mag-provenance](https://github.com/EBI-Metagenomics/mag-provenance/tree/dev).
- Always assume ENA/NCBI metadata may contain errors.

## Known routes (from awesome-mag's access notes; confirm per study)

| Route | Resources |
|---|---|
| `genomes_in_insdc` | MGnify Genomes (ENA) · OceanDNA (DDBJ PRJDB11811) · Glacier-fed Streams MAGs (NCBI) · GROWdb (NCBI PRJNA946291) · RGMGC (NCBI Assembly) · SPIRE (ENA MAG layer, actually bins) · Human Gut Archaeome (check) |
| `contigs_or_reads_in_insdc` | LakePulse (ENA reads PRJEB29238) · RUG2 (ENA reads + assemblies) · Meta320 (SRA) · GEM (IMG) · SMAG (Zenodo) · Bin Chicken RBG (Zenodo) · cFMD (Zenodo) · Microbiome Datahub (NIG) · mOTUs DB · OMDB · HRGM / HROM / MRGM |
| `not_insdc` | GOMC · TPMC / TPMC-S · QXLSG · PIGC · Caprinae gut · UHSG · ICRGGC |
| Aggregators (only for finding source studies) | MAGdb (offline) · gcMeta · Anammox catalog |
| Out of scope (viral) | GPD · MGV · UHGV · MetaVR |

## Steps

1. **Resolve the paper** in Europe PMC:
   ```bash
   curl -s "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:<doi>&format=json&resultType=core"
   ```
   Note `pmcid`, `isOpenAccess` and the title.

   **Starting from an accession instead of a DOI**, look for the paper in this order:
   1. ENA xref:
      ```bash
      curl -s "https://www.ebi.ac.uk/ena/xref/rest/json/search?accession=<PRJ…>&source=EuropePMC"
      ```
      ENA lists both "Publications in which the generation of this record is described"
      and "Publications citing this record". The second kind includes papers that
      merely cite the study, so check that the paper actually generated the data.
   2. A web search for the PRJ accession and for the ERP/SRP/DRP accession.
   3. If there's no paper, the ENA study description is often enough.

2. **Read the Data-availability and Methods sections.**
   - Open access:
     ```bash
     curl -s "https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML"
     ```
   - Not open access: look in the user's papers folder for the file named after the DOI,
     lowercased, with `/` replaced by `_`
     (`10.1038/s41564-023-01435-6` → `10.1038_s41564-023-01435-6.pdf`). If it isn't
     there, ask for it. Don't guess accessions from the abstract.

3. **Extract the accessions** from the text and the supplements:
   - Studies and projects: `PRJ(EB|NA|DB)\d+`, `[EDS]RP\d+`
   - Samples: `SAM(EA|N|D)\d+`, `[EDS]RS\d+`
   - Runs: `[EDS]RR\d+`
   - Assemblies and genomes: `ERZ\d+`, `GCA_\d+\.\d+`, WGS set prefixes
   - Non-INSDC IDs, which should be reported: CNP/CNS, CRA/PRJCA, GWH, and Zenodo or
     Figshare DOIs

4. **Check each INSDC study in ENA.** For each, find what exists, how much of it is
   public, and the `assembly_type`:
   ```bash
   B=https://www.ebi.ac.uk/ena/portal/api/search
   curl -s "$B?result=wgs_set&query=study_accession=<PRJ>&fields=accession,assembly_type,sample_accession&format=tsv&limit=0"
   curl -s "$B?result=analysis&query=study_accession=<PRJ>&fields=analysis_accession,assembly_type,sample_accession&format=tsv&limit=0"
   curl -s "$B?result=read_run&query=study_accession=<PRJ>&fields=run_accession,sample_accession,library_strategy,library_layout,instrument_platform,read_count&format=tsv&limit=0"
   ```
   - **Reading `assembly_type`:**
     - `Metagenome-Assembled Genome (MAG)` → MAGs
     - `binned metagenome` → bins
     - `primary metagenome` → contigs
     - `clone or isolate` → isolate candidates
   - **MAGs are often in a different study than the reads** (a TPA/assembly study).
     Follow "Related ENA Records" both for the study and for one of its samples (ENA
     xref).
   - **NCBI-only genomes:** look up the GCA or BioSample in
     `https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt` or with
     E-utilities. Map **by BioSample, not by project**: projects can span several biomes.
   - **Usability notes for reads:**
     - single-end short reads yield few MAGs, unless they are merged pairs or
       high-coverage
     - long reads are tracked but not routed for now
     - a usable study needs at least one run with ≥ 200k reads (millions preferred)

5. **Classify** the dataset with one route. Several routes can apply to a mixed case.
   - `genomes_in_insdc`: MAG, bin or isolate genomes public in INSDC. Give counts by
     `assembly_type`, and flag SPIRE.
   - `contigs_or_reads_in_insdc`: genomes are elsewhere, but ERZ assemblies or WGS runs
     are public. Give counts, the paired / single / long split, and the read-count range.
   - `not_insdc`: nothing usable in INSDC. Say where the data actually is.
   - `unclear`: the accessions don't resolve, or the data isn't public yet.

## Output to the user

- The paper: title, DOI, open access or not.
- The route.
- A table: accession → type → public count → note.
- The biome as the paper states it, as a verbatim quote with its section.
- Caveats:
  - controls, mock or synthetic samples
  - samples from other biomes (mixed studies)
  - spiked or enriched samples
  - whether it's really WGS
  - co-assembly
  - long reads

Don't write any file unless the user asks.
