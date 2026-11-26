"""Constants shared by the commands.

- ENA Portal API clauses for each data layer (PROPOSAL.md Appendix B). MAGs are `wgs_set`
  records; unbinned assemblies and bins (how MGnify deposits its genomes) are `analysis`.
- Organism scope. MGnify builds one catalogue per kingdom (the pipeline's `--kingdom`), so every
  biome config sets `kingdom`; configs extend the exclusions with `exclude_taxa_extra`.
  Assemblies and reads have no organism, so the scope applies to genomes only.
"""

import re

SAMPLE_ID = re.compile(r"\b(SAM[END][A-Z]?\d+|[EDS]RS\d+)\b")  # BioSample, primary or secondary
LAYER_TYPES = {  # data type → the ENA assembly_type its layer is searched by
    "mag": "metagenome-assembled genome (mag)",
    "contigs": "primary metagenome",
    "bin": "binned metagenome",
}
MAG_LAYER = f'assembly_type="{LAYER_TYPES["mag"]}"'
CONTIG_LAYER = f'assembly_type="{LAYER_TYPES["contigs"]}"'
BIN_LAYER = f'assembly_type="{LAYER_TYPES["bin"]}"'
METAGENOMIC_WGS = 'library_source="METAGENOMIC" AND (library_strategy="WGS" OR library_strategy="OTHER")'

LONG_READ_PLATFORMS = {"OXFORD_NANOPORE", "PACBIO_SMRT"}
MIN_READS = 1_000_000  # the survey's "runs ≥ 1 M"

KINGDOM_TAXA = {"prokaryotes": [2, 2157], "eukaryotes": [2759]}  # Bacteria + Archaea; Eukaryota
EXCLUDED_TAXA = [10239, 33208, 58023, 2870]  # Viruses, Metazoa, Tracheophyta, Phaeophyceae
METAGENOMES = 408169  # the 'metagenomes' node: every '<environment> metagenome' taxon sits under it
