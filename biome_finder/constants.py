"""Every constant of biome-finder, by the module that uses it.

- ENA Portal API clauses for each data layer (PROPOSAL.md Appendix B). MAGs are `wgs_set`
  records; unbinned assemblies and bins (how MGnify deposits its genomes) are `analysis`.
- Organism scope. MGnify builds one catalogue per kingdom (the pipeline's `--kingdom`), so every
  biome config sets `kingdom`; configs extend the exclusions with `exclude_taxa_extra`.
  Assemblies and reads have no organism, so the scope applies to genomes only.
- Then each module's own: table columns, scoring weights and patterns, service endpoints, and
  the dashboard's labels.
"""

import os
import re
from pathlib import Path

SAMPLE_ID = re.compile(r"\b(SAM[END][A-Z]?\d+|[EDS]RS\d+)\b")  # BioSample, primary or secondary
LAYER_TYPES = {  # data type → the ENA assembly_type of its layer
    "mag": "metagenome-assembled genome (mag)",
    "contigs": "primary metagenome",
    "bin": "binned metagenome",
    "sag": "environmental single-cell amplified genome (sag)",  # left out of search's genome query
}
GENOME_FIELDS = [  # sample fields of an ENA assembly where a submitter states the environment
    "isolation_source",
    "environment_biome",
    "environment_feature",
    "environment_material",
    "broad_scale_environmental_context",
    "local_environmental_context",
    "environmental_medium",
]
MAG_LAYER = f'assembly_type="{LAYER_TYPES["mag"]}"'
CONTIG_LAYER = f'assembly_type="{LAYER_TYPES["contigs"]}"'
BIN_LAYER = f'assembly_type="{LAYER_TYPES["bin"]}"'
METAGENOMIC_WGS = 'library_source="METAGENOMIC" AND (library_strategy="WGS" OR library_strategy="OTHER")'

LONG_READ_PLATFORMS = {"OXFORD_NANOPORE", "PACBIO_SMRT"}
MIN_READS = 1_000_000  # the survey's "runs ≥ 1 M"

KINGDOM_TAXA = {"prokaryotes": [2, 2157], "eukaryotes": [2759]}  # Bacteria + Archaea; Eukaryota
EXCLUDED_TAXA = [10239, 33208, 58023, 2870]  # Viruses, Metazoa, Tracheophyta, Phaeophyceae
METAGENOMES = 408169  # the 'metagenomes' node: every '<environment> metagenome' taxon sits under it


# --- NCBI ranks (models.py) ---

NCBI_RANKS = ("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")


# --- Scoring (score.py): curation actions, flags, weights ---

# Action 0, hard filters: a sample whose own text names one of these is not a natural sample of
# any biome. Patterns were checked against freshwater-water-prok (2026-11): bare "control" there
# is an experiment's untreated arm (a real sample, left to the biome evidence), and "mock" is
# mostly ENA's "mock sample for coassembly", the virtual sample of a co-assembly. An enrichment
# culture is fine if no organisms were added while culturing, which only a curator can tell, so
# `enrichment` sends the row to curation (CURATION_FLAGS) rather than excluding it.
HARD_FILTERS = {
    "control": re.compile(
        r"\b(?:negative|blank|extraction|kit|reagent|no[- ]template)[ ._-]?controls?\b|\bblank|\bntc\b", re.I
    ),
    "mock": re.compile(r"\bmock\b(?! sample for co-?assembly)", re.I),  # see CO_ASSEMBLY_SAMPLE
    "synthetic": re.compile(r"\bsynthetic\b", re.I),
    "enrichment": re.compile(r"\benrichments?\b|\benriched (?:microbiome|culture|communit)", re.I),
    "spiked": re.compile(r"\bspiked?\b|\bspike-?ins?\b", re.I),
    "virome": re.compile(  # the viral fraction (virions, VLPs): no cellular genomes to recover
        r"\bviromes?\b|\bviri?ons?\b|\bvirus-like particles?\b|\bVLPs?\b"
        r"|\bviral (?:fraction|particles?|concentrates?|metagenomes?)\b",
        re.I,
    ),
    "single_cell": re.compile(r"\bsorted cells?\b|\bsorted cell/s\b", re.I),  # JGI's "Sorted cell/s from …": SAG data
    # an experimental incubation of natural water or sediment: fine for its untreated arm, not after
    # a treatment, which only the paper tells
    "manipulated": re.compile(r"\b(?:micro|meso)-?cosms?\b", re.I),
    # a culture, not the community itself; "culture-independent" says the opposite (NARMS)
    "culture": re.compile(r"\bcultur(?:e|es|ed|ing)\b(?![- ]independent)|\bnon-?axenic\b", re.I),
}
CO_ASSEMBLY_SAMPLE = re.compile(r"\bmock sample for co-?assembly\b", re.I)  # ENA's text for a virtual sample
CURATION_FLAGS = {"enrichment", "culture", "manipulated", "type_strain", "surveillance_isolate"}  # → tier C
EXCLUDING_FLAGS = {  # any of these sends a row to tier D
    *HARD_FILTERS,
    "out_of_scope_taxon",
    "excluded_envo",
    "off_biome",
    "amplicon_only",  # a run's library is amplicon (strategy AMPLICON, or selection PCR)
    "sample_low_reads",  # none of the sample's runs reaches 200k reads, so its reads can't be routed
} - CURATION_FLAGS
# Actions 1-5: biome evidence. Base weights are PROPOSAL.md §3.7's, still to calibrate.
W_THIRD_PARTY = 0.9  # action 1, a curated biome lineage (MGnify study)
UNINFORMATIVE_LINEAGES = {"root:Mixed"}  # MGnify's node for studies spanning several biomes
W_ENVO = 0.8  # action 2, an env_* term of the config (ENVO, subclasses included)
W_SAMPLE_KEYWORD = 0.4  # action 2, a config keyword in the sample's text: generic
W_SAMPLE_EXCLUDE = 0.6  # action 2, an exclude keyword in the sample's text: a conflict
W_STUDY = 0.6  # action 4, a config keyword in the study's title or description
W_STUDY_EXCLUDE = 0.3  # study texts often name other habitats in passing, so below the conflict threshold
W_METAGENOME = 0.3  # action 5, the sample's organism is one of the config's metagenome taxa: often a default
SIBLINGS = 0.1  # action 3, score modifier from the other samples of the same studies
CONFLICT = 0.6  # an opposing item this heavy sets needs_curation
TIER_A, TIER_B = 80, 60
ENV_FIELDS = {
    "env_broad_scale": "broad_scale_environmental_context",
    "env_local_scale": "local_environmental_context",
    "env_medium": "environmental_medium",
}
ENVO_ID = re.compile(r"ENVO[:_ ]?(\d{7,8})", re.I)


# --- search: samples.tsv ---

SAMPLES_COLUMNS = [  # PROPOSAL.md §3.6; columns are added, never renamed
    "record_accession",
    "data_type",
    "data_type_evidence",
    "assembly_level",
    "source_db",
    "sample_accession",
    "related_accessions",
    "study_accession",
    "secondary_study_accession",
    "bioproject",
    "reported_tax_id",
    "reported_tax_name",
    "tax_lineage",
    "domain",
    "taxdump_date",
    "publication",
    "paper_access",
    "paper_pdf",
    "reported_biome",
    "third_party_biome",
    "curated_biome",
    "score",
    "tier",
    "needs_curation",
    "evidence",
    "flags",
    "read_stats",
    "first_public",
    "curation_note",
    "decision",
    "destination",
    "local_path",
    "prepare_status",
]
KEPT_COLUMNS = ("curation_note", "decision", "destination", "local_path", "prepare_status")
MIN_RUN_READS = 200_000  # a sample's reads are routed only if one of its runs reaches this
REPORTED_BIOME = [  # (reported_biome key, ENA sample field)
    ("env_broad_scale", "broad_scale_environmental_context"),
    ("env_local_scale", "local_environmental_context"),
    ("env_medium", "environmental_medium"),
    ("isolation_source", "isolation_source"),
    ("organism", "scientific_name"),
    ("host", "host"),
    ("body_site", "host_body_site"),
]
SURVEILLANCE_PACKAGES = ("One Health Enteric", "Pathogen.")  # BioSample packages of public-health surveillance


# --- survey: tables and kingdom terms ---

KINGDOM_TERMS = {
    "prokaryotes": re.compile(r"bacteri|archae|prokaryot|\b16s\b", re.I),
    "eukaryotes": re.compile(
        r"eukaryot|protist|fung(?:i|al)|\balga|diatom|phytoplankton|dinoflagell|ciliat"
        r"|haptophyt|chlorophyt|oomycet|\byeast|amoeb|foraminifer|radiolar|\b18s\b",
        re.I,
    ),
}
SURVEY_COLUMNS = [
    "biome",
    "kingdom",
    "metagenome_taxon",
    "tax_id",
    "catalogue",
    "catalogue_updated",
    "catalogue_species",
    "catalogue_genomes",
    "catalogue_mags",
    "mags",
    "mags_new",
    "mags_mgnify",
    "mags_to_add",
    "mag_studies",
    "mgnify_bins",
    "mgnify_bins_new",
    "mgnify_bins_to_add",
    "assemblies",
    "assemblies_mgnify",
    "assemblies_kingdom",
    "assemblies_to_bin",
    "samples_to_bin",
    "assembly_studies",
    "runs",
    "runs_1m",
    "runs_1m_kingdom",
    "runs_1m_new",
    "runs_1m_kingdom_new",
    "runs_to_assemble",
    "samples_to_assemble",
    "runs_1m_paired",
    "runs_1m_long",
    "gbases_1m",
    "gbases_1m_kingdom",
    "gbases_to_assemble",
    "read_studies_1m",
]
RECORD_COLUMNS = [
    "data_type",
    "accession",
    "study_accession",
    "from_mgnify",
    "sample_accession",
    "metagenome_taxon",
    "first_public",
    "is_new",
    "read_count",
    "base_count",
    "library_layout",
    "instrument_platform",
    "center_name",
    "assembly_software",
    "kingdom_evidence",
    "kingdom_match",
    "kingdom_terms",
    "in_catalogue",
    "route",
    "covered_by",
    "scope_taxon",
    "tax_id",
    "scientific_name",
    "source_samples",
    "host",
    "country",
    "lat",
    "lon",
    "collection_date",
    "assembly_type",
]
STUDY_COLUMNS = [
    "study_accession",
    "from_mgnify",
    "metagenome_taxa",
    "first_public",
    "mags_to_add",
    "assemblies_to_bin",
    "samples_to_bin",
    "assemblies_with_mags",
    "mags_studies",
    "runs_to_assemble",
    "samples_to_assemble",
]


# --- Survey dashboard (report.py) ---

SUBMITTED, BY_MGNIFY = "MAGs (submitters)", "MAGs (MGnify)"
ASSEMBLIES, RUNS = "Assemblies", "Runs ≥1M"
KINDS = [SUBMITTED, BY_MGNIFY, ASSEMBLIES, RUNS]
ROUTES = {
    "MAGs to add": "→ genomes-catalogue-pipeline",
    "Samples to bin": "→ genomes-generation",
    "Samples to assemble": "→ miassembler, then genomes-generation",
}
TOP_STUDIES = 15
PROJECTS = (
    "Projects are ENA/INSDC BioProjects (PRJEB, PRJNA, PRJDB): one submission grouping the samples "
    "and data of one research effort, usually one paper. Counted per route: the projects holding "
    "the MAGs, the MAGs to add, and the assemblies or runs of the samples to bin or assemble. "
    "Many projects mean more varied data, and "
    "more to curate, since curation reads each project's paper and metadata. For MAGs (MGnify), "
    "the projects are MGnify's own analysis projects, not the original ones."
)
ACCENT = "#0072B2"
COLORS = {
    SUBMITTED: "#0072B2",
    BY_MGNIFY: "#56B4E9",
    ASSEMBLIES: "#E69F00",
    RUNS: "#009E73",
    "MAGs": "#0072B2",
    "MAGs to add": "#0072B2",
    "Assemblies to bin": "#E69F00",
    "Runs to assemble": "#009E73",
    "Samples to bin": "#E69F00",
    "Samples to assemble": "#009E73",
    "short reads": "#009E73",
    "long reads": "#CC79A7",
    "in catalogue": "#BBBBBB",
    "Bacteria": "#0072B2",
    "Archaea": "#D55E00",
    "Eukaryota": "#009E73",
}
# column: (title, tooltip)
REPORT_COLUMNS = {
    "biome": (
        "Catalogue",
        "The biome config: one catalogue per biome and kingdom, as in MGnify. *: MGnify already has "
        "this catalogue, so the survey counts an update",
    ),
    "metagenome_taxon": ("Metagenome taxon", "ENA metagenome taxon the records were found under"),
    "tax_id": ("Tax ID", "NCBI tax ID of the metagenome taxon"),
    "kingdom": ("Kingdom", "prokaryotes (Bacteria + Archaea) or eukaryotes"),
    "mags": (
        "MAGs",
        "MAGs of the catalogue's kingdom in ENA: deposited by their authors, or produced by MGnify from "
        "this biome's public data (ENA stores most of those as 'binned metagenome' records)",
    ),
    "to_add": (
        "MAGs to add",
        "MAGs (both kinds) not in the existing MGnify catalogue (all of them, for a new "
        "catalogue): ready for genomes-catalogue-pipeline",
    ),
    "samples_bin": (
        "Samples to bin",
        "Samples with a metagenome assembly but no genome yet: input for genomes-generation. Only "
        "assemblies with evidence for the kingdom count (a genome of that kingdom from the same sample, "
        "or a project that names it). A sample can have several assemblies, e.g. by its submitter and "
        "by MGnify; hover the routes chart for the assembly count",
    ),
    "samples_assemble": (
        "Samples to assemble",
        "Samples with only reads (no genome, no assembly): input for miassembler, then "
        "genomes-generation. Only metagenomic WGS runs with at least 1 million reads and evidence for "
        "the kingdom count. A sample can have several runs (lanes, re-sequencing); hover the routes "
        "chart for the run count",
    ),
    "tbases": ("Runs Tb", "Terabases (10^12 bases, not bytes) in the runs to assemble"),
}
FIGURE_LAYOUT = {
    "template": "plotly_white",
    "height": 380,
    "margin": {"l": 10, "r": 10, "t": 36, "b": 10},
    "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02, "x": 0, "title": None},
    "font": {"size": 12},
}
TABLE_CSS = """
:host .tabulator-col-title { white-space: normal; word-break: normal; overflow-wrap: normal; text-overflow: clip; }
:host .tabulator-col { font-size: 12px; }
:host .tabulator-col-sorter, :host .tabulator-arrow { display: none; }  /* click the title to sort */
"""
ABOUT = """
Every number is for the catalogue's kingdom: its genomes, and the assemblies and runs with
evidence for that kingdom. Counts are proxies: a record belongs to a catalogue through the
metagenome taxon its submitter chose (plus the host, for host-associated biomes), not through
curation. Hover a column title for its definition.

**Where new genomes would come from.** Each sample enters at its most processed public data:
*MAGs to add* are ready for the catalogue; *samples to bin* have assemblies that need binning
first; *samples to assemble* have only reads, which need assembly and binning. Many MAGs to add
means the catalogue can grow now; many samples to bin or assemble means it would pay to generate
new MAGs. Samples, not records, because a sample can have several assemblies or runs.
"""


# --- HTTP cache ---

CACHE_DIR = Path(os.environ.get("BIOME_FINDER_CACHE", "~/.cache/biome-finder")).expanduser()
CACHE_MAX_AGE = 7 * 24 * 3600  # ponytail: fixed one-week TTL, add a --refresh flag if stale data bites
USER_AGENT = "biome-finder/0.1"


# --- ENA ---

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api"
ENA_TAXONOMY = "https://www.ebi.ac.uk/ena/taxonomy/rest"
ENA_BROWSER = "https://www.ebi.ac.uk/ena/browser/api"
DERIVED_FROM = {"derived from", "derived_from", "sample derived from"}


# --- JGI GOLD export ---

GOLD_LEVELS = ["ECOSYSTEM", "ECOSYSTEM CATEGORY", "ECOSYSTEM TYPE", "ECOSYSTEM SUBTYPE", "SPECIFIC ECOSYSTEM"]


# --- NCBI ---

TAXDUMP_ARCHIVE = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive"
ASSEMBLY_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"


# --- OLS ---

OLS_API = "https://www.ebi.ac.uk/ols4/api/ontologies"


# --- MGnify ---

MGNIFY_API_V2 = "https://www.ebi.ac.uk/metagenomics/api/v2"
