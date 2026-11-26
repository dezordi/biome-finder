"""Typed queries and responses for every third-party service biome-finder calls, by service.

Each client in `clients/` takes a query model (or a plain identifier, when the query is a
single one) and returns `Response[T]`: the HTTP status plus the data models its callers use.
A data model lists only the fields biome-finder uses; for ENA, those are also the fields
requested from the Portal API.

The last section holds biome-finder's own tables (a survey record, a samples.tsv row), so
commands pass typed rows between their steps.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field, fields
from types import NoneType
from typing import Any, ClassVar, Generic, Protocol, TypeVar, get_args, get_type_hints
from .constants import NCBI_RANKS


class DataModel(Protocol):
    """Any dataclass."""

    __dataclass_fields__: ClassVar[dict[str, Any]]


T = TypeVar("T")
M = TypeVar("M", bound=DataModel)


@dataclass(frozen=True)
class Response(Generic[T]):
    """`status`: the HTTP status, 200 for a cached response too, 404 for "not found" (any other
    error raises). A response built from several requests (pages, batches) carries the first
    status that isn't 200. `data`: what the call returns."""

    status: int
    data: T

    @property
    def ok(self) -> bool:
        return self.status == 200


def column(name: str) -> Any:
    """A field read from a source column with another name (e.g. MGnify's `Genome_type`)."""
    return field(metadata={"column": name})


def columns(model: type[DataModel]) -> list[str]:
    """The source columns of a data model, in field order."""
    return [f.metadata.get("column", f.name) for f in fields(model)]


def from_row(model: type[M], row: Mapping[str, Any]) -> M:
    """A data model from one source row (TSV or JSON). A missing column is empty; `bool` ("True" /
    "False"), `int` and `float` fields are converted, and an empty value becomes None when the
    field allows it (False for a plain `bool`)."""
    hints = get_type_hints(model)
    values = {}
    for f in fields(model):
        value, hint = row.get(f.metadata.get("column", f.name)), hints[f.name]
        kinds = get_args(hint) or (hint,)
        if value in ("", None):
            value = None if NoneType in kinds else False if bool in kinds else ""
        elif bool in kinds:
            value = value in (True, "True", "true")
        elif int in kinds:
            value = int(value)
        elif float in kinds:
            value = float(value)
        values[f.name] = value
    return model(**values)


# --- ENA Portal API (https://www.ebi.ac.uk/ena/portal/api) ---


@dataclass(frozen=True)
class EnaSearch:
    """One Portal API search. `query` is in the portal's syntax (`tax_eq`, `tax_tree`, …);
    `accessions` restricts it to these primary IDs (sent as `includeAccessions`); `fields`
    overrides the row model's fields, e.g. to fetch only the accessions."""

    RESULTS: ClassVar[set[str]] = {"wgs_set", "analysis", "read_run", "sample", "study", "assembly"}

    result: str
    query: str = ""
    accessions: tuple[str, ...] = ()
    fields: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.result not in self.RESULTS:
            raise ValueError(f"unknown ENA result {self.result!r}; expected one of {sorted(self.RESULTS)}")
        if not self.query and not self.accessions:  # would fetch the whole result table
            raise ValueError(f"ENA {self.result} search needs a query or accessions")


@dataclass(frozen=True)
class EnaMag:
    """A MAG (`wgs_set`, assembly_type "metagenome-assembled genome (mag)")."""

    accession: str
    study_accession: str
    sample_accession: str
    first_public: str  # YYYY-MM-DD
    tax_id: str
    scientific_name: str
    country: str
    collection_date: str


@dataclass(frozen=True)
class EnaAssembly:
    """A metagenome assembly (`analysis`, assembly_type "primary metagenome")."""

    analysis_accession: str
    study_accession: str
    sample_accession: str
    first_public: str
    center_name: str
    assembly_software: str
    country: str
    lat: str  # decimal degrees; ';'-joined when ENA merges several samples
    lon: str
    collection_date: str


@dataclass(frozen=True)
class EnaBin:
    """A bin (`analysis`, assembly_type "binned metagenome"), how MGnify deposits its genomes."""

    analysis_accession: str
    study_accession: str
    sample_accession: str
    first_public: str
    tax_id: str
    scientific_name: str


def is_community(assembly_type: str, scientific_name: str) -> bool:
    """An unbinned metagenome assembly rather than a genome: typed `primary metagenome`, or
    untyped with a metagenome as its organism (how NCBI-submitted ones look in ENA)."""
    kind = assembly_type.lower()
    return kind == "primary metagenome" or (not kind and "metagenome" in scientific_name.lower())


@dataclass(frozen=True)
class EnaLinkedWgs:
    """A WGS set linked to a sample, of any type (`wgs_set`)."""

    accession: str
    study_accession: str
    sample_accession: str
    first_public: str
    assembly_type: str  # "Metagenome-Assembled Genome (MAG)", empty, …
    tax_id: str
    scientific_name: str  # the organism, or the metagenome itself
    country: str
    collection_date: str

    @property
    def community(self) -> bool:
        return is_community(self.assembly_type, self.scientific_name)


@dataclass(frozen=True)
class EnaLinkedAnalysis:
    """An assembly analysis linked to a sample, of any type (`analysis`, SEQUENCE_ASSEMBLY)."""

    accession: str = column("analysis_accession")
    study_accession: str = column("study_accession")
    sample_accession: str = column("sample_accession")
    first_public: str = column("first_public")
    assembly_type: str = column("assembly_type")
    tax_id: str = column("tax_id")
    scientific_name: str = column("scientific_name")
    center_name: str = column("center_name")
    assembly_software: str = column("assembly_software")
    country: str = column("country")
    lat: str = column("lat")
    lon: str = column("lon")
    collection_date: str = column("collection_date")

    @property
    def community(self) -> bool:
        return is_community(self.assembly_type, self.scientific_name)


@dataclass(frozen=True)
class EnaRun:
    """A read run (`read_run`)."""

    run_accession: str
    study_accession: str
    sample_accession: str
    first_public: str
    read_count: int | None
    base_count: int | None
    library_layout: str  # PAIRED, SINGLE
    instrument_platform: str  # ILLUMINA, OXFORD_NANOPORE, PACBIO_SMRT, …
    country: str
    lat: str  # decimal degrees; ';'-joined when ENA merges several samples
    lon: str
    collection_date: str


@dataclass(frozen=True)
class EnaGenome:
    """A single-organism genome assembly (`assembly`: GCA), with the sample fields that tell its
    data type. NCBI-submitted MAGs come untyped here, so `assembly_type` alone can't tell a MAG
    from an isolate: the BioSample package and `environmental_sample` do (PROPOSAL.md §3.4)."""

    accession: str  # GCA_…, versionless
    assembly_type: str  # mostly empty; "clone or isolate", …
    assembly_level: str  # contig, scaffold, chromosome, complete genome
    environmental_sample: str  # "true" for metagenome-derived
    ncbi_reporting_standard: str  # BioSample package: MIMAG.*, MISAG.*, MIGS.ba.*, Microbe.1.0, …
    strain: str
    culture_collection: str
    isolate: str  # free text, often a MAG or SAG name (e.g. "UBA10892")
    sample_accession: str
    study_accession: str
    tax_id: str
    scientific_name: str
    first_created: str
    description: str


@dataclass(frozen=True)
class EnaRunLibrary:
    """How a run's library was made (`read_run`): amplicon data hides under strategy OTHER with
    a PCR selection, which the survey's WGS/OTHER query lets through."""

    run_accession: str
    library_strategy: str  # WGS, OTHER, AMPLICON, …
    library_selection: str  # RANDOM, PCR, size fractionation, …


@dataclass(frozen=True)
class EnaSample:
    """A sample's IDs and host (`sample`)."""

    sample_accession: str  # primary: SAMEA/SAMN/SAMD
    secondary_sample_accession: str  # ERS/SRS/DRS
    host: str
    host_tax_id: str

    @property
    def ids(self) -> list[str]:
        return [a for a in (self.sample_accession, self.secondary_sample_accession) if a]


@dataclass(frozen=True)
class EnaSampleMetadata:
    """A sample's descriptive fields (`sample`), the raw biome evidence `search` reads. ENA's
    `environment_biome/feature/material` repeat the three MIxS env_* fields, so only those are
    requested."""

    sample_accession: str  # primary: SAMEA/SAMN/SAMD
    secondary_sample_accession: str  # ERS/SRS/DRS
    sample_title: str
    sample_description: str
    scientific_name: str  # the organism; a metagenome taxon for community samples
    broad_scale_environmental_context: str  # MIxS env_broad_scale
    local_environmental_context: str  # env_local_scale
    environmental_medium: str  # env_medium
    isolation_source: str
    host: str
    host_body_site: str
    lat: str
    lon: str
    country: str

    @property
    def ids(self) -> list[str]:
        return [a for a in (self.sample_accession, self.secondary_sample_accession) if a]


@dataclass(frozen=True)
class EnaStudy:
    """A study's accession and text (`study`)."""

    study_accession: str
    study_title: str
    study_description: str

    @property
    def text(self) -> str:
        return f"{self.study_title} {self.study_description}"


# --- ENA Taxonomy REST (https://www.ebi.ac.uk/ena/taxonomy/rest) ---


@dataclass(frozen=True)
class EnaTaxon:
    tax_id: str = column("taxId")
    scientific_name: str = column("scientificName")


# --- MGnify API v2 (https://www.ebi.ac.uk/metagenomics/api/v2): genome catalogues ---


@dataclass(frozen=True)
class MgnifyCatalogue:
    catalogue_id: str
    catalogue_type: str  # prokaryotes | eukaryotes
    updated_at: str
    genome_count: int  # species
    unclustered_genome_count: int  # all genomes
    ftp_url: str


@dataclass(frozen=True)
class MgnifyStudy:
    """A study MGnify analysed (`/studies/`), with the biome its curators gave it."""

    accession: str  # MGYS…
    ena_accessions: tuple[str, ...]  # the ENA study's PRJ and ERP/SRP/DRP IDs
    lineage: str  # e.g. root:Environmental:Aquatic:Freshwater:Lake


# --- MGnify FTP: a catalogue's genomes-all_metadata.tsv ---


@dataclass(frozen=True)
class MgnifyCatalogueGenome:
    genome: str = column("Genome")
    genome_type: str = column("Genome_type")  # MAG | Isolate
    genome_accession: str = column("Genome_accession")
    sample_accession: str = column("Sample_accession")


# --- NCBI Taxonomy: monthly new_taxdump archive (ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive) ---


@dataclass(frozen=True)
class NcbiTaxon:
    """A taxon of one taxdump release (`nodes.dmp`, `rankedlineage.dmp`, `taxidlineage.dmp`).
    Each rank field holds the name at that rank, the taxon's own rank included, and None where
    the lineage has no taxon at that rank (e.g. no kingdom for most Bacteria). For viruses,
    `domain` is their realm."""

    tax_id: str
    name: str
    rank: str  # the taxon's own rank: species, genus, strain, no rank, …
    domain: str | None
    kingdom: str | None
    phylum: str | None
    class_: str | None  # `class` is a Python keyword
    order: str | None
    family: str | None
    genus: str | None
    species: str | None
    ancestors: frozenset[str]  # tax IDs above it

    @property
    def ranks(self) -> dict[str, str | None]:
        """NCBI_RANKS → name or None, root first."""
        return {r: getattr(self, "class_" if r == "class" else r) for r in NCBI_RANKS}


@dataclass(frozen=True)
class NcbiAssembly:
    """A GenBank assembly's row in NCBI's assembly_summary_genbank.txt, with the fields that
    tell its data type and type material (PROPOSAL.md §3.4)."""

    accession: str  # GCA_…, versionless
    biosample: str
    excluded_from_refseq: str  # e.g. "derived from metagenome", "derived from single cell"
    relation_to_type_material: str  # e.g. "assembly from type material"


@dataclass(frozen=True)
class NcbiTaxdump:
    """The part of a taxdump release a set of tax IDs needs."""

    release: str  # YYYY-MM-01
    merged: dict[str, str]  # old tax ID → current one (`merged.dmp`)
    deleted: frozenset[str]  # `delnodes.dmp`, among the requested IDs
    taxa: dict[str, NcbiTaxon]  # by current tax ID


# --- OLS4 (https://www.ebi.ac.uk/ols4/api): ontology terms (ENVO) ---


@dataclass(frozen=True)
class OlsTerm:
    """An ontology term with the names a submitter may write it as."""

    obo_id: str  # e.g. ENVO:00000021
    label: str
    synonyms: tuple[str, ...]


# --- JGI GOLD: the public data export (goldData.xlsx), downloaded by the curator ---


@dataclass(frozen=True)
class GoldBiosample:
    """A BioSample's GOLD ecosystem: from a sequencing project's NCBI BioSample to its GOLD
    Biosample (metagenomes) or Organism (isolates)."""

    biosample_accession: str  # SAMN/SAMEA/SAMD
    gold_id: str  # Gb… or Go…
    ecosystem_path: str  # "Environmental > Aquatic > Freshwater > Lake > Unclassified"


# --- biome-finder's own tables ---


class Table:
    """A row of one of biome-finder's own TSV tables: `to_tsv` writes it as text, column by
    column (a field with `column` metadata None isn't written)."""

    def to_tsv(self) -> dict[str, str]:
        out = {}
        for f in fields(self):  # type: ignore[arg-type]
            if f.metadata.get("column", f.name) is None:
                continue
            value = getattr(self, f.name)
            if isinstance(value, list):
                value = ";".join(value)
            out[f.name] = "" if value is None else str(value)
        return out


@dataclass
class SurveyRecord(Table):
    """One survey record, a row of `biomes/<name>/records.tsv` (RECORD_COLUMNS): an ENA record
    plus the survey's own columns, filled step by step."""

    data_type: str  # mag, bin, contigs, reads
    accession: str
    study_accession: str = ""
    from_mgnify: bool = False  # in one of MGnify's own ENA studies
    sample_accession: str = ""  # the record's own sample (a MAG's sample, not its source)
    metagenome_taxon: str = ""
    first_public: str = ""
    is_new: bool | None = None  # public on or after the catalogue cutoff; None without a catalogue
    read_count: int | None = None
    base_count: int | None = None
    library_layout: str = ""
    instrument_platform: str = ""
    center_name: str = ""
    assembly_software: str = ""
    kingdom_evidence: str = ""  # genome, study_text, none, other_kingdom (assemblies and runs)
    kingdom_match: bool | None = None  # None for genomes
    kingdom_terms: str = ""  # `;`-joined
    in_catalogue: bool | None = None  # genomes only
    route: str = ""  # bin, assemble, covered (assemblies and runs)
    covered_by: str = ""  # `;`-joined studies
    scope_taxon: str = ""
    tax_id: str = ""
    scientific_name: str = ""
    source_samples: str = ""  # `;`-joined
    host: str = ""
    country: str = ""
    lat: str = ""
    lon: str = ""
    collection_date: str = ""
    assembly_type: str = ""


@dataclass
class SurveyCounts:
    """The counts of one set of survey records (`survey.summarize`)."""

    mags: int
    mags_new: int | None
    mags_mgnify: int
    mags_to_add: int
    mag_studies: int
    mgnify_bins: int
    mgnify_bins_new: int | None
    mgnify_bins_to_add: int
    assemblies: int
    assemblies_mgnify: int
    assemblies_kingdom: int
    assemblies_to_bin: int
    samples_to_bin: int
    assembly_studies: int
    runs: int
    runs_1m: int
    runs_1m_kingdom: int
    runs_1m_new: int | None
    runs_1m_kingdom_new: int | None
    runs_to_assemble: int
    samples_to_assemble: int
    runs_1m_paired: int
    runs_1m_long: int
    gbases_1m: float
    gbases_1m_kingdom: float
    gbases_to_assemble: float
    read_studies_1m: int


@dataclass
class SurveyRow(Table):
    """One row of survey.tsv (SURVEY_COLUMNS): a metagenome taxon of a biome, or its TOTAL."""

    biome: str
    kingdom: str
    metagenome_taxon: str
    tax_id: str
    catalogue: str
    catalogue_updated: str
    catalogue_species: int | None
    catalogue_genomes: int | None
    catalogue_mags: int | None
    counts: SurveyCounts = field(metadata={"column": None})

    def to_tsv(self) -> dict[str, str]:
        counts = {k: "" if v is None else str(v) for k, v in self.counts.__dict__.items()}
        return {**super().to_tsv(), **counts}


@dataclass
class StudyRow(Table):
    """One row of a biome's studies.tsv (STUDY_COLUMNS): a study's records on a route."""

    study_accession: str
    from_mgnify: bool
    metagenome_taxa: str  # `;`-joined
    first_public: str
    mags_to_add: int
    assemblies_to_bin: int
    samples_to_bin: int
    assemblies_with_mags: int
    mags_studies: str  # `;`-joined
    runs_to_assemble: int
    samples_to_assemble: int


@dataclass
class SampleRow(Table):
    """One samples.tsv row (PROPOSAL.md §3.6). Its fields are the table's columns in order
    (SAMPLES_COLUMNS), plus `read_count`, which search uses and doesn't write. `sample_accession`,
    `flags` and `evidence` are lists, `;`-joined in the table (`; ` for evidence)."""

    record_accession: str
    data_type: str
    data_type_evidence: str = ""
    assembly_level: str = ""
    source_db: str = ""
    sample_accession: list[str] = field(default_factory=list)
    related_accessions: str = ""
    study_accession: str = ""
    secondary_study_accession: str = ""
    bioproject: str = ""
    reported_tax_id: str = ""
    reported_tax_name: str = ""
    tax_lineage: str = ""
    domain: str = ""
    taxdump_date: str = ""
    publication: str = ""
    paper_access: str = ""
    paper_pdf: str = ""
    reported_biome: str = ""
    third_party_biome: str = ""
    curated_biome: str = ""
    score: int | None = None
    tier: str = ""
    needs_curation: bool | None = None
    evidence: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    read_stats: str = ""
    first_public: str = ""
    curation_note: str = ""
    decision: str = ""
    destination: str = ""
    local_path: str = ""
    prepare_status: str = ""
    read_count: int = field(default=0, metadata={"column": None})  # not written

    def flag(self, *names: str) -> None:
        """Add flags, once each, keeping their order."""
        self.flags += [n for n in names if n not in self.flags]

    def to_tsv(self) -> dict[str, str]:
        """The row as samples.tsv text, column by column."""
        out = {}
        for f in fields(self):
            if f.metadata.get("column", f.name) is None:
                continue
            value = getattr(self, f.name)
            if isinstance(value, list):
                value = ("; " if f.name == "evidence" else ";").join(value)
            out[f.name] = "" if value is None else str(value)
        return out

    @classmethod
    def from_tsv(cls, row: Mapping[str, str]) -> "SampleRow":
        """A row read back from samples.tsv (e.g. to merge a curator's decisions)."""
        columns = {f.name for f in fields(cls) if f.metadata.get("column", f.name) is not None}
        values: dict[str, Any] = {k: v for k, v in row.items() if k in columns}
        for name, sep in (("sample_accession", ";"), ("flags", ";"), ("evidence", "; ")):
            values[name] = [x for x in values.get(name, "").split(sep) if x]
        values["score"] = int(values["score"]) if values.get("score") else None
        values["needs_curation"] = {"True": True, "False": False}.get(values.get("needs_curation", ""))
        return cls(**values)
