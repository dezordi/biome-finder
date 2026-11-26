"""Typed queries and responses for every third-party service biome-finder calls, by service.

Each client in `clients/` takes a query model (or a plain identifier, when the query is a
single one) and returns `Response[T]`: the HTTP status plus the data models its callers use.
A data model lists only the fields biome-finder uses; for ENA, those are also the fields
requested from the Portal API.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field, fields
from types import NoneType
from typing import Any, ClassVar, Generic, Protocol, TypeVar, get_args, get_type_hints


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
    """A data model from one source row (TSV or JSON). A missing column is empty; `int` and
    `float` fields are converted, and an empty value becomes None when the field allows it."""
    hints = get_type_hints(model)
    values = {}
    for f in fields(model):
        value, hint = row.get(f.metadata.get("column", f.name)), hints[f.name]
        kinds = get_args(hint) or (hint,)
        if value in ("", None):
            value = None if NoneType in kinds else ""
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

    RESULTS: ClassVar[set[str]] = {"wgs_set", "analysis", "read_run", "sample", "study"}

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


# --- MGnify FTP: a catalogue's genomes-all_metadata.tsv ---


@dataclass(frozen=True)
class MgnifyCatalogueGenome:
    genome: str = column("Genome")
    genome_type: str = column("Genome_type")  # MAG | Isolate
    genome_accession: str = column("Genome_accession")
    sample_accession: str = column("Sample_accession")


# --- NCBI Taxonomy: monthly new_taxdump archive (ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive) ---


NCBI_RANKS = ("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")


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
