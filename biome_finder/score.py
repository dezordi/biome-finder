"""Scoring: the curation actions of PROPOSAL.md §3.7, cheapest first."""

import re

from .models import EnaSampleMetadata

# Action 0, hard filters: a sample whose own text names one of these is not a natural sample of
# any biome. Patterns were checked against freshwater-water-prok (2026-11): bare "control" there
# is an experiment's untreated arm (a real sample, left to the biome evidence), and "mock" is
# mostly ENA's "mock sample for coassembly", the virtual sample of a co-assembly.
HARD_FILTERS = {
    "control": re.compile(
        r"\b(?:negative|blank|extraction|kit|reagent|no[- ]template)[ ._-]?controls?\b|\bblank|\bntc\b", re.I
    ),
    "mock": re.compile(r"\bmock\b(?! sample for co-?assembly)", re.I),
    "synthetic": re.compile(r"\bsynthetic\b", re.I),
    "enrichment": re.compile(r"\benrichments?\b|\benriched (?:microbiome|culture|communit)", re.I),
    "spiked": re.compile(r"\bspiked?\b|\bspike-?ins?\b", re.I),
}
EXCLUDING_FLAGS = {*HARD_FILTERS, "out_of_scope_taxon"}  # any of these sends a row to tier X


def sample_text(s: EnaSampleMetadata) -> str:
    """The sample's free text and environment fields, where its submitter describes it."""
    return " ".join(
        [
            s.sample_title,
            s.sample_description,
            s.isolation_source,
            s.broad_scale_environmental_context,
            s.local_environmental_context,
            s.environmental_medium,
        ]
    )


def hard_flags(s: EnaSampleMetadata) -> list[str]:
    """The HARD_FILTERS a sample's text matches."""
    text = sample_text(s)
    return [flag for flag, pattern in HARD_FILTERS.items() if pattern.search(text)]
