"""Scoring: the curation actions of PROPOSAL.md §3.7, cheapest first."""

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from .models import EnaSampleMetadata, OlsTerm
from .constants import (
    CONFLICT,
    CO_ASSEMBLY_SAMPLE,
    ENVO_ID,
    ENV_FIELDS,
    HARD_FILTERS,
    SIBLINGS,
    TIER_A,
    TIER_B,
    UNINFORMATIVE_LINEAGES,
    W_ENVO,
    W_METAGENOME,
    W_SAMPLE_EXCLUDE,
    W_SAMPLE_KEYWORD,
    W_STUDY,
    W_STUDY_EXCLUDE,
    W_THIRD_PARTY,
)


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
    """The HARD_FILTERS a sample's text matches. An enrichment is a culture, so `enrichment` alone
    stands for both."""
    text = sample_text(s)
    found = [flag for flag, pattern in HARD_FILTERS.items() if pattern.search(text)]
    return [f for f in found if f != "culture" or "enrichment" not in found]


@dataclass(frozen=True)
class Evidence:
    """One evidence item (CLAUDE.md "Scoring conventions"). `claim` is `in` (the config's biome)
    or `out` (an excluded or other biome); `quote` is the text that matched."""

    action: str
    source: str
    claim: str
    specificity: float
    weight: float
    quote: str

    def __str__(self) -> str:
        return f"{self.action}:{self.source}:{self.claim}:{self.quote}"


@dataclass(frozen=True)
class Terms:
    """A config's ENVO terms with their subclasses, as IDs and as lowercase labels and synonyms,
    since most submitters write the name without the ID."""

    ids: frozenset[str]
    labels: frozenset[str]
    synonyms: frozenset[str]

    @classmethod
    def of(cls, terms: Iterable[OlsTerm]) -> "Terms":
        terms = list(terms)
        return cls(
            frozenset(t.obo_id for t in terms),
            frozenset(t.label.lower() for t in terms),
            frozenset(n.lower() for t in terms for n in t.synonyms),
        )

    def match(self, value: str) -> int:
        """How strongly `value` (an env_* field) names one of the terms: 3 by ENVO ID, 2 by exact
        label once any bracketed ID is removed ("freshwater lake biome [ENVO:01000252]"), 1 by
        exact synonym, 0 not at all. Synonyms are weakest because ENVO shares them: "lake" is a
        synonym of pond, reservoir and saline evaporation pond."""
        if any(f"ENVO:{m.zfill(8)}" in self.ids for m in ENVO_ID.findall(value)):
            return 3
        name = re.sub(r"\[[^]]*\]", "", value).strip().lower()
        return 2 if name in self.labels else 1 if name in self.synonyms else 0


def keyword(words: Iterable[str]) -> re.Pattern[str] | None:
    """Config keywords as one pattern; each matches at a word start ("estuar" → "estuarine")."""
    words = [w for w in words if w]
    return re.compile(r"\b(?:" + "|".join(re.escape(w) for w in words) + ")", re.I) if words else None


def without(text: str, names: Iterable[str]) -> str:
    """`text` without these names (case-insensitive) and ENA's co-assembly boilerplate: the
    metagenome taxa, which default titles repeat ("Metagenome or environmental sample from
    freshwater metagenome") and action 5 already counts, and "mock sample for coassembly",
    which is no mock community."""
    text = CO_ASSEMBLY_SAMPLE.sub(" ", text)
    for n in names:
        text = re.sub(re.escape(n), " ", text, flags=re.I)
    return text


def sample_items(
    s: EnaSampleMetadata,
    include: Terms,
    exclude: Terms,
    kw_in: re.Pattern[str] | None,
    kw_out: re.Pattern[str] | None,
    taxa: set[str],
) -> list[Evidence]:
    """Actions 2 and 5 for one sample: its env_* terms, the config keywords in its own text, and
    its metagenome taxon (`taxa`: the config's, lowercase). The keywords skip the organism and
    `taxa` names, so a default title doesn't count twice. An env_* term under `envo_exclude` gives
    an `excluded_envo` item (the sample is ruled out, CLAUDE.md) when it matches at least as
    strongly as an included one (`Terms.match`)."""
    items = []
    for key, field in ENV_FIELDS.items():
        value = getattr(s, field)
        out, in_ = (exclude.match(value), include.match(value)) if value else (0, 0)
        if out and out >= in_:
            items.append(Evidence("sample_metadata", key, "excluded_envo", 1.0, 1.0, value))
        elif in_:
            items.append(Evidence("sample_metadata", key, "in", 1.0, W_ENVO, value))
    text = without(sample_text(s), {*taxa, s.scientific_name})
    if kw_out and (m := kw_out.search(text)):
        items.append(Evidence("sample_metadata", "text", "out", 1.0, W_SAMPLE_EXCLUDE, m.group(0)))
    if kw_in and (m := kw_in.search(text)):
        items.append(Evidence("sample_metadata", "text", "in", 1.0, W_SAMPLE_KEYWORD, m.group(0)))
    organism = s.scientific_name.lower()
    if organism in taxa:
        items.append(Evidence("metagenome_name", "organism", "in", 1.0, W_METAGENOME, s.scientific_name))
    elif organism.endswith("metagenome"):
        items.append(Evidence("metagenome_name", "organism", "out", 1.0, W_METAGENOME, s.scientific_name))
    return items


def study_items(
    study: str, text: str, kw_in: re.Pattern[str] | None, kw_out: re.Pattern[str] | None, taxa: set[str]
) -> list[Evidence]:
    """Action 4: the config keywords in a study's title and description, skipping the `taxa` names."""
    text = without(text, taxa)
    items = []
    if kw_in and (m := kw_in.search(text)):
        items.append(Evidence("project_description", study, "in", 1.0, W_STUDY, m.group(0)))
    if kw_out and (m := kw_out.search(text)):
        items.append(Evidence("project_description", study, "out", 1.0, W_STUDY_EXCLUDE, m.group(0)))
    return items


def under(path: str, parent: str, sep: str = ":") -> bool:
    """`path` is `parent` or one of its sub-paths."""
    return path == parent or path.startswith(parent + sep)


def tree_claim(path: str, targets: list[str], excludes: list[str], sep: str) -> tuple[str, float]:
    """How a curated biome path relates to the config's: (`in` or `out`, specificity). Under one
    of `targets` and none of `excludes`: in, 1. An ancestor of a target (e.g. `…:Freshwater` for
    `…:Freshwater:Lake`) supports it only partly: in, its depth over the target's. Anything else,
    excluded subtrees included: out, 1."""
    if any(under(path, x, sep) for x in excludes):
        return "out", 1.0
    if any(under(path, t, sep) for t in targets):
        return "in", 1.0
    depth = len(path.split(sep))
    partial = [depth / len(t.split(sep)) for t in targets if under(t, path, sep)]
    return ("in", max(partial)) if partial else ("out", 1.0)


def lineage_claim(lineage: str, cfg: Mapping[str, Any]) -> tuple[str, float]:
    """tree_claim for an MGnify lineage against `mgnify_lineages` / `mgnify_lineages_exclude`,
    after mapping it through `lineage_equivalents` (a lineage and its sub-lineages)."""
    for src, dst in (cfg.get("lineage_equivalents") or {}).items():
        if under(lineage, src):
            lineage = dst + lineage[len(src) :]
    return tree_claim(lineage, cfg.get("mgnify_lineages") or [], cfg.get("mgnify_lineages_exclude") or [], ":")


def gold_claim(path: str, cfg: Mapping[str, Any]) -> tuple[str, float]:
    """tree_claim for a GOLD ecosystem path against `gold_ecosystems` / `gold_ecosystems_exclude`.
    Trailing `Unclassified` levels say nothing, so they are dropped first."""
    levels = path.split(" > ")
    while levels and levels[-1] in ("Unclassified", ""):
        levels.pop()
    return tree_claim(
        " > ".join(levels), cfg.get("gold_ecosystems") or [], cfg.get("gold_ecosystems_exclude") or [], " > "
    )


def third_party_items(source: str, lineage: str, cfg: Mapping[str, Any]) -> list[Evidence]:
    """Action 1: a curated biome lineage (e.g. the MGnify study's) as one evidence item; none
    for an UNINFORMATIVE_LINEAGES one, which says nothing about a single sample."""
    if not lineage or lineage in UNINFORMATIVE_LINEAGES:
        return []
    claim, specificity = lineage_claim(lineage, cfg)
    return [Evidence("third_party_biome", source, claim, round(specificity, 2), W_THIRD_PARTY, lineage)]


def gold_items(gold_id: str, path: str, cfg: Mapping[str, Any]) -> list[Evidence]:
    """Action 1: the sample's GOLD ecosystem as one evidence item, if the config has GOLD paths."""
    if not path or not cfg.get("gold_ecosystems"):
        return []
    claim, specificity = gold_claim(path, cfg)
    return [Evidence("third_party_biome", f"gold:{gold_id}", claim, round(specificity, 2), W_THIRD_PARTY, path)]


def combine(items: list[Evidence], siblings: float = 0.0) -> tuple[float, bool]:
    """Score 0-100 and conflict for one sample. Agreeing items combine as 1 − Π(1 − w·s);
    opposing ones scale that down the same way; `siblings` (−1..1) adds up to ±SIBLINGS. An
    opposing item of weight ≥ CONFLICT is a conflict."""

    def pooled(claim: str) -> float:
        p = 1.0
        for e in items:
            if e.claim == claim:
                p *= 1 - e.weight * e.specificity
        return 1 - p

    score = pooled("in") * (1 - pooled("out")) + SIBLINGS * siblings
    conflict = any(e.claim == "out" and e.weight >= CONFLICT for e in items)
    return round(100 * min(max(score, 0.0), 1.0)), conflict


def off_biome(sample_items: list[Evidence]) -> bool:
    """The sample's own evidence (its metadata, organism and GOLD ecosystem: not its study's)
    points to another biome at conflict weight, and none of it points to the config's at that
    weight: the sample belongs to another biome, whatever its study says, since a study can mix
    biomes. A generic keyword or a default metagenome name doesn't rescue it (e.g. GOLD
    `… > Lake > Sediment` against "freshwater" in the isolation source); a strong item both
    ways (an ENVO term against GOLD) is a conflict for a curator.
    A keyword in free text never excludes on its own (it can be incidental, or a manipulation
    such as a microcosm or a culture, which a curator checks in the paper): it stays a conflict,
    tier C (`off_biome_keyword`)."""
    against = [e for e in sample_items if e.claim in ("out", "excluded_envo") and e.source != "text"]
    out = max((e.weight for e in against), default=0.0)
    in_ = max((e.weight for e in sample_items if e.claim == "in"), default=0.0)
    return out >= CONFLICT > in_


def tier(score: float, conflict: bool) -> str:
    """A ≥ TIER_A with no conflict, B ≥ TIER_B with no conflict, C otherwise."""
    if conflict:
        return "C"
    return "A" if score >= TIER_A else "B" if score >= TIER_B else "C"
