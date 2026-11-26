"""ENA: Portal API (searches), Browser API (sample XML) and Taxonomy REST."""

import csv
import io
import json
import urllib.parse
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from typing import Any

from ..log import warn
from ..models import (
    M,
    EnaGenome,
    EnaLinkedAnalysis,
    EnaLinkedWgs,
    EnaRunLibrary,
    EnaSample,
    EnaSampleMetadata,
    EnaSearch,
    EnaStudy,
    EnaTaxon,
    Response,
    columns,
    from_row,
)
from ..constants import DERIVED_FROM, ENA_BROWSER, ENA_PORTAL, ENA_TAXONOMY, SAMPLE_ID
from .http import batches, first_status, get

csv.field_size_limit(10**9)


# --- Portal API ---


def count(request: EnaSearch) -> int:
    """Rows a query would return. A long query (e.g. hundreds of sample IDs) goes by POST, since it
    would not fit in a URL; short ones stay GET, keeping their cache keys."""
    form: dict[str, str] = {"result": request.result, "query": request.query}
    res = get(f"{ENA_PORTAL}/count", data=dict(form)) if len(request.query) > 1500 else get(f"{ENA_PORTAL}/count", form)
    text = res.data
    return int(text.split()[-1]) if text else 0


def parse_tsv(text: str, fields: list[str]) -> list[dict[str, str]] | None:
    """Rows of an ENA TSV page; None if the page is truncated (ENA ends it with 'Premature EOF')."""
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if any(None in r.values() or None in r for r in rows) or "Premature EOF" in text[-200:]:
        return None
    return rows


def search(request: EnaSearch, model: type[M]) -> Response[list[M]]:
    """All rows of an ENA Portal search, as `model`s. Uses POST with limit=0 (large GET
    responses get truncated with 'Premature EOF'). A query's row total is checked against the
    count endpoint: a mismatch means records published or suppressed between the two calls, or
    a silent cut. The form's key order is part of the cache key, so keep it stable."""
    fields = list(request.fields) or columns(model)
    form: dict[str, object] = {"result": request.result}
    if request.query:
        form["query"] = request.query
    if request.accessions:
        form["includeAccessions"] = ",".join(request.accessions)
    form |= {"fields": ",".join(fields), "format": "tsv", "limit": 0}
    total = None if request.accessions else count(request)
    if total == 0:
        return Response(200, [])
    res = get(f"{ENA_PORTAL}/search", data=form, check=lambda t: parse_tsv(t, fields) is not None)
    rows = (parse_tsv(res.data, fields) if res.data else None) or []
    if total is not None and len(rows) != total:
        warn(f"ENA {request.result}: {len(rows)} rows, count endpoint says {total} ({request.query})")
    return Response(res.status, [from_row(model, r) for r in rows])


def samples(accessions: Iterable[str], model: type[M]) -> Response[list[M]]:
    """Samples by primary or secondary ID, as `model`. Primary IDs go through `includeAccessions`,
    which rejects secondary (ERS/SRS/DRS) ones; those go through an OR query."""
    accessions = set(accessions)
    primary = {a for a in accessions if a.startswith("SAM")}
    responses = [search(EnaSearch("sample", accessions=tuple(b)), model) for b in batches(primary, 500)]
    responses += [
        search(EnaSearch("sample", " OR ".join(f'secondary_sample_accession="{a}"' for a in b)), model)
        for b in batches(accessions - primary, 100)
    ]
    return Response(first_status(responses), [s for r in responses for s in r.data])


def genomes(query: str) -> Response[list[EnaGenome]]:
    """Genome assemblies (`assembly`) matching a Portal query, one per accession (the result has
    a row per assembly version)."""
    res = search(EnaSearch("assembly", query), EnaGenome)
    return Response(res.status, list({g.accession: g for g in reversed(res.data)}.values()))


def run_libraries(accessions: Iterable[str]) -> Response[list[EnaRunLibrary]]:
    """The library strategy and selection of each run, 500 runs per request."""
    responses = [search(EnaSearch("read_run", accessions=tuple(b)), EnaRunLibrary) for b in batches(accessions, 500)]
    return Response(first_status(responses), [r for res in responses for r in res.data])


def sample_hosts(accessions: Iterable[str]) -> Response[list[EnaSample]]:
    """IDs and host of each sample."""
    return samples(accessions, EnaSample)


def sample_metadata(accessions: Iterable[str]) -> Response[list[EnaSampleMetadata]]:
    """The descriptive fields of each sample (title, description, env_*, isolation source, …)."""
    return samples(accessions, EnaSampleMetadata)


def sample_assemblies(samples: Iterable[str], clause: str = "") -> Response[list[EnaLinkedWgs | EnaLinkedAnalysis]]:
    """Every assembly or genome ENA links to these samples, of any type, taxon or submitter, as a
    sample's page in the ENA Browser lists them: assembly analyses and WGS sets. `clause` narrows
    them (e.g. ` AND NOT tax_tree(2759)`). 200 samples per query, by their primary IDs."""
    statuses: list[Response[Any]] = []
    found: list[EnaLinkedWgs | EnaLinkedAnalysis] = []
    for b in batches(samples, 200):
        ids = " OR ".join(f'sample_accession="{s}"' for s in b)
        analyses = search(
            EnaSearch("analysis", f"analysis_type=SEQUENCE_ASSEMBLY AND ({ids}){clause}"), EnaLinkedAnalysis
        )
        wgs = search(EnaSearch("wgs_set", f"({ids}){clause}"), EnaLinkedWgs)
        statuses += [analyses, wgs]
        found += [*analyses.data, *wgs.data]
    return Response(first_status(statuses), found)


def studies(accessions: Iterable[str]) -> Response[list[EnaStudy]]:
    """Title and description of each study, fetched by accession list."""
    responses = [search(EnaSearch("study", accessions=tuple(b)), EnaStudy) for b in batches(accessions, 300)]
    return Response(first_status(responses), [s for r in responses for s in r.data])


def emg_studies() -> Response[set[str]]:
    """ENA studies registered by MGnify: where it deposits the assemblies (TPA ERZ) and genomes
    it produces from public data. Center names vary ("EMG", "EMG;EMBL-EBI"), so any that
    contains EMG matches."""
    res = search(EnaSearch("study", 'center_name="*EMG*"', fields=("study_accession",)), EnaStudy)
    return Response(res.status, {s.study_accession for s in res.data})


# --- Browser API: sample XML ---


def derived_from(sample_accessions: Iterable[str]) -> Response[dict[str, list[str]]]:
    """{sample: [source sample IDs]} from each sample's `derived from` attribute (ENA sample
    XML, 100 samples per request). Keys are the primary IDs (SAMEA/SAMN/SAMD)."""
    out: dict[str, list[str]] = {}
    responses = [get(f"{ENA_BROWSER}/xml/{','.join(b)}") for b in batches(sample_accessions, 100)]
    for res in responses:
        for sample in ET.fromstring(res.data).iter("SAMPLE") if res.data else []:
            values = [
                a.findtext("VALUE") or ""
                for a in sample.iter("SAMPLE_ATTRIBUTE")
                if (a.findtext("TAG") or "").strip().lower() in DERIVED_FROM
            ]
            out[sample.get("accession", "")] = SAMPLE_ID.findall(" ".join(values))
    return Response(first_status(responses), out)


def _taxon(res: Response[str | None]) -> Response[EnaTaxon | None]:
    hits = json.loads(res.data) if res.data else []
    hits = hits if isinstance(hits, list) else [hits]
    return Response(res.status, from_row(EnaTaxon, hits[0]) if hits else None)


def taxon_by_name(scientific_name: str) -> Response[EnaTaxon | None]:
    """The taxon with this exact ENA scientific name, or None."""
    return _taxon(get(f"{ENA_TAXONOMY}/scientific-name/{urllib.parse.quote(scientific_name)}"))


def taxon(tax_id: int | str) -> Response[EnaTaxon | None]:
    """The taxon with this ENA/NCBI tax ID, or None."""
    return _taxon(get(f"{ENA_TAXONOMY}/tax-id/{tax_id}"))
