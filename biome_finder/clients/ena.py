"""ENA: Portal API (searches), Browser API (sample XML) and Taxonomy REST."""

import csv
import io
import json
import re
import urllib.parse
import xml.etree.ElementTree as ET
from collections.abc import Iterable

from ..log import warn
from ..models import M, EnaSample, EnaSearch, EnaStudy, EnaTaxon, Response, columns, from_row
from .http import batches, first_status, get

PORTAL = "https://www.ebi.ac.uk/ena/portal/api"
TAXONOMY = "https://www.ebi.ac.uk/ena/taxonomy/rest"
BROWSER = "https://www.ebi.ac.uk/ena/browser/api"

csv.field_size_limit(10**9)


# --- Portal API ---


def count(request: EnaSearch) -> int:
    text = get(f"{PORTAL}/count", {"result": request.result, "query": request.query}).data
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
    res = get(f"{PORTAL}/search", data=form, check=lambda t: parse_tsv(t, fields) is not None)
    rows = (parse_tsv(res.data, fields) if res.data else None) or []
    if total is not None and len(rows) != total:
        warn(f"ENA {request.result}: {len(rows)} rows, count endpoint says {total} ({request.query})")
    return Response(res.status, [from_row(model, r) for r in rows])


def sample_hosts(accessions: Iterable[str]) -> Response[list[EnaSample]]:
    """IDs and host of each sample. Primary IDs go through `includeAccessions`, which rejects
    secondary (ERS/SRS/DRS) ones; those go through an OR query."""
    accessions = set(accessions)
    primary = {a for a in accessions if a.startswith("SAM")}
    responses = [search(EnaSearch("sample", accessions=tuple(b)), EnaSample) for b in batches(primary, 500)]
    responses += [
        search(EnaSearch("sample", " OR ".join(f'secondary_sample_accession="{a}"' for a in b)), EnaSample)
        for b in batches(accessions - primary, 100)
    ]
    return Response(first_status(responses), [s for r in responses for s in r.data])


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

SAMPLE_ID = re.compile(r"\b(SAM[END][A-Z]?\d+|[EDS]RS\d+)\b")
DERIVED_FROM = {"derived from", "derived_from", "sample derived from"}


def derived_from(sample_accessions: Iterable[str]) -> Response[dict[str, list[str]]]:
    """{sample: [source sample IDs]} from each sample's `derived from` attribute (ENA sample
    XML, 100 samples per request). Keys are the primary IDs (SAMEA/SAMN/SAMD)."""
    out: dict[str, list[str]] = {}
    responses = [get(f"{BROWSER}/xml/{','.join(b)}") for b in batches(sample_accessions, 100)]
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
    return _taxon(get(f"{TAXONOMY}/scientific-name/{urllib.parse.quote(scientific_name)}"))


def taxon(tax_id: int | str) -> Response[EnaTaxon | None]:
    """The taxon with this ENA/NCBI tax ID, or None."""
    return _taxon(get(f"{TAXONOMY}/tax-id/{tax_id}"))
