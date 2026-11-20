"""Thin clients for the public APIs biome-finder reads (ENA, MGnify), with a disk cache."""

import csv
import hashlib
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api"
ENA_TAXONOMY = "https://www.ebi.ac.uk/ena/taxonomy/rest"
ENA_BROWSER = "https://www.ebi.ac.uk/ena/browser/api"
MGNIFY = "https://www.ebi.ac.uk/metagenomics/api/v1"  # biomes
MGNIFY_V2 = "https://www.ebi.ac.uk/metagenomics/api/v2"  # genome catalogues (v1 lists outdated versions)

CACHE_DIR = Path(os.environ.get("BIOME_FINDER_CACHE", "~/.cache/biome-finder")).expanduser()
CACHE_MAX_AGE = 7 * 24 * 3600  # ponytail: fixed one-week TTL, add a --refresh flag if stale data bites
USER_AGENT = "biome-finder/0.1"

csv.field_size_limit(10**9)


def get(url, params=None, retries=5, check=None, data=None):
    """GET a URL (text), cached on disk; retries with backoff on 429/5xx/network errors.

    `check(text) -> bool` rejects a bad response (e.g. truncated): it is retried and never cached.
    `data` (dict) makes it a form-encoded POST.
    """
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    body = urllib.parse.urlencode(data).encode() if data else None
    cache = CACHE_DIR / hashlib.sha1(url.encode() + (body or b"")).hexdigest()
    if cache.exists() and time.time() - cache.stat().st_mtime < CACHE_MAX_AGE:
        return cache.read_text()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as r:
                text = r.read().decode()
            if check is None or check(text):
                break
            if attempt == retries - 1:
                raise RuntimeError(f"bad response after {retries} attempts: {url}")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                raise
        except urllib.error.URLError:
            if attempt == retries - 1:
                raise
        time.sleep(2**attempt)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(text)
    return text


def ena_count(result, query):
    text = get(f"{ENA_PORTAL}/count", {"result": result, "query": query})
    return int(text.split()[-1])


def parse_tsv(text, fields):
    """Rows of an ENA TSV page; None if the page is truncated (ENA ends it with 'Premature EOF')."""
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if any(None in r.values() or None in r for r in rows) or "Premature EOF" in text[-200:]:
        return None
    return rows


def ena_search(result, query, fields):
    """All rows of an ENA Portal search, as dicts. Uses POST with limit=0 (large GET responses
    get truncated with 'Premature EOF'); the row total is checked against the count endpoint."""
    total = ena_count(result, query)
    if total == 0:
        return []
    form = {"result": result, "query": query, "fields": ",".join(fields), "format": "tsv", "limit": 0}
    text = get(f"{ENA_PORTAL}/search", data=form, check=lambda t: parse_tsv(t, fields) is not None)
    rows = parse_tsv(text, fields)
    if len(rows) != total:  # records published/suppressed between the two calls, or a silent cut
        warn(f"ENA {result}: {len(rows)} rows, count endpoint says {total} ({query})")
    return rows


SAMPLE_ID = re.compile(r"\b(SAM[END][A-Z]?\d+|[EDS]RS\d+)\b")
DERIVED_FROM = {"derived from", "derived_from", "sample derived from"}


def batches(items, size):
    items = sorted(set(items))
    return [items[i:i + size] for i in range(0, len(items), size)]


def derived_from(sample_accessions):
    """{sample: [source sample IDs]} from each sample's `derived from` attribute (ENA sample
    XML, 100 samples per request). Keys are the primary IDs (SAMEA/SAMN/SAMD)."""
    out = {}
    for batch in batches(sample_accessions, 100):
        text = get(f"{ENA_BROWSER}/xml/{','.join(batch)}")
        for sample in ET.fromstring(text).iter("SAMPLE") if text else []:
            values = [a.findtext("VALUE") or "" for a in sample.iter("SAMPLE_ATTRIBUTE")
                      if (a.findtext("TAG") or "").strip().lower() in DERIVED_FROM]
            out[sample.get("accession")] = SAMPLE_ID.findall(" ".join(values))
    return out


def sample_hosts(accessions):
    """{sample ID: ENA sample row with host fields}, keyed by both primary and secondary IDs.
    Primary IDs go through `includeAccessions`, which rejects secondary (ERS/SRS/DRS) ones."""
    fields = ["sample_accession", "secondary_sample_accession", "host", "host_tax_id"]
    primary = [a for a in accessions if a.startswith("SAM")]
    rows = []
    for batch in batches(primary, 500):
        form = {"result": "sample", "includeAccessions": ",".join(batch), "fields": ",".join(fields),
                "format": "tsv", "limit": 0}
        rows += parse_tsv(get(f"{ENA_PORTAL}/search", data=form,
                              check=lambda t: parse_tsv(t, fields) is not None), fields)
    for batch in batches(set(accessions) - set(primary), 100):
        rows += ena_search("sample", " OR ".join(f'secondary_sample_accession="{a}"' for a in batch), fields)
    return {acc: r for r in rows for acc in (r["sample_accession"], r["secondary_sample_accession"]) if acc}


def study_texts(accessions):
    """{study accession: title + description}, fetched by accession list (POST)."""
    fields = ["study_accession", "study_title", "study_description"]
    out = {}
    for batch in batches(accessions, 300):
        form = {"result": "study", "includeAccessions": ",".join(batch), "fields": ",".join(fields),
                "format": "tsv", "limit": 0}
        for r in parse_tsv(get(f"{ENA_PORTAL}/search", data=form,
                               check=lambda t: parse_tsv(t, fields) is not None), fields):
            out[r["study_accession"]] = f'{r["study_title"]} {r["study_description"]}'
    return out


def mgnify_catalogue_genomes(catalogue):
    """Every genome of an MGnify catalogue (its FTP `genomes-all_metadata.tsv`: Genome,
    Genome_type, Genome_accession, Sample_accession, …), or [] if the file isn't there."""
    url = catalogue["ftp_url"].replace("http://", "https://").rstrip("/") + "/genomes-all_metadata.tsv"
    text = get(url)
    if not text:
        warn(f"no genome list for {catalogue['catalogue_id']} at {url}")
        return []
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def mgnify_ena_studies():
    """ENA studies registered by MGnify (center EMG): where it deposits the assemblies (TPA ERZ)
    and MAGs it produces from public data."""
    # center names vary ("EMG", "EMG;EMBL-EBI"), so match any that contains EMG
    return {r["study_accession"] for r in ena_search("study", 'center_name="*EMG*"', ["study_accession"])}


def ena_taxid(scientific_name):
    """Tax ID for an exact ENA scientific name, or None."""
    text = get(f"{ENA_TAXONOMY}/scientific-name/{urllib.parse.quote(scientific_name)}")
    hits = json.loads(text) if text else []
    return hits[0]["taxId"] if hits else None


def ena_taxname(tax_id):
    """Scientific name for an ENA/NCBI tax ID, or the ID itself if unknown."""
    text = get(f"{ENA_TAXONOMY}/tax-id/{tax_id}")
    return json.loads(text)["scientificName"] if text else str(tax_id)


def mgnify_catalogue(catalogue_id):
    """An MGnify genome catalogue (catalogue_type, updated_at, genome_count, …), or None."""
    text = get(f"{MGNIFY_V2}/genomes/catalogues/{catalogue_id}")
    return json.loads(text) if text else None


def mgnify_subtree(lineage):
    """{lineage: samples-count} for a biome lineage and every descendant, or None if it doesn't
    exist. A biome's own `samples-count` covers only samples attached to that exact node."""
    url = f"{MGNIFY}/biomes/{urllib.parse.quote(lineage, safe=':')}/children?page_size=100"
    nodes = {}
    while url:  # the children listing includes the node itself and every descendant
        text = get(url)
        if not text:
            return None
        page = json.loads(text)
        nodes.update({d["id"]: d["attributes"]["samples-count"] for d in page["data"]})
        url = page["links"]["next"]
    return nodes


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)
