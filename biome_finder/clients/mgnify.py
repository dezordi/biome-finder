"""MGnify: API v2 (genome catalogues, studies) and the genomes FTP. Catalogues come from v2 because
v1 lists outdated versions."""

import csv
import io
import json

from ..log import warn
from ..models import MgnifyCatalogue, MgnifyCatalogueGenome, MgnifyStudy, Response, from_row
from .http import first_status, get
from ..constants import MGNIFY_API_V2

# --- API v2: genome catalogues ---


def catalogue(catalogue_id: str) -> Response[MgnifyCatalogue | None]:
    """An MGnify genome catalogue, or None if it doesn't exist."""
    res = get(f"{MGNIFY_API_V2}/genomes/catalogues/{catalogue_id}")
    return Response(res.status, from_row(MgnifyCatalogue, json.loads(res.data)) if res.data else None)


# --- API v2: studies ---


def studies() -> Response[list[MgnifyStudy]]:
    """Every study MGnify has analysed, with its curated biome (about 5,800, 100 per page)."""
    found, responses, page = [], [], 1
    while True:
        res = get(f"{MGNIFY_API_V2}/studies/", {"page": str(page), "page_size": "100"})
        responses.append(res)
        data = json.loads(res.data) if res.data else {}
        for s in data.get("items", []):
            found.append(
                MgnifyStudy(
                    s["accession"], tuple(s.get("ena_accessions") or ()), (s.get("biome") or {}).get("lineage") or ""
                )
            )
        if not data.get("items") or len(found) >= data.get("count", 0):
            break
        page += 1
    return Response(first_status(responses), found)


# --- FTP: catalogue genome lists ---


def catalogue_genomes(cat: MgnifyCatalogue) -> Response[list[MgnifyCatalogueGenome]]:
    """Every genome of an MGnify catalogue, from its FTP `genomes-all_metadata.tsv`; status 404
    and no genomes if the file isn't there."""
    url = cat.ftp_url.replace("http://", "https://").rstrip("/") + "/genomes-all_metadata.tsv"
    res = get(url)
    if not res.data:
        warn(f"no genome list for {cat.catalogue_id} at {url}")
        return Response(res.status, [])
    rows = csv.DictReader(io.StringIO(res.data), delimiter="\t")
    return Response(res.status, [from_row(MgnifyCatalogueGenome, r) for r in rows])
