"""MGnify: API v2 (genome catalogues) and the genomes FTP. Catalogues come from v2 because
v1 lists outdated versions."""

import csv
import io
import json

from ..log import warn
from ..models import MgnifyCatalogue, MgnifyCatalogueGenome, Response, from_row
from .http import get

API_V2 = "https://www.ebi.ac.uk/metagenomics/api/v2"


# --- API v2: genome catalogues ---


def catalogue(catalogue_id: str) -> Response[MgnifyCatalogue | None]:
    """An MGnify genome catalogue, or None if it doesn't exist."""
    res = get(f"{API_V2}/genomes/catalogues/{catalogue_id}")
    return Response(res.status, from_row(MgnifyCatalogue, json.loads(res.data)) if res.data else None)


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
