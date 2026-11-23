"""MGnify: API v2 (genome catalogues), API v1 (biomes) and the genomes FTP.

Catalogues come from v2 because v1 lists outdated versions. Biomes stay on v1: v2 lists the
tree without sample counts, and its samples endpoint holds only part of MGnify's samples."""

import csv
import io
import json
import urllib.parse

from ..log import warn
from ..models import MgnifyBiome, MgnifyCatalogue, MgnifyCatalogueGenome, Response, from_row
from .http import get

API_V1 = "https://www.ebi.ac.uk/metagenomics/api/v1"
API_V2 = "https://www.ebi.ac.uk/metagenomics/api/v2"


# --- API v2: genome catalogues ---


def catalogue(catalogue_id: str) -> Response[MgnifyCatalogue | None]:
    """An MGnify genome catalogue, or None if it doesn't exist."""
    res = get(f"{API_V2}/genomes/catalogues/{catalogue_id}")
    return Response(res.status, from_row(MgnifyCatalogue, json.loads(res.data)) if res.data else None)


# --- API v1: biomes ---


def biome_subtree(lineage: str) -> Response[list[MgnifyBiome]]:
    """A biome lineage and every descendant, from its paged children listing (which includes
    the node itself), each with its own samples count: a node's count covers only samples
    attached to that exact node. Status 404 if the lineage doesn't exist."""
    url: str | None = f"{API_V1}/biomes/{urllib.parse.quote(lineage, safe=':')}/children?page_size=100"
    nodes: list[MgnifyBiome] = []
    while url:
        res = get(url)
        if not res.data:
            return Response(res.status, [])
        page = json.loads(res.data)
        nodes += [MgnifyBiome(d["id"], d["attributes"]["samples-count"]) for d in page["data"]]
        url = page["links"]["next"]
    return Response(200, nodes)


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
