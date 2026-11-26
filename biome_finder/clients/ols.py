"""OLS4: ontology terms and their subclasses (for ENVO)."""

import json
import urllib.parse
from typing import Any

from ..models import OlsTerm, Response
from .http import first_status, get
from ..constants import OLS_API


def term_url(obo_id: str) -> str:
    """A term's OLS_API URL; OLS_API wants its IRI URL-encoded twice."""
    iri = f"http://purl.obolibrary.org/obo/{obo_id.replace(':', '_')}"
    ontology = obo_id.split(":")[0].lower()
    return f"{OLS_API}/{ontology}/terms/{urllib.parse.quote(urllib.parse.quote(iri, safe=''), safe='')}"


def parse(t: dict[str, Any]) -> OlsTerm:
    return OlsTerm(t["obo_id"], t["label"], tuple(t.get("synonyms") or ()))


def with_descendants(obo_id: str) -> Response[list[OlsTerm]]:
    """A term and everything under it (`hierarchicalDescendants`: subclasses and parts, the way
    ENVO nests e.g. 'freshwater lake biome' under 'freshwater lake'). Empty if OLS_API doesn't know
    the term."""
    url = term_url(obo_id)
    res = get(url)
    if not res.data:
        return Response(res.status, [])
    terms, responses, page = [parse(json.loads(res.data))], [res], 0
    while True:
        r = get(f"{url}/hierarchicalDescendants", {"size": "500", "page": str(page)})
        responses.append(r)
        data = json.loads(r.data) if r.data else {}
        terms += [parse(t) for t in data.get("_embedded", {}).get("terms", [])]
        page += 1
        if page >= data.get("page", {}).get("totalPages", 0):
            break
    return Response(first_status(responses), terms)
