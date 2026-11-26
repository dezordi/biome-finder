"""NCBI: the monthly taxdump archive (later E-utilities and the genomes FTP)."""

import io
import re
import zipfile
from collections.abc import Iterable, Iterator

from ..models import NCBI_RANKS, NcbiTaxdump, NcbiTaxon, Response
from .http import download, get

TAXDUMP_ARCHIVE = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive"


def latest_taxdump() -> Response[str]:
    """The newest release in the archive (YYYY-MM-01), from its directory listing."""
    res = get(f"{TAXDUMP_ARCHIVE}/")
    return Response(res.status, max(re.findall(r"new_taxdump_(\d{4}-\d{2}-01)\.zip", res.data or ""), default=""))


def dmp(z: zipfile.ZipFile, name: str) -> Iterator[list[str]]:
    """The rows of one .dmp file, split on its `\\t|\\t` separator."""
    with z.open(name) as f:
        for line in io.TextIOWrapper(f, encoding="utf-8"):
            yield line.rstrip("\n").removesuffix("\t|").split("\t|\t")


def taxdump(release: str, tax_ids: Iterable[str]) -> Response[NcbiTaxdump | None]:
    """What `tax_ids` need from the release `new_taxdump_<release>.zip` (YYYY-MM-01): merges,
    deletions and the current taxa with their ranks. Only those rows are kept, so the full
    tree (2.7 M taxa) is never held in memory. None if the release isn't in the archive."""
    res = download(f"{TAXDUMP_ARCHIVE}/new_taxdump_{release}.zip")
    if not res.data:
        return Response(res.status, None)
    wanted = {t for t in tax_ids if t}
    with zipfile.ZipFile(res.data) as z:
        merged = {old: new for old, new in dmp(z, "merged.dmp") if old in wanted}
        deleted = frozenset(t for (t,) in dmp(z, "delnodes.dmp") if t in wanted)
        current = {merged.get(t, t) for t in wanted}
        ancestors = {t: frozenset(a.split()) for t, a in dmp(z, "taxidlineage.dmp") if t in current}
        own_rank = {r[0]: r[2] for r in dmp(z, "nodes.dmp") if r[0] in current}
        taxa = {}
        for t, name, *names in dmp(z, "rankedlineage.dmp"):  # species, genus, …, domain: above the taxon only
            if t in current:
                ranks = dict(zip(NCBI_RANKS, reversed([n or None for n in names])))
                if own_rank.get(t) in ranks:
                    ranks[own_rank[t]] = name
                ranks["class_"] = ranks.pop("class")
                taxa[t] = NcbiTaxon(t, name, own_rank.get(t, ""), **ranks, ancestors=ancestors.get(t, frozenset()))
    return Response(200, NcbiTaxdump(release, merged, deleted, taxa))
