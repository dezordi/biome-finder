"""JGI GOLD: the public data export ("Public Studies/Biosamples/SPs/APs/Organisms Excel").

GOLD's download page is behind a browser check, so the curator downloads the file (see
docs/search.md). Reading it takes minutes, so the BioSample → ecosystem table is cached as a
TSV next to the HTTP cache, keyed by the file's path, size and modification time."""

import csv
import hashlib
import warnings
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from ..log import info
from ..models import GoldBiosample, Response
from .http import CACHE_DIR

LEVELS = ["ECOSYSTEM", "ECOSYSTEM CATEGORY", "ECOSYSTEM TYPE", "ECOSYSTEM SUBTYPE", "SPECIFIC ECOSYSTEM"]


def sheet(workbook: Any, name: str, columns: list[str]) -> Iterator[list[Any]]:
    """The `columns` of a sheet, row by row. The export doesn't declare its sheet sizes, which
    read-only mode needs, hence reset_dimensions; a row ending in empty cells comes back short."""
    ws = workbook[name]
    ws.reset_dimensions()
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    idx = [header.index(c) for c in columns]
    for r in rows:
        yield [r[i] if i < len(r) else None for i in idx]


def read_export(path: Path) -> list[GoldBiosample]:
    """Every BioSample accession of a GOLD sequencing project with its ecosystem path."""
    import openpyxl  # only needed here

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # "Workbook contains no default style"
        wb = openpyxl.load_workbook(path, read_only=True)
        paths = {}
        for kind in ("BIOSAMPLE", "ORGANISM"):
            for gid, *levels in sheet(wb, kind.title(), [f"{kind} GOLD ID", *(f"{kind} {lv}" for lv in LEVELS)]):
                paths[gid] = " > ".join(str(v or "") for v in levels)
        found = {}
        for acc, bs, org in sheet(
            wb, "Sequencing Project", ["NCBI BIOSAMPLE ACCESSION", "BIOSAMPLE GOLD ID", "ORGANISM GOLD ID"]
        ):
            gid = bs if bs in paths else org
            if acc and gid in paths:
                found[acc] = GoldBiosample(acc, gid, paths[gid])
    return list(found.values())


def biosamples(path: str | Path) -> Response[list[GoldBiosample]]:
    """The export's BioSample → ecosystem table, from the cache when the file hasn't changed."""
    path = Path(path).expanduser().resolve()
    stat = path.stat()
    key = hashlib.sha1(f"{path}:{stat.st_size}:{stat.st_mtime_ns}".encode()).hexdigest()
    cache = CACHE_DIR / "files" / f"gold-{key}.tsv"
    fields = ["biosample_accession", "gold_id", "ecosystem_path"]
    if cache.exists():
        with open(cache, newline="") as f:
            return Response(200, [GoldBiosample(**r) for r in csv.DictReader(f, delimiter="\t")])
    info(f"reading the GOLD export {path} (a few minutes, once; cached as {cache})")
    found = read_export(path)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with open(cache, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows({k: getattr(g, k) for k in fields} for g in found)
    return Response(200, found)
