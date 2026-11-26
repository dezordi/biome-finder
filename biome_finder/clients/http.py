"""HTTP GET/POST with a disk cache and retries, shared by every client."""

import hashlib
import os
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable
from pathlib import Path

from ..log import debug, info
from ..models import Response

CACHE_DIR = Path(os.environ.get("BIOME_FINDER_CACHE", "~/.cache/biome-finder")).expanduser()
CACHE_MAX_AGE = 7 * 24 * 3600  # ponytail: fixed one-week TTL, add a --refresh flag if stale data bites
USER_AGENT = "biome-finder/0.1"


def get(
    url: str,
    params: dict[str, str] | None = None,
    retries: int = 5,
    check: Callable[[str], bool] | None = None,
    data: dict[str, object] | None = None,
) -> Response[str | None]:
    """GET a URL (text), cached on disk; retries with backoff on 429/5xx/network errors.
    404 gives `Response(404, None)`; other HTTP errors raise.

    `check(text) -> bool` rejects a bad response (e.g. truncated): it is retried and never cached.
    `data` makes it a form-encoded POST.
    """
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    body = urllib.parse.urlencode(data).encode() if data else None
    cache = CACHE_DIR / hashlib.sha1(url.encode() + (body or b"")).hexdigest()
    if cache.exists() and time.time() - cache.stat().st_mtime < CACHE_MAX_AGE:
        debug(f"cached {url}")
        return Response(200, cache.read_text())
    debug(f"{'POST' if body else 'GET'} {url}")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as r:
                status, text = r.status, r.read().decode()
            if check is None or check(text):
                break
            if attempt == retries - 1:
                raise RuntimeError(f"bad response after {retries} attempts: {url}")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return Response(404, None)
            if e.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                raise
            info(f"HTTP {e.code}, retrying in {2**attempt} s: {url}")
        except urllib.error.URLError as e:
            if attempt == retries - 1:
                raise
            info(f"{e.reason}, retrying in {2**attempt} s: {url}")
        time.sleep(2**attempt)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(text)
    return Response(status, text)


def download(url: str, retries: int = 5) -> Response[Path | None]:
    """A file saved once under the cache's `files/` folder, by its name, for files that never
    change once published (e.g. a monthly taxdump release). 404 gives `Response(404, None)`."""
    path = CACHE_DIR / "files" / url.rsplit("/", 1)[-1]
    if path.exists():
        return Response(200, path)
    info(f"downloading {url} to {path} (once; later runs reuse it)")
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + ".part")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as r, open(part, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            part.rename(path)
            return Response(200, path)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return Response(404, None)
            if e.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == retries - 1:
                raise
        time.sleep(2**attempt)
    raise RuntimeError(f"download failed: {url}")


def first_status(responses: Iterable[Response]) -> int:
    """The status of a response built from several requests: the first that isn't 200."""
    return next((r.status for r in responses if r.status != 200), 200)


def batches(items: Iterable[str], size: int) -> list[list[str]]:
    """Unique items, sorted, in chunks of `size` (one request each)."""
    items = sorted(set(items))
    return [items[i : i + size] for i in range(0, len(items), size)]
