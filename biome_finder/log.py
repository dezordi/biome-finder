"""Progress and warnings on stderr, through the `biome_finder` logger (set up by the CLI)."""

import logging
import sys

logger = logging.getLogger("biome_finder")


def setup(verbose: bool = False) -> None:
    """Timestamped lines on stderr: INFO by default, DEBUG (every HTTP request) with `verbose`."""
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S"))
    logger.handlers = [handler]
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)


def info(msg: str) -> None:
    logger.info(msg)


def debug(msg: str) -> None:
    logger.debug(msg)


def warn(msg: str) -> None:
    logger.warning(msg)


def counts(items: dict[str, int]) -> str:
    """`a 3, b 1`, largest first, thousands separated, zeros left out."""
    return ", ".join(f"{k} {v:,}" for k, v in sorted(items.items(), key=lambda kv: -kv[1]) if v) or "none"
