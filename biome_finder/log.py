"""Warnings on stderr."""

import sys


def warn(msg: str) -> None:
    print(f"warning: {msg}", file=sys.stderr)
