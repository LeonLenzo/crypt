"""Shared caption plumbing for the figure directories.

Each figure owns its own caption.py, beside the panel it describes. What they share is only
the mechanics: read a TSV, wrap a paragraph. The words live with the figure.

Wrapping happens at WRITE TIME, not in the source, so a number changing width cannot leave a
ragged line behind. break_on_hyphens is off because it splits "no-op" and "PHI-base" across
lines. See crypt-figure-caption-discipline: every number in a caption is a claim about the
data, and a claim nobody can recompute is one nobody can check.
"""

import csv
import textwrap
from pathlib import Path

WIDTH = 78
csv.field_size_limit(10 ** 7)


def rd(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def para(*chunks):
    """Wrap each paragraph. "" is a blank line; a chunk with a newline passes through."""
    out = []
    for c in chunks:
        if c == "":
            out.append("")
        elif "\n" in c:
            out.append(c)
        else:
            out.append(textwrap.fill(c, WIDTH, break_on_hyphens=False))
    return "\n".join(out)


def write(path, text):
    Path(path).write_text(text.rstrip("\n") + "\n", encoding="utf-8")
    print(f"  -> {path}")
