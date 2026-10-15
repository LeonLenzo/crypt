#!/usr/bin/env python3
"""_seeds.py — read the seeds/ tables, and refuse rather than return nothing.

Written 2026-10-15 because the first loader for `exclusions.tsv` reported **0 taxids blocked**
and looked like it had worked. Four of the six seeds tables open with a `#` comment preamble
(9 lines in `exclusions.tsv`, 15 in `hosts.tsv` and `residents.tsv`, 8 in `sources.tsv`), so a
bare `csv.DictReader` takes the first comment line as the header. Every row then has one key
whose name is a sentence, `row.get("ncbi_taxid")` is None for all of them, and the loop selects
nothing.

That failure is worse than a crash in the one way that matters here: **it fails permissive.**
All 44 exclusion decisions would have been ignored silently, including the fungi-only oomycete
exclusion and the *Capsicum annuum* catch (a PLANT that USDA's fungus column legitimately
records as recorded-on, which was being selected as a reference with 6 annotated assemblies).
`seeds/README.md` says nothing may be dropped from the build without a row in
`exclusions.tsv`; the inverse has to hold too, and a row there has to actually bite.

So `read_seed_tsv` does two things a DictReader does not:

  * skips the `#` preamble before looking for the header, and
  * takes the columns the caller needs and REFUSES if any is missing, naming what it found.

`residents.tsv` has no header at all, which is why `read_seed_rows` exists alongside.

The same discipline as `01_literature/_resolver.py`, which writes a join table and then proves
it joined: a loader that can return empty for a structural reason must say so, not shrug.
"""

import csv
from pathlib import Path


def _data_lines(path: Path):
    """Yield the lines of a seeds table with the `#` preamble and blank lines removed.

    Only a comment at the START of a line is a comment. A `#` inside a value is data: the
    reason columns carry prose and `exclusions.tsv` cites things like "CLAUDE.md active task 5".
    """
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield line


def read_seed_tsv(path: Path, require=()) -> list:
    """Rows of a headered seeds table as dicts, after the preamble.

    `require` names the columns the caller actually reads. Any missing one raises, because the
    alternative is a loop that quietly matches nothing. Pass it: the whole point of this
    module is that an empty result has to be a fact about the data and not about the parse.
    """
    rows = list(csv.DictReader(_data_lines(path), delimiter="\t"))
    if not rows:
        raise ValueError(f"{path} has no data rows below its preamble")
    missing = [c for c in require if c not in rows[0]]
    if missing:
        raise ValueError(
            f"{path}: missing column(s) {missing}. Found {list(rows[0])}. "
            "If the first key looks like a sentence, the `#` preamble was read as the header")
    return rows


def read_seed_rows(path: Path, min_cols: int = 1) -> list:
    """Rows of a HEADERLESS seeds table as lists of fields, after the preamble.

    `residents.tsv` is the case: its 15-line preamble documents the columns and the data
    starts immediately, so there is no header row to name them. Rows with fewer than
    `min_cols` fields are returned as-is; the caller decides what a short row means.
    """
    rows = [line.rstrip("\n").split("\t") for line in _data_lines(path)]
    if not rows:
        raise ValueError(f"{path} has no data rows below its preamble")
    short = sum(1 for r in rows if len(r) < min_cols)
    if short == len(rows):
        raise ValueError(
            f"{path}: every row has fewer than {min_cols} tab-separated fields. "
            "Check the file is tab-separated and the preamble is `#`-prefixed")
    return rows
