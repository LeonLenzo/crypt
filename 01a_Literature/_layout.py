#!/usr/bin/env python3
"""Where each kind of file lives in 01a_Literature, and which kind may be regenerated.

Named `_layout` and not `_paths` because the repo root already has a `_paths.py`, and a
second module of that name shadows it for whichever script imports first. That bit
immediately: review.py resolved `_paths` to the root module and lost SEEDS.

Added 2026-10-06 after three separate incidents in one day, all the same root cause: nothing
in the layout distinguished a file that may be rebuilt from a file that must never be
overwritten. `papers.tsv` was hand-annotated and the next `triage.py` run silently discarded
the annotations; an imported accession list survived only as a shell argument; and the file
created to fix that was gitignored, because tracking depended on whether someone remembered
`git add -f` rather than on what kind of file it was.

    SEEDS   HAND-WRITTEN and irreplaceable. The module's whole intellectual asset: every
            curated value traces to a rule, a join and a quote in here. Tracked in git by
            default. NEVER written by a script except by appending a reviewed row.
    BRONZE  immutable captures of somebody else's data, under DATA for now: runinfo/,
            biosample_attrs/, geo/, ena/. Deleting one costs a refetch and nothing else.
    DATA    generated tables (silver): runs.tsv, bioprojects.tsv, papers.tsv,
            paper_bioproject.tsv, registry.tsv, worklist.tsv. Rebuildable from BRONZE +
            SEEDS, so gitignored. Hand-editing one is always a mistake.
    GOLD    small human-readable outputs: host_summary.tsv, offtarget.tsv, the gap analyses.
            Generated, but tracked anyway because they are summary numbers a reader wants to
            see change over the life of the project, and they are a few kilobytes each.

The rule in one line: **if losing a file costs refetching, it is BRONZE or DATA; if it costs
re-reading papers, it is SEEDS.**
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SEEDS = HERE / "seeds"
GOLD = HERE / "gold"
STUDIES = HERE / "studies"

# Bronze lives under DATA until stage two of the layout split.
RUNINFO = DATA / "runinfo"
BIOSAMPLE = DATA / "biosample_attrs"
GEO = DATA / "geo"
ENA = DATA / "ena"
