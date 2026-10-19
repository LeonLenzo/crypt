#!/usr/bin/env python3
"""_reads.py — what is on disk, and what counts as present.

Module-level support, beside `_assembly.py` and `_taxonomy.py`, because three steps need the
same answer to the same question: `03_cohort/migrate_reads.py` to know what to move,
`03_cohort/fetch_reads.py` to know what to skip, `04_classify` to know what it may classify.
That question was answered three separate ways in 03_kraken and the three disagreed.

**PRESENCE IS A NON-EMPTY FILE IN ITS FINAL LOCATION.** Never a directory, never a name in a
table. The 2026-09 CDS tree had 1,997 accession DIRECTORIES on scratch against 1,162 real
FASTAs, and a directory-existence check counted all 1,997 as held, so the build ran on a third
less sequence than it reported. The flat reads tree carries the same hazard in a second form:
39 orphaned `_<RUN>_sra` prefetch directories from downloads that died midway, which an
`iterdir()` loop will happily yield as runs.

THE LAYOUT, one directory per BioProject:

    reads/<BioProject>/<Run>_1.fastq.gz      paired
    reads/<BioProject>/<Run>_2.fastq.gz
    reads/<BioProject>/<Run>.fastq.gz        single-end, or a paired run whose second mate
                                             died in fasterq-dump -- which is recurring and
                                             not an edge case. 459 of the 3,224 runs already
                                             on scratch are single-file.
    reads_tmp/<Run>/                         in flight; NOT presence

A run is present if at least one of its final-location files clears MIN_FASTQ_BYTES. It is
`complete` only if the mate count matches what the archive says, which this module cannot know
on its own, so `n_files` is reported and the judgment is left to the caller.
"""

import csv
import re
import time
from collections import defaultdict
from pathlib import Path

# A gzipped FASTQ smaller than this is a truncation or an empty archive, not a library. The
# smallest real file among the 5,989 on scratch is comfortably above it; the floor exists to
# catch the 224-file silent truncation that a flat kraken2 timeout caused in 2026-09.
MIN_FASTQ_BYTES = 100_000

RUN_RE = re.compile(r"^([A-Z]RR\d+)")
# `_<RUN>_sra` is prefetch's own working directory. Matched explicitly so it can be reported
# as debris rather than silently skipped: an orphan means a download died, and that is a fact
# worth surfacing, not tidying away.
ORPHAN_RE = re.compile(r"^_([A-Z]RR\d+)_sra$")

MANIFEST_COLS = ["Run", "BioProject", "n_files", "bytes", "files", "atime_days"]


def run_of(name: str) -> str | None:
    m = RUN_RE.match(name)
    return m.group(1) if m else None


def fastq_files(d: Path) -> list[Path]:
    """Non-empty .fastq.gz directly in `d`. Not recursive, by design: a nested file is debris."""
    if not d.is_dir():
        return []
    return sorted(p for p in d.iterdir()
                  if p.is_file() and p.name.endswith(".fastq.gz")
                  and p.stat().st_size >= MIN_FASTQ_BYTES)


def scan_flat(d: Path) -> tuple[dict, list]:
    """The OLD layout: one directory, every run. Returns (run -> [paths], orphan dirs)."""
    out, orphans = defaultdict(list), []
    if not d.is_dir():
        return {}, []
    for p in sorted(d.iterdir()):
        if p.is_dir():
            if ORPHAN_RE.match(p.name):
                orphans.append(p)
            continue
        if not p.name.endswith(".fastq.gz") or p.stat().st_size < MIN_FASTQ_BYTES:
            continue
        run = run_of(p.name)
        if run:
            out[run].append(p)
    return dict(out), orphans


def scan_tree(reads: Path) -> dict:
    """The NEW layout. Returns run -> {BioProject, files, bytes, atime_days}.

    Keyed on run and not on (project, run), because a run belongs to exactly one BioProject and
    finding the same run under two would mean the plan and the disk disagree. The caller is
    given the first and told about the second rather than having one silently win.
    """
    out, seen = {}, {}
    if not reads.is_dir():
        return out
    now = time.time()
    for proj in sorted(p for p in reads.iterdir() if p.is_dir()):
        for f in fastq_files(proj):
            run = run_of(f.name)
            if not run:
                continue
            if run in seen and seen[run] != proj.name:
                raise RuntimeError(
                    f"{run} is present under BOTH {seen[run]} and {proj.name}. The tree and "
                    f"the fetch plan disagree about which project it belongs to; resolve it "
                    f"before classifying, or the run is counted twice.")
            seen[run] = proj.name
            e = out.setdefault(run, {"BioProject": proj.name, "files": [], "bytes": 0,
                                     "atime_days": 0.0})
            st = f.stat()
            e["files"].append(f)
            e["bytes"] += st.st_size
            e["atime_days"] = max(e["atime_days"], (now - st.st_atime) / 86400)
    return out


def write_manifest(tree: dict, path: Path) -> None:
    """One row per run present. This is the FACT table; run_list.tsv is the plan."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MANIFEST_COLS)
        for run in sorted(tree):
            e = tree[run]
            w.writerow([run, e["BioProject"], len(e["files"]), e["bytes"],
                        ";".join(f.name for f in sorted(e["files"])),
                        f"{e['atime_days']:.1f}"])
    tmp.replace(path)
