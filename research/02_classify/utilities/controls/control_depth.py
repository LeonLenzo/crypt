#!/usr/bin/env python3
"""control_depth.py -- per-run depth for the controls, in run_depth.tsv's schema.

    python 02_classify/utilities/controls/control_depth.py

WHY THIS EXISTS. `05_detect/detect.py` keys every run against `silver/run_depth.tsv` and skips
anything absent from it. A control is not in `run_list.tsv` and never will be -- that is the
whole point of keeping it outside the cohort tree -- so every control run was dropped at that
join, silently. Recorded as a gap on 2026-10-26 and closed here rather than with another
throwaway script that re-implements detect.py's arithmetic.

WHAT detect.py ACTUALLY NEEDS is two fields: `read_len` and `layout`, which it turns into
k-mers per spot. Both are already in the ENA filereport captures that `fetch_controls.py`
writes to `bronze/controls/meta/`, so nothing new is fetched:

    read_len = base_count / read_count     bases PER SPOT, both mates summed for PAIRED
    layout   = PAIRED if the fastq_ftp list has two or more URLs, else SINGLE

The remaining columns are filled so the file is a drop-in for the same reader, not because
anything downstream reads them.

ENA's `base_count` is the spot total, which is why `kmers_per_spot()` halves it for PAIRED
rather than doubling. Getting that backwards doubles redundancy and halves saturation.
"""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import CONTROL_META, CONTROL_REPORTS, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

OUT = SILVER / "control_depth.tsv"
COLS = ["Run", "BioProject", "layout", "strategy", "selection", "reads", "bases",
        "read_len", "platform", "model", "n_fastq", "fastq_bytes", "fastq_urls"]


def main() -> None:
    # Only accessions that have actually been classified: a capture exists for every
    # CANDIDATE study, including the six not yet fetched, and a depth row for a run with no
    # report is harmless but misleading when the file is read as "what we have".
    have = {d.name for d in CONTROL_REPORTS.iterdir() if d.is_dir()} \
        if CONTROL_REPORTS.exists() else set()
    classified = {r.name.split(".")[0]
                  for d in CONTROL_REPORTS.iterdir() if d.is_dir()
                  for r in d.glob("*.report.tsv")} if CONTROL_REPORTS.exists() else set()

    rows, skipped = [], []
    for cap in sorted(CONTROL_META.glob("*.tsv")):
        for r in csv.DictReader(open(cap), delimiter="\t"):
            run = r["run_accession"]
            if run not in classified:
                continue
            urls = [u for u in (r.get("fastq_ftp") or "").split(";") if u]
            try:
                reads, bases = int(r["read_count"]), int(r["base_count"])
            except (ValueError, KeyError, TypeError):
                skipped.append((run, "no read_count/base_count")); continue
            if not reads:
                skipped.append((run, "read_count = 0")); continue
            rows.append({
                "Run": run, "BioProject": r.get("study_accession") or cap.stem,
                "layout": "PAIRED" if len(urls) >= 2 else "SINGLE",
                "strategy": r.get("library_strategy", ""),
                "selection": r.get("library_selection", ""),
                "reads": reads, "bases": bases,
                "read_len": round(bases / reads, 1),
                "platform": r.get("instrument_platform", ""),
                "model": r.get("instrument_model", ""),
                "n_fastq": len(urls),
                "fastq_bytes": r.get("fastq_bytes", ""),
                "fastq_urls": r.get("fastq_ftp", ""),
            })

    with atomic_write(OUT) as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["BioProject"], r["Run"])))

    print(f"-> {OUT}: {len(rows)} runs over {len({r['BioProject'] for r in rows})} accessions")
    by = {}
    for r in rows:
        by.setdefault(r["BioProject"], []).append(r)
    for acc, rs in sorted(by.items()):
        gk = [r["bases"] / 1e9 for r in rs]
        print(f"   {acc:14s} {len(rs):3d} runs  {rs[0]['layout']:6s} "
              f"{rs[0]['read_len']:6.1f} bp/spot  {sum(gk):7.2f} Gbase total")
    # A classified run with no capture row is a real inconsistency, not a tidy-up.
    orphan = classified - {r["Run"] for r in rows}
    if orphan:
        print(f"   WARNING: {len(orphan)} classified run(s) have no ENA capture row: "
              f"{', '.join(sorted(orphan)[:6])}", file=sys.stderr)
    for run, why in skipped:
        print(f"   skipped {run}: {why}", file=sys.stderr)
    print(f"   accessions with reports: {len(have)}")


if __name__ == "__main__":
    main()
