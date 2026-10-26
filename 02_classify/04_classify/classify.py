#!/usr/bin/env python3
"""classify.py — run db_v4 over the cohort reads on disk, one BioProject at a time.

    python 02_classify/04_classify/classify.py --list
    python 02_classify/04_classify/classify.py --project PRJEB39201 --db ... --threads 16
    python 02_classify/04_classify/classify.py --project PRJEB39201 --force

Reads come from `reads/<BioProject>/`, reports go to `reports/<BioProject>/<Run>.report.tsv`.
One report per run, and a run with a report is DONE -- which is the check that makes "is this
project finished" a directory listing, and the check that licenses deleting its reads.

## No depth filter, and no selection of any kind

Every run on disk is classified. That is the clean-run rule, and after the 2026-10-19
rarefaction it is also the statistically correct choice: **a detection is valid at any depth;
only ABSENCE is depth-dependent.** A pathogen recovered from a 0.019 Gkmer library had to be
abundant to appear at all, so it is a stronger observation than the same call from a deep one,
not a weaker one. Nothing is gained by withholding shallow libraries from classification and a
great deal is lost -- a filter here would be a researcher degree of freedom applied before the
data is seen.

The depth caveat attaches to the DENOMINATOR of a prevalence statement, never to a call, and
that is 05_detect's business. `run_depth.tsv` carries Gkmer per run so 05_detect can report
each call beside the sensitivity that produced it. The measured relationship is in
`report/figures/rarefaction_caption.txt`: retention of established detections is 23% at 0.019
Gkmer, 55% at 0.21, 73% at 1.0.

## Resumable, and the resume check is a NON-EMPTY REPORT

A report that exists but is empty means kraken2 died partway, which happened to 443 of 3,225
runs in the 2026-09-25 pass when $MYSCRATCH pointed at the decommissioned fl3. Existence alone
would mark those done forever. Same rule as everywhere else in this module: presence is a
non-empty file.

## Why one project per job, and why that is not arbitrary

db_v4's hash is 41 GB and the memory request sets the CPU bill, so the expensive thing is
LOADING the database, not classifying against it. A per-run array task pays one 41 GB Lustre
load per run. One long-lived job per project pays `--workers` loads and amortises them over
every run in the project. For the big projects it is worth staging the hash into /dev/shm
first and passing --mmap, exactly as the rarefaction job does; for a 20-run project the stage
costs more than it saves.

**--memory-mapping against LUSTRE is catastrophic** (measured: 13 minutes, zero reports,
kraken2 at 1.6% CPU, because every hash lookup becomes a filesystem RPC). It is only correct
against /dev/shm. The flag is therefore not defaulted on and the SLURM script decides.
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import DB, READS, READS_MANIFEST, REPORTS, SILVER  # noqa: E402
from _reads import scan_tree, write_manifest  # noqa: E402
from _reports import parse_report  # noqa: E402

# db_v3's published operating point, carried across unchanged so a db_v4 number differs
# because the DATABASE changed and not because the settings did.
CONFIDENCE = "0.15"
MIN_HIT_GROUPS = "3"

# A kraken2 report smaller than this is a crash, not a sample: even a run with no fungal reads
# writes the unclassified/root lines.
MIN_REPORT_BYTES = 200

STATUS = SILVER / "classify_status.tsv"
STATUS_COLS = ["Run", "BioProject", "status", "seconds", "report_bytes", "species_rows",
               "reads_classified"]


def manifest(reads: Path = READS, manifest_path: Path = READS_MANIFEST) -> dict:
    """What is on disk, by SCANNING, never by reading a table.

    This read `reads_manifest.tsv` until 2026-10-21, and that was wrong in the specific way
    `_reads.py` opens by warning about: presence is a non-empty file in its final location,
    never a name in a table. The manifest is written by `migrate_reads.py`, which is a one-shot
    that ran once against the tree as it stood that day. Every run fetched afterwards was
    invisible here, so a classify job over six freshly downloaded projects reported
    "no reads on disk" and exited 0 against 120 runs and 305 GB that were sitting right there.

    A stale snapshot fails SILENTLY and in the safe-looking direction: it under-reports, so the
    job succeeds at doing nothing. Scanning costs one stat per file and cannot drift.

    The manifest is still WRITTEN, as a record of what the tree held at classify time, but
    nothing reads it to decide what exists.
    """
    tree = scan_tree(reads)
    if not tree:
        sys.exit(f"no reads found under {reads}\n"
                 "Fetch some: python 02_classify/03_cohort/fetch_reads.py --project <PRJ>")
    try:
        write_manifest(tree, manifest_path)
    except OSError as e:                        # a read-only checkout must not stop a classify
        print(f"  (could not refresh {manifest_path.name}: {e})", file=sys.stderr)
    # Return the scan entries UNCHANGED. The three keys the rest of this file reads
    # (BioProject, files, bytes) are exactly what scan_tree already provides, and `files` is a
    # list of absolute Paths rather than a ";"-joined string of basenames that has to be
    # rejoined to READS/<proj>. Rebuilding a narrower dict here is what broke the first attempt
    # at this fix: it dropped `files` and every worker raised KeyError.
    return tree


def done(report: Path) -> bool:
    """Presence is a NON-EMPTY report. See the module docstring."""
    return report.exists() and report.stat().st_size >= MIN_REPORT_BYTES


def classify_run(run: str, m: dict, db: Path, threads: int, mmap: bool,
                 tmp_dir: Path, reports: Path = REPORTS) -> dict:
    proj = m["BioProject"]
    out_dir = reports / proj
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / f"{run}.report.tsv"
    if done(report):
        return {"Run": run, "BioProject": proj, "status": "cached", "seconds": 0,
                "report_bytes": report.stat().st_size, "species_rows": "",
                "reads_classified": ""}

    files = list(m["files"])
    # Re-checked rather than trusted even though the scan just ran: /scratch is under a 21-day
    # purge, and a job that queues behind a long one can start hours after the scan.
    missing = [f.name for f in files if not f.exists()]
    if missing:
        return {"Run": run, "BioProject": proj, "status": f"missing:{missing[0]}",
                "seconds": 0, "report_bytes": 0, "species_rows": "", "reads_classified": ""}

    # Write to a temp name and rename on success, so an interrupted job never leaves a
    # partial report that the resume check would accept.
    tmp = tmp_dir / f"{run}.report.partial"
    cmd = ["kraken2", "--db", str(db), "--threads", str(threads),
           "--confidence", CONFIDENCE, "--minimum-hit-groups", MIN_HIT_GROUPS,
           # Without this the report has 6 columns and no distinct-minimizer count, so every
           # coverage in 05_detect would be None. It is not optional.
           "--report-minimizer-data", "--report", str(tmp), "--output", "/dev/null"]
    if mmap:
        cmd.append("--memory-mapping")
    if len(files) == 2:
        cmd.append("--paired")
    cmd += [str(f) for f in files]

    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as e:                      # noqa: BLE001 - must never kill the pool
        return {"Run": run, "BioProject": proj, "status": f"error:{type(e).__name__}",
                "seconds": round(time.time() - t0), "report_bytes": 0,
                "species_rows": "", "reads_classified": ""}
    el = round(time.time() - t0)
    if r.returncode != 0 or not tmp.exists() or tmp.stat().st_size < MIN_REPORT_BYTES:
        tmp.unlink(missing_ok=True)
        return {"Run": run, "BioProject": proj, "status": "failed", "seconds": el,
                "report_bytes": 0, "species_rows": "",
                "reads_classified": r.stderr.strip()[-200:] if r.stderr else ""}
    tmp.replace(report)

    n_sp, reads = 0, 0
    for rec in parse_report(report):
        if rec["rank"] == "S" and rec["reads_clade"]:
            n_sp += 1
            reads += rec["reads_clade"]
    print(f"  [done ] {run} {el:5d}s  {n_sp:4d} species  {reads:>12,d} reads", flush=True)
    return {"Run": run, "BioProject": proj, "status": "ok", "seconds": el,
            "report_bytes": report.stat().st_size, "species_rows": n_sp,
            "reads_classified": reads}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", help="BioProject to classify; omit with --list")
    ap.add_argument("--list", action="store_true", help="what is on disk and what is done")
    ap.add_argument("--db", default=str(DB / "db_v4"))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--mmap", action="store_true",
                    help="--memory-mapping. Correct ONLY against /dev/shm, never Lustre")
    ap.add_argument("--force", action="store_true", help="reclassify runs that have a report")
    ap.add_argument("--tmp-dir", default=os.environ.get("MYSCRATCH", "/tmp") + "/kraken_tmp")
    # The NEGATIVE CONTROLS live in their own reads tree (`_layout.CONTROL_READS`) so they
    # cannot enter the cohort's denominator -- see the comment beside that constant. The tree
    # has the identical shape, so pointing these two at it is the whole integration; nothing
    # about the classification changes, which is the point, since a control classified under
    # different settings measures nothing.
    #   `--project` then names the control ACCESSION, and the manifest is written beside the
    # reads rather than over the cohort's.
    ap.add_argument("--reads", default=str(READS),
                    help="reads root; use CONTROL_READS for the negative controls")
    ap.add_argument("--reports", default=str(REPORTS), help="reports root; pair it with --reads")
    args = ap.parse_args()

    reads_root, reports_root = Path(args.reads), Path(args.reports)
    man = manifest(reads_root, reads_root.parent / "reads_manifest.tsv")
    by_proj = {}
    for run, m in man.items():
        by_proj.setdefault(m["BioProject"], []).append(run)

    if args.list or not args.project:
        print(f"{'BioProject':14s} {'onDisk':>7s} {'done':>6s} {'todo':>6s} {'TiB':>7s}")
        tot = [0, 0, 0]
        for proj in sorted(by_proj, key=lambda p: -len(by_proj[p])):
            runs = by_proj[proj]
            d = sum(1 for r in runs if done(reports_root / proj / f"{r}.report.tsv"))
            b = sum(int(man[r]["bytes"]) for r in runs) / 2**40
            print(f"{proj:14s} {len(runs):7d} {d:6d} {len(runs)-d:6d} {b:7.3f}")
            tot[0] += len(runs); tot[1] += d; tot[2] += len(runs) - d
        print(f"{'TOTAL':14s} {tot[0]:7d} {tot[1]:6d} {tot[2]:6d}")
        if not args.project:
            return

    runs = sorted(by_proj.get(args.project, []))
    if not runs:
        sys.exit(f"no reads on disk for {args.project}")
    if args.force:
        for r in runs:
            (reports_root / args.project / f"{r}.report.tsv").unlink(missing_ok=True)
    todo = [r for r in runs if not done(reports_root / args.project / f"{r}.report.tsv")]
    print(f"\n{args.project}: {len(runs)} runs on disk, {len(todo)} to classify")
    if not todo:
        print("nothing to do")
        return

    tmp_dir = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(classify_run, r, man[r], Path(args.db), args.threads,
                          args.mmap, tmp_dir, reports_root): r for r in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            rows.append(fut.result())
            if i % 25 == 0:
                el = time.time() - t0
                print(f"  [{i}/{len(todo)}] {el/60:.1f} min, "
                      f"eta {el/i*(len(todo)-i)/60:.0f} min", flush=True)

    # Append, because a project is classified in its own job and the table is cumulative.
    SILVER.mkdir(parents=True, exist_ok=True)
    new = not STATUS.exists()
    with open(STATUS, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_COLS, delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
        w.writerows(rows)

    ok = sum(1 for r in rows if r["status"] == "ok")
    bad = [r for r in rows if r["status"] not in ("ok", "cached")]
    print(f"\n{args.project}: {ok} classified, {len(bad)} failed, "
          f"{(time.time()-t0)/60:.1f} min")
    for r in bad[:10]:
        print(f"  FAILED {r['Run']}: {r['status']} {r['reads_classified'][:120]}")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
