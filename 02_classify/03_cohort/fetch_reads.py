#!/usr/bin/env python3
"""fetch_reads.py — pull a BioProject's reads from ENA into reads/<BioProject>/.

    python 02_classify/03_cohort/fetch_reads.py --project PRJEB83049 --pilot-only
    python 02_classify/03_cohort/fetch_reads.py --project PRJEB83049 --dry-run

## ENA, not prefetch + fasterq-dump, and the reason is CPU

The old route cost real compute: jobs 48206159-48206164 ran ~140 hours at 72 allocated cores,
roughly 10,000 SU for 14 TB, because `fasterq-dump` decompressing an .sra and `pigz`
recompressing the output are both CPU-bound. ENA serves the same libraries as ready-made
`fastq.gz`, so fetching is pure network I/O and needs one core.

Checked before adopting it: **all 2,459 pilot runs have an ENA fastq.gz. 100%, no exceptions.**
Validated against 400 runs already on scratch, ENA's bytes are a uniform 0.937 of ours (range
0.907-0.948, no scatter) with ZERO mate-count disagreements. A consistent 6% offset across six
projects is an encoding difference -- compression level and header length -- not a content
difference; different data would scatter. Kraken2 reads sequence, not headers.

URLs come from `run_depth.tsv`, which already cached them one request per BioProject. No
further ENA API traffic happens here.

## Three file shapes, and one of them needs a rule

Across the pilot set: 1,687 runs paired, 771 single, and exactly ONE (SRR22142971 in
PRJNA896566) serving `<RUN>.fastq.gz` alongside `_1`/`_2`. That third file is ENA's orphan
reads -- mates whose partner was discarded. Taking it as well would give the run three files
and make `n_files` meaningless downstream, so `_1`/`_2` win and the orphan file is skipped
with a line in the status table saying so.

## Temp lives OUTSIDE the project tree

Downloads land in `reads_tmp/<Run>/` and are renamed into `reads/<BioProject>/` only once every
file is complete and `gzip -t` passes. The flat 03_kraken tree accumulated 39 orphaned
`_<RUN>_sra` prefetch directories exactly because temp sat beside the output, and in a
per-project tree a half-finished download looks like progress. Presence is a non-empty file in
its FINAL location, always.

## Verification is on bytes, then on gzip integrity

ENA publishes `fastq_bytes` per file, so a truncated download is detectable without
decompressing. A size mismatch is a retry, not a warning: the 2026-09 pass silently truncated
224 of 3,790 files because a flat timeout cut large runs and the success test was existence
only.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import READS, READS_TMP, RUN_LIST, SILVER  # noqa: E402
from _reads import MIN_FASTQ_BYTES, scan_tree  # noqa: E402

DEPTH = SILVER / "run_depth.tsv"
STATUS = SILVER / "fetch_status.tsv"
STATUS_COLS = ["Run", "BioProject", "status", "seconds", "bytes", "n_files", "note"]

# TWENTY-FOUR, and the "8 not 16" rule does NOT apply here. That rule came from NCBI, where 16
# parallel prefetch streams caused connection contention with no throughput gain. The
# constraint against EBI is different and was measured on 2026-10-20: every worker gets its own
# ~0.7 MB/s regardless of how many are running, which is latency-bound TCP (Perth to Hinxton is
# ~300 ms RTT, and one stream cannot fill that window). Contention is not the limit, so
# throughput scales with streams instead of collapsing. At 8 workers a 3.2 GB ENA file takes
# ~75 minutes and the 0.37 TiB pilot set takes ~8 hours.
# If this ever needs to go further, the next lever is lftp's `pget -n` (present on Setonix,
# aria2c is not), which splits ONE file across parallel segments and helps where a project has
# fewer files than workers.
DEFAULT_WORKERS = 24
RETRIES = 5

# ── THE RETRY POLICY, AND THE 623-RUN FAILURE IT IS FIXING ──────────────────────────────────
# Job 51341245 (PRJEB83049, 2026-10-22) marked 623 of 820 runs `failed`. Almost none of them
# were real failures. The median time-to-failure across all 623 was 13 SECONDS, and so were
# the minimum, the 25th and the 75th percentile. Thirteen seconds is the three backoff sleeps
# (2+4+6 = 12 s) and essentially nothing else, so curl exited INSTANTLY on all three attempts.
# The files were fine: ENA still served them with matching `fastq_bytes` the next morning, and
# the same curl flags fetched 5 MB in under a second. EBI had refused connections for a while,
# after ~13 h of 24 sustained streams.
#
# Twelve seconds of total backoff cannot ride out anything a public archive does. The schedule
# below is sized to survive a ~15 minute interruption, which is the realistic unit:
#
#     30s, 60s, 120s, 240s   =  7.5 min of waiting across 5 attempts
#
# and codes that mean "come back later" wait DOUBLE that, so a rate-limit episode gets ~15 min.
# This is cheap insurance: the alternative is re-queueing a 14-hour job.
BACKOFF_BASE = 30
BACKOFF_CAP = 600

# curl exit 22 means HTTP >= 400 under --fail, which lumps "slow down" in with "gone".
# `-w %{http_code}` separates them, and the separation is the whole point: retrying a 404 is
# waste, and abandoning a 503 after 12 s is how 623 runs were lost. The old code could not
# tell the two apart because `capture_output=True` swallowed curl's stderr and nothing ever
# printed it, which is why `fetch_51341245.err` was EMPTY for a 623-run failure and the cause
# had to be inferred from timing. Every terminal failure now prints the reason.
HTTP_RETRY = {408, 429, 500, 502, 503, 504}   # transient; come back later
HTTP_FATAL = {400, 401, 403, 404, 410}        # retrying cannot help
# 416 is NOT an error here. `-C -` on a file that is ALREADY COMPLETE asks for a range past the
# end of it, and the server is right to refuse. The previous comment claimed `--fail` "handled"
# this; it does the opposite, turning a finished download into exit 22. Treated as "the bytes
# are already here" and allowed to fall through to verification.
HTTP_COMPLETE = 416

# ── THE TIMEOUT, AND THE BUG IT IS FIXING ────────────────────────────────────────────────────
# This divisor was 2 MB/s until 2026-10-20 and that number was wrong, in direct contradiction
# of the DEFAULT_WORKERS comment eight lines above, which already recorded the MEASURED
# per-stream rate as ~0.7 MB/s. Both were written in the same sitting. The measurement was
# right; the budget was a guess that never got reconciled with it.
#
# What it cost: 54 of 80 pilot runs failed, and the failure was perfectly size-dependent,
# because that is what a too-tight timeout looks like.
#
#     largest file   ok  fail        budget   needed at 0.7 MB/s
#        <0.5 GB     51     0          900s        408s
#       0.5-1 GB     10     1          900s       1021s   over
#       1-1.5 GB      1     9          900s       1702s   over
#       1.5-2 GB      6    36         1134s       2384s   over
#          >2 GB      0    28         1415s       3187s   over
#
# Re-measured against ENA the same day: 0.69 MB/s on a single stream. Not a transient. The
# files were fine throughout (HTTP 200, Content-Length matching the cached byte count, range
# requests served), so nothing was wrong except the clock.
#
# ── WHY THE FIX IS NOT JUST A BIGGER NUMBER ──────────────────────────────────────────────────
# A max-time sized from an assumed rate fails the same way whenever the assumption drifts, and
# it cannot tell a slow transfer from a dead one. curl can tell them apart directly:
#
#   --speed-limit/--speed-time   abort only if throughput stays under the floor for that long.
#                                This is the REAL timeout now: a genuine stall dies in 3 min
#                                instead of 24, and a slow-but-progressing transfer never dies.
#   --max-time                   a backstop only, sized at 0.25 MB/s. That is ~3x below the
#                                measured rate, so it should never be what fires. If it does,
#                                something is wrong that a retry will not fix.
#
# 50 KB/s is a stall floor, not a rate estimate: it is 1/14th of measured, so only a transfer
# that has effectively stopped trips it.
STALL_BYTES = 50 * 1024      # --speed-limit
STALL_SECONDS = 180          # --speed-time
BACKSTOP_BPS = 0.25 * 1024 * 1024


def budget(nbytes: int) -> int:
    """Backstop max-time. Deliberately generous; stall detection is the real timeout."""
    return max(1800, int(nbytes / BACKSTOP_BPS) + 300)


def backoff(attempt: int, http: int) -> int:
    """Seconds to wait before the next attempt. See THE RETRY POLICY above."""
    wait = min(BACKOFF_CAP, BACKOFF_BASE * (2 ** attempt))
    return wait * 2 if http in HTTP_RETRY else wait


def pick_files(urls: list, run: str) -> tuple[list, str]:
    """Choose which ENA files to take. Returns (urls, note)."""
    paired = [u for u in urls if u.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    if len(paired) == 2:
        note = ""
        if len(urls) > 2:
            note = "skipped ENA orphan-reads file alongside _1/_2"
        return sorted(paired), note
    return urls, ""


def fetch_run(run: str, proj: str, urls: list, sizes: list, dry: bool) -> dict:
    dest = READS / proj
    chosen, note = pick_files(urls, run)
    names = [u.rsplit("/", 1)[-1] for u in chosen]
    final = [dest / n for n in names]
    if all(f.exists() and f.stat().st_size >= MIN_FASTQ_BYTES for f in final):
        return {"Run": run, "BioProject": proj, "status": "cached", "seconds": 0,
                "bytes": sum(f.stat().st_size for f in final), "n_files": len(final),
                "note": note}
    if dry:
        return {"Run": run, "BioProject": proj, "status": "would-fetch", "seconds": 0,
                "bytes": sum(sizes[:len(chosen)]), "n_files": len(chosen), "note": note}

    want = {u.rsplit("/", 1)[-1]: s for u, s in zip(urls, sizes)}
    tmp = READS_TMP / run
    tmp.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        got = []
        for u, name in zip(chosen, names):
            out = tmp / name
            exp = want.get(name, 0)
            url = u if u.startswith("http") else "https://" + u
            verified, why, http = False, "", 0
            for attempt in range(RETRIES):
                # RESUME on a transport failure. ENA sends Accept-Ranges: bytes, so a retry
                # continues from the partial instead of re-downloading from zero. At 0.7 MB/s a
                # 2.3 GB file is ~55 minutes, so restarting a retry is most of an hour thrown
                # away. On a fresh attempt there is no partial, so `-C -` is safe
                # unconditionally; on an already-complete one it yields 416, see HTTP_COMPLETE.
                cmd = ["curl", "-sSL", "--fail",
                       "--speed-limit", str(STALL_BYTES),
                       "--speed-time", str(STALL_SECONDS),
                       "--max-time", str(budget(exp)),
                       "-w", "%{http_code}",
                       "-C", "-", "-o", str(out), url]
                r = subprocess.run(cmd, capture_output=True, text=True)
                # %{http_code} is "000" when curl never got a response at all (DNS, refused
                # connection, stall), which is exactly the case that must stay retryable.
                http = int((r.stdout or "").strip()[-3:] or 0)

                if r.returncode == 0 or http == HTTP_COMPLETE:
                    # Size first: ENA publishes the byte count, so truncation is detectable
                    # without decompressing the whole file.
                    size = out.stat().st_size if out.exists() else 0
                    if exp and abs(size - exp) > max(1024, exp * 0.001):
                        why = f"size {size} != expected {exp}"
                        # Wrong SIZE is not a stall: resuming would append to a bad file, so
                        # this one restarts clean. Only transport failures resume.
                        out.unlink(missing_ok=True)
                    elif subprocess.run(["gzip", "-t", str(out)],
                                        capture_output=True).returncode != 0:
                        why = "gzip -t failed"
                        out.unlink(missing_ok=True)   # corrupt: restart clean, do not resume
                    else:
                        verified = True
                        break
                else:
                    why = (f"curl exit {r.returncode} http {http}: "
                           f"{(r.stderr or '').strip()[:120]}")
                    if http in HTTP_FATAL:
                        break                         # retrying cannot help
                    # else: transport failure. KEEP the partial so the next attempt resumes.

                if attempt < RETRIES - 1:
                    time.sleep(backoff(attempt, http))

            if not verified:
                # The partial is deliberately LEFT in reads_tmp so a resubmit resumes it
                # instead of starting from zero. That is safe because presence is defined as a
                # non-empty file in its FINAL location, so a parked partial can never be
                # mistaken for a finished run; main() prints how much is sitting there.
                print(f"  [fail ] {run} {name} after {attempt + 1}x: {why}", flush=True)
                return {"Run": run, "BioProject": proj, "status": "failed", "n_files": 0,
                        "seconds": round(time.time() - t0), "bytes": 0,
                        "note": f"{name}: {why}"[:200]}
            got.append(out)

        # Publish: rename into the project tree only once EVERY file is verified.
        dest.mkdir(parents=True, exist_ok=True)
        for p in got:
            os.rename(p, dest / p.name)
        shutil.rmtree(tmp, ignore_errors=True)        # empty now; only the WIN cleans up
        nb = sum((dest / p.name).stat().st_size for p in got)
        el = round(time.time() - t0)
        print(f"  [done ] {run} {el:5d}s {nb/1e9:6.2f} GB "
              f"{nb/1e6/max(el,1):6.1f} MB/s", flush=True)
        return {"Run": run, "BioProject": proj, "status": "ok", "seconds": el,
                "bytes": nb, "n_files": len(got), "note": note}
    except Exception as e:                      # noqa: BLE001 - never kill the pool
        return {"Run": run, "BioProject": proj, "status": f"error:{type(e).__name__}",
                "seconds": round(time.time() - t0), "bytes": 0, "n_files": 0, "note": str(e)[:120]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, help="BioProject to fetch")
    ap.add_argument("--pilot-only", action="store_true",
                    help="only runs marked pilot=yes in run_list.tsv. The pilot is a SEEDED "
                         "RANDOM sample of whole BioSamples; never substitute the runs already "
                         "on disk, which are the STAT-selected ones and are biased toward "
                         "pathogen presence")
    ap.add_argument("--limit", type=int,
                    help="stop after this many runs. SMOKE TESTS ONLY: it takes the first N "
                         "in accession order, which is not a random subset -- accession order "
                         "is submission order. To change the sample SIZE, redraw with "
                         "resolve.py --pilot N, which keeps BioSamples whole and stays seeded")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    plan = [r for r in csv.DictReader(open(RUN_LIST), delimiter="\t")
            if r["BioProject"] == args.project and r["fetch_eligible"] == "yes"]
    if not plan:
        sys.exit(f"nothing fetch-eligible for {args.project}")
    if args.pilot_only:
        plan = [r for r in plan if r["pilot"] == "yes"]
    depth = {r["Run"]: r for r in csv.DictReader(open(DEPTH), delimiter="\t")}

    todo = []
    for r in plan:
        d = depth.get(r["Run"])
        if not d or not d["fastq_urls"]:
            print(f"  {r['Run']}: no ENA fastq, skipped")
            continue
        urls = [u for u in d["fastq_urls"].split(";") if u]
        sizes = [int(s) for s in (d["fastq_bytes"] or "").split(";") if s.strip().isdigit()]
        sizes += [0] * (len(urls) - len(sizes))
        todo.append((r["Run"], urls, sizes))
    if args.limit:
        todo = todo[:args.limit]

    tb = sum(sum(s) for _, _, s in todo) / 2**40
    print(f"{args.project}: {len(todo)} runs, {tb:.3f} TiB"
          f"{'  (DRY RUN)' if args.dry_run else ''}")
    READS_TMP.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(fetch_run, run, args.project, urls, sizes, args.dry_run)
                for run, urls, sizes in todo]
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % 20 == 0 and not args.dry_run:
                el = time.time() - t0
                done = sum(r["bytes"] for r in rows) / 2**30
                print(f"  [{i}/{len(todo)}] {done:.1f} GiB, {el/60:.1f} min, "
                      f"eta {el/i*(len(todo)-i)/60:.0f} min", flush=True)

    if not args.dry_run:
        new = not STATUS.exists()
        with open(STATUS, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=STATUS_COLS, delimiter="\t", lineterminator="\n")
            if new:
                w.writeheader()
            w.writerows(rows)

    ok = sum(1 for r in rows if r["status"] in ("ok", "cached"))
    bad = [r for r in rows if r["status"] not in ("ok", "cached", "would-fetch")]
    gb = sum(r["bytes"] for r in rows) / 2**30
    print(f"\n{args.project}: {ok}/{len(rows)} ok, {len(bad)} failed, {gb:.1f} GiB, "
          f"{(time.time()-t0)/60:.1f} min")
    for r in bad[:10]:
        print(f"  FAILED {r['Run']}: {r['status']} {r['note']}")
    notes = [r for r in rows if r["note"]]
    if notes:
        print(f"  {len(notes)} runs carry a note, e.g. {notes[0]['Run']}: {notes[0]['note']}")
    parked = sum(f.stat().st_size for f in READS_TMP.rglob("*") if f.is_file())
    if parked:
        print(f"  {parked / 2**30:.1f} GiB of partials parked in {READS_TMP.name}/ -- a "
              f"resubmit RESUMES these, it does not refetch them")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
