#!/usr/bin/env python3
"""Inventory the cohort's sequencing reads, and fetch the ones not already on disk.

    python release/01_cohort/scripts/fetch_reads.py                         # inventory only
    python release/01_cohort/scripts/fetch_reads.py --reads-dir DIR         # mark what DIR already holds
    python release/01_cohort/scripts/fetch_reads.py --reads-dir DIR --fetch # pull the missing ones

Reads  data/curation/screen_samples.tsv   the cohort: one row per run
Writes data/raw/reads_ena/<project>.tsv   ENA's file listing per project (fetched, cached)
       data/raw/reads_manifest.tsv        one row per run: expected files, md5, bytes, local status

WHY THE FETCH LIVES HERE, NOT IN 02_classify. The reads are public SRA/ENA data the cohort
points to; materialising them is part of building the cohort, not part of the classifier (the
database is the method). So release owns the read manifest and research consumes it.

ENA, NOT prefetch + fasterq-dump. ENA serves every library as ready-made fastq.gz, so fetching
is network I/O on one core rather than CPU-bound .sra decompression; and ENA publishes an md5
per file, so integrity is checked without hashing terabytes. Every cohort run has an ENA
fastq.gz (verified across the pilot set).

PRESENCE IS A FILE ON DISK, never a row in this manifest: the manifest is the thing being
built, so trusting it to decide what exists would let a stale row suppress a real fetch.
`--reads-dir` is scanned run by run; `--inventory` seeds presence from a prior manifest when no
tree is reachable (a bootstrap, not a substitute for the scan). Nothing already complete and
md5-clean is refetched.
"""

import argparse
import csv
import hashlib
import collections
import concurrent.futures as cf
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
DATA = MODULE / "data"
SCREEN = DATA / "curation" / "screen_samples.tsv"
ENADIR = DATA / "raw" / "reads_ena"
OUT = DATA / "raw" / "reads_manifest.tsv"
PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
UA = {"User-Agent": "crypt/1.0 (+mailto:leon.lenzo@curtin.edu.au)"}
COLS = ["run", "BioProject", "n_files", "files", "ena_bytes", "ena_md5",
        "local_bytes", "status"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_manifest(rows):
    """Write the manifest. Called BEFORE the fetch, periodically during it, and at the end.

    It used to be written only after the pool drained, which meant a job killed at its wall
    clock left the manifest describing the PREVIOUS run. That is how `reads_manifest.tsv`
    came to be dated 2026-11-13 while 875 runs sat on disk, and it is why the state had to
    be counted off the filesystem instead. Written through a temp file and renamed, so a
    reader never sees a half-written table.
    """
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".tsv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(OUT)


def ena_files(project, refresh=False):
    """ENA's fastq listing for a project: {run: [(filename, url, md5, bytes), ...]}, cached."""
    cache = ENADIR / f"{project}.tsv"
    if cache.exists() and not refresh:
        rows = read(cache)
    else:
        url = (f"{PORTAL}?accession={project}&result=read_run"
               "&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv")
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as fh:
            text = fh.read().decode()
        ENADIR.mkdir(parents=True, exist_ok=True)
        cache.write_text(text, encoding="utf-8")
        rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))

    out = {}
    for r in rows:
        ftp = [x for x in (r.get("fastq_ftp") or "").split(";") if x]
        md5 = (r.get("fastq_md5") or "").split(";")
        size = (r.get("fastq_bytes") or "").split(";")
        files = [(u.rsplit("/", 1)[-1], "https://" + u, md5[i] if i < len(md5) else "",
                  int(size[i]) if i < len(size) and size[i] else 0)
                 for i, u in enumerate(ftp)]
        # ENA's orphan-mate file (<run>.fastq.gz beside _1/_2) would make n_files meaningless;
        # when a pair is present the paired files win and the bare file is dropped.
        if any("_1.fastq" in f or "_2.fastq" in f for f, *_ in files):
            files = [f for f in files if "_1.fastq" in f[0] or "_2.fastq" in f[0]]
        out[r["run_accession"]] = files
    return out


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


class ShortRead(Exception):
    """The stream ended early. The bytes on disk are GOOD, there are just not enough yet."""


def fetch_file(url, dst, md5, size, verify, tries=8):
    """Stream one fastq.gz to a .part, RESUMING it, size- and md5-check it, rename into place.

    THE .part SURVIVES A FAILED ATTEMPT, WHICH IS THE WHOLE POINT. EBI closes these
    connections cleanly partway through a large file -- measured on 2026-11-19 at 201 MB,
    316 MB, 402 MB, 432 MB and 469 MB into files of 1.7-2.5 GB, so the cut point is arbitrary
    and always well short. The old code deleted the .part and restarted from zero on every
    failure, which means a file EBI will not serve in one go can never finish no matter how
    many retries it gets: three attempts from zero is three truncations. Job 53279097 spent
    eight hours proving that 21 times.

    EBI advertises `Accept-Ranges: bytes` and answers 206, so each attempt resumes from where
    the last one stopped and every truncation becomes progress instead of a wasted transfer.

    A SHORT READ AND A BAD CHECKSUM GET OPPOSITE TREATMENT. Missing bytes are resumed; wrong
    bytes cannot be, so an md5 failure over a file of the RIGHT length deletes the .part and
    starts clean. Keeping it would resume forever on top of corruption.

    `tries` is 8 rather than 4 because an attempt now advances the file instead of replacing
    it, so more attempts are worth having.
    """
    tmp = dst.with_suffix(dst.suffix + ".part")
    for attempt in range(1, tries + 1):
        try:
            have = tmp.stat().st_size if tmp.exists() else 0
            if size and have > size:
                # Longer than ENA says: it cannot be a prefix of the right file.
                tmp.unlink(missing_ok=True)
                have = 0
            req = urllib.request.Request(url, headers=dict(UA))
            if have:
                req.add_header("Range", f"bytes={have}-")
            with urllib.request.urlopen(req, timeout=600) as r:
                # A server that ignores Range replies 200 with the WHOLE file; appending that
                # to what we have would silently produce a longer, corrupt file that still
                # fails only at the md5. Append only on an explicit 206.
                resuming = bool(have) and getattr(r, "status", r.getcode()) == 206
                if have and not resuming:
                    have = 0
                with open(tmp, "ab" if resuming else "wb") as fh:
                    while True:
                        block = r.read(1 << 20)
                        if not block:
                            break
                        fh.write(block)
            got = tmp.stat().st_size
            if size and got != size:
                raise ShortRead(f"got {got} of {size} bytes ({got / size * 100:.0f}%)")
            if verify and md5 and md5sum(tmp) != md5:
                tmp.unlink(missing_ok=True)
                raise ValueError("md5 mismatch over a file of the right length: "
                                 "corrupt, restarting clean")
            tmp.rename(dst)
            return dst.stat().st_size
        except Exception as e:
            if attempt == tries:
                raise
            # Say so. A silent retry made a job that was coping look like a job that was hung,
            # and EBI cutting a stream is the normal case here, not a fault: it is why `tries`
            # exists. Report what is banked so a resume reads as progress.
            banked = tmp.stat().st_size if tmp.exists() else 0
            where = f"  [{banked / size * 100:.0f}% banked]" if size and banked else ""
            print(f"    retry {attempt}/{tries - 1}  {dst.name}  "
                  f"{type(e).__name__}: {e}{where}", flush=True)
            time.sleep(5 * attempt)


def _refresh_status(rows, todo, errs, final=True):
    """Recompute status for every queued run from what is NOW on disk. Returns (ok, bad).

    Split out of the end-of-job block so the manifest can be refreshed mid-flight. A run is
    `fetched` only when every one of its files is present.

    `final` separates two states a single on-disk check cannot tell apart. A run missing a
    file because that file ERRORED is `failed: <reason>`. A run missing a file that simply
    has not been reached yet is `pending` mid-flight and only becomes `failed: incomplete`
    once the pool has drained. Without the distinction an interim manifest reports every
    run not yet started as a failure, which is the opposite of informative.
    """
    by_run = {row["run"]: row for row in rows}
    ok = bad = 0
    for run, (dirp, fs) in todo.items():
        row = by_run[run]
        got = [n for n, _, _ in fs if (dirp / n).exists()]
        if len(got) == len(fs) and not errs[run]:
            row["local_bytes"] = sum((dirp / n).stat().st_size for n, _, _ in fs)
            row["status"] = "fetched"
            ok += 1
        elif errs[run]:
            row["status"] = f"failed: {errs[run][0]}"
            bad += 1
        elif final:
            row["status"] = "failed: incomplete"
            bad += 1
        else:
            row["local_bytes"] = sum((dirp / n).stat().st_size for n in got)
            row["status"] = f"pending ({len(got)}/{len(fs)} files)"
    return ok, bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads-dir", type=Path, help="reads tree to scan (reads-dir/<project>/)")
    ap.add_argument("--inventory", type=Path,
                    help="a prior manifest to seed presence from when no tree is reachable")
    ap.add_argument("--fetch", action="store_true", help="download the missing files")
    ap.add_argument("--workers", type=int, default=24,
                    help="concurrent download streams (EBI scales ~linearly to 24)")
    ap.add_argument("--no-verify", action="store_true", help="skip md5 check on fetch")
    ap.add_argument("--refresh", action="store_true", help="re-request ENA listings")
    args = ap.parse_args()

    cohort = read(SCREEN)
    # a prior manifest may spell the run column "run" or "Run"; accept either
    known = {}
    if args.inventory and args.inventory.exists():
        for r in read(args.inventory):
            known[r.get("run") or r.get("Run")] = r

    rows, by_proj = [], collections.defaultdict(list)
    for r in cohort:
        by_proj[r["BioProject"]].append(r["run"])
    listings = {p: ena_files(p, args.refresh) for p in sorted(by_proj)}

    fetched = failed = 0
    todo = {}   # run -> (dirp, [(name,url,md5)...]) for missing runs when --fetch
    for r in cohort:
        run, proj = r["run"], r["BioProject"]
        files = listings.get(proj, {}).get(run, [])
        names = [f[0] for f in files]
        row = {"run": run, "BioProject": proj, "n_files": len(files),
               "files": ";".join(names), "ena_bytes": sum(f[3] for f in files),
               "ena_md5": ";".join(f[2] for f in files), "local_bytes": "", "status": ""}

        dirp = (args.reads_dir / proj) if args.reads_dir else None
        # Count what is actually there rather than testing all-or-nothing, so a run holding
        # one mate of two is reported as PARTIAL instead of being lumped in with runs that
        # have nothing. The two need different readings: a partial run is unusable for
        # analysis but nearly paid for, and on 2026-11-17 there were 175 of them hidden
        # inside a single "missing 685".
        got = [n for n in names if (dirp / n).exists() and (dirp / n).stat().st_size] if dirp \
            else []
        present = bool(dirp) and names and len(got) == len(names)
        if present:
            row["local_bytes"] = sum((dirp / n).stat().st_size for n in names)
            row["status"] = "present"
        elif args.fetch and dirp and names:
            # queued; overwritten by _refresh_status once the pool has run
            row["status"] = f"partial ({len(got)}/{len(names)} files)" if got else "missing"
            if got:
                row["local_bytes"] = sum((dirp / n).stat().st_size for n in got)
            todo[run] = (dirp, [(n, url, md5, sz) for n, url, md5, sz in files])
        elif got:
            row["status"] = f"partial ({len(got)}/{len(names)} files)"
            row["local_bytes"] = sum((dirp / n).stat().st_size for n in got)
        elif not files:
            row["status"] = "no ENA fastq"
        elif run in known:
            row["status"] = "present (inventory)"
            row["local_bytes"] = known[run].get("bytes") or known[run].get("local_bytes") or ""
        elif args.reads_dir or args.inventory:
            row["status"] = "missing"  # a tree or a prior inventory was given and lacks it
        else:
            row["status"] = "unscanned"
        rows.append(row)

    # Download the queued files concurrently: EBI gives each stream its own bandwidth, so
    # throughput scales with workers. One file per task; a run counts fetched only if every one
    # of its files lands. Nothing already on disk is re-requested, so a re-run resumes.
    # Write it once BEFORE fetching anything, so even a job killed in its first minute
    # leaves a manifest describing the CURRENT disk rather than a previous run's.
    write_manifest(rows)
    if todo:
        for dirp, _ in todo.values():
            dirp.mkdir(parents=True, exist_ok=True)
        tasks = [(run, dirp / n, url, md5, sz) for run, (dirp, fs) in todo.items()
                 for n, url, md5, sz in fs if not (dirp / n).exists()]

        # HALF-COMPLETE RUNS FIRST. The commit point is the FILE, not the run, so whenever a
        # file fails its retries or a job is killed at the wall clock, the sibling mate that
        # already landed stays on disk and the run sits unusable with one of two files. On
        # 2026-11-17 that was 148 runs. Ordering by how much of the run is already there
        # finishes those first, which converts half-runs into analysable runs instead of
        # starting new ones. It does not change WHAT gets fetched, only the order.
        have = {run: sum(1 for n, _, _, _ in fs if (dirp / n).exists())
                for run, (dirp, fs) in todo.items()}
        tasks.sort(key=lambda t: -have[t[0]])
        n_partial = sum(1 for t in tasks if have[t[0]])
        print(f"fetching {len(tasks)} files for {len(todo)} runs with {args.workers} workers",
              flush=True)
        if n_partial:
            print(f"  {n_partial} of those complete a run that is already part-downloaded; "
                  f"queued first", flush=True)

        errs = collections.defaultdict(list)
        t0 = time.time()
        done_n = done_b = 0
        with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(fetch_file, url, dst, md5, sz, not args.no_verify): (run, dst)
                    for run, dst, url, md5, sz in tasks}
            for fut in cf.as_completed(futs):
                run, dst = futs[fut]
                try:
                    nb = fut.result() or 0
                    done_n += 1
                    done_b += nb
                    el = max(time.time() - t0, 1e-9)
                    rate = done_b / el
                    left = len(tasks) - done_n
                    eta = (left * (done_b / done_n) / rate / 3600) if done_n and rate else 0
                    print(f"  [{done_n}/{len(tasks)}] {dst.name}  {nb/1e9:.2f} GB  "
                          f"{done_b/1e12:.3f} TB total  {rate/1e6:.1f} MB/s  "
                          f"ETA {eta:.1f} h", flush=True)
                except Exception as e:
                    errs[run].append(str(e))
                    done_n += 1
                    print(f"  [{done_n}/{len(tasks)}] FAILED {dst.name}  {e}", flush=True)
                # Keep the manifest honest while the job runs, not only if it survives to
                # the end. Every 25 files is cheap against ~1.5 GB per file.
                if done_n % 25 == 0:
                    _refresh_status(rows, todo, errs, final=False)
                    write_manifest(rows)
        fetched, failed = _refresh_status(rows, todo, errs)

    write_manifest(rows)

    st = collections.Counter(r["status"].split(":")[0].split(" (")[0] for r in rows)
    tb = sum(int(r["ena_bytes"] or 0) for r in rows) / 1e12
    print(f"-> {OUT.relative_to(MODULE)}  {len(rows)} runs over {len(by_proj)} projects, "
          f"{tb:.2f} TB on ENA")
    for s in ("present", "present (inventory)", "fetched", "partial", "pending", "missing",
              "unscanned", "no ENA fastq", "failed"):
        if st[s]:
            print(f"   {s:<20}{st[s]:>5}")
    if args.fetch:
        print(f"   fetched {fetched}, failed {failed} this run")


if __name__ == "__main__":
    main()
