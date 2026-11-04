#!/usr/bin/env python3
"""Inventory the cohort's sequencing reads, and fetch the ones not already on disk.

    python release/01_cohort/fetch_reads.py                         # inventory only
    python release/01_cohort/fetch_reads.py --reads-dir DIR         # mark what DIR already holds
    python release/01_cohort/fetch_reads.py --reads-dir DIR --fetch # pull the missing ones

Reads  curation/screen_samples.tsv   the cohort: one row per run
Writes raw/reads_ena/<project>.tsv   ENA's file listing per project (fetched, cached)
       raw/reads_manifest.tsv        one row per run: expected files, md5, bytes, local status

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
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCREEN = HERE / "curation" / "screen_samples.tsv"
ENADIR = HERE / "raw" / "reads_ena"
OUT = HERE / "raw" / "reads_manifest.tsv"
PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
UA = {"User-Agent": "crypt/1.0 (+mailto:leon.lenzo@curtin.edu.au)"}
COLS = ["run", "BioProject", "n_files", "files", "ena_bytes", "ena_md5",
        "local_bytes", "status"]


def read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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


def fetch_file(url, dst, md5, verify):
    tmp = dst.with_suffix(dst.suffix + ".part")
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=600) as r, \
            open(tmp, "wb") as fh:
        while True:
            block = r.read(1 << 20)
            if not block:
                break
            fh.write(block)
    if verify and md5 and md5sum(tmp) != md5:
        tmp.unlink()
        raise ValueError("md5 mismatch")
    tmp.rename(dst)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads-dir", type=Path, help="reads tree to scan (reads-dir/<project>/)")
    ap.add_argument("--inventory", type=Path,
                    help="a prior manifest to seed presence from when no tree is reachable")
    ap.add_argument("--fetch", action="store_true", help="download the missing files")
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
    for r in cohort:
        run, proj = r["run"], r["BioProject"]
        files = listings.get(proj, {}).get(run, [])
        names = [f[0] for f in files]
        row = {"run": run, "BioProject": proj, "n_files": len(files),
               "files": ";".join(names), "ena_bytes": sum(f[3] for f in files),
               "ena_md5": ";".join(f[2] for f in files), "local_bytes": "", "status": ""}

        dirp = (args.reads_dir / proj) if args.reads_dir else None
        present = bool(dirp) and names and all((dirp / n).exists() and (dirp / n).stat().st_size
                                               for n in names)
        if present:
            row["local_bytes"] = sum((dirp / n).stat().st_size for n in names)
            row["status"] = "present"
        elif args.fetch and dirp and names:
            dirp.mkdir(parents=True, exist_ok=True)
            try:
                for n, url, md5, _ in files:
                    if not (dirp / n).exists():
                        fetch_file(url, dirp / n, md5, not args.no_verify)
                row["local_bytes"] = sum((dirp / n).stat().st_size for n in names)
                row["status"], fetched = "fetched", fetched + 1
            except Exception as e:
                row["status"], failed = f"failed: {e}", failed + 1
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

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    st = collections.Counter(r["status"].split(":")[0] for r in rows)
    tb = sum(int(r["ena_bytes"] or 0) for r in rows) / 1e12
    print(f"-> {OUT.relative_to(HERE)}  {len(rows)} runs over {len(by_proj)} projects, "
          f"{tb:.2f} TB on ENA")
    for s in ("present", "present (inventory)", "fetched", "missing", "unscanned",
              "no ENA fastq", "failed"):
        if st[s]:
            print(f"   {s:<20}{st[s]:>5}")
    if args.fetch:
        print(f"   fetched {fetched}, failed {failed} this run")


if __name__ == "__main__":
    main()
