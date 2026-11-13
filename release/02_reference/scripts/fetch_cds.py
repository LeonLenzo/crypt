#!/usr/bin/env python3
"""Download CDS (coding sequences only) for the assemblies in a build tier.

    python release/02_reference/scripts/fetch_cds.py --tier db0_baseline
    python release/02_reference/scripts/fetch_cds.py --tier db0_baseline --have acacia.txt
    python release/02_reference/scripts/fetch_cds.py --tier db0_baseline,db1_congener --dry-run

Reads  data/outputs/build_plan.tsv       the tiered assembly selection from plan_build.py
       --have FILE (optional)            accessions already available (one per line); skipped
Writes data/raw/cds/<accession>/cds_from_genomic.fna   the CDS, one dir per accession
       data/raw/cds_fetch_manifest.tsv   accession -> status, bytes, TRACKED

CDS ONLY. `datasets download genome --include cds` pulls cds_from_genomic.fna and nothing
else, which is the extraction of the submitter's annotation, not the genome. An assembly with
no gene model has no CDS to pull; plan_build only lists annotated assemblies, so that case
does not arise here.

Idempotent: an accession whose cds_from_genomic.fna is already present and non-empty is
skipped, so a killed run resumes by re-running. --have takes the accessions already on Acacia
(from `datasets`-built db_v4) so only the genuine gap is fetched; produce it with
`rclone lsf acacia:bucket/db/search/cds/pathogen` or from cds_manifest.tsv.

--workers caps concurrency at 8. NCBI throttles and the recurring failure in this project is
a dropped connection mid-package, not rate-limiting; 8 parallel single-accession packages
retries cleanly where one 700-accession package would restart whole.
"""

import argparse
import csv
import shutil
import subprocess
import sys
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
PLAN = HERE / "data/outputs/build_plan.tsv"
CDS = HERE / "data/raw/cds"
MANIFEST = HERE / "data/raw/cds_fetch_manifest.tsv"


def read_have(path):
    if not path:
        return set()
    return {ln.strip().split()[0].split("/")[0] for ln in open(path)
            if ln.strip() and not ln.startswith("#")}


def cds_path(acc):
    return CDS / acc / "cds_from_genomic.fna"


def fetch_one(acc, retries=3):
    """Download one accession's CDS into CDS/<acc>/. Returns (acc, status, bytes)."""
    dest = cds_path(acc)
    if dest.exists() and dest.stat().st_size > 0:
        return acc, "present", dest.stat().st_size
    dest.parent.mkdir(parents=True, exist_ok=True)
    zpath = dest.parent / "cds.zip"
    for attempt in range(1, retries + 1):
        try:
            r = subprocess.run(
                ["datasets", "download", "genome", "accession", acc,
                 "--include", "cds", "--filename", str(zpath)],
                capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                if attempt == retries:
                    return acc, f"datasets_error:{r.stderr.strip()[:60]}", 0
                continue
            with zipfile.ZipFile(zpath) as z:
                member = next((m for m in z.namelist()
                               if m.endswith("cds_from_genomic.fna")), None)
                if not member:
                    return acc, "no_cds_in_package", 0
                with z.open(member) as src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)
            zpath.unlink(missing_ok=True)
            return acc, "fetched", dest.stat().st_size
        except subprocess.TimeoutExpired:
            if attempt == retries:
                return acc, "timeout", 0
        except (zipfile.BadZipFile, OSError) as e:
            if attempt == retries:
                return acc, f"error:{type(e).__name__}", 0
    return acc, "failed", 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tier", default="db0_baseline",
                    help="comma-separated db tiers to fetch (default db0_baseline)")
    ap.add_argument("--have", help="file of accessions already available, to skip")
    ap.add_argument("--workers", type=int, default=8, help="parallel downloads (<=8)")
    ap.add_argument("--dry-run", action="store_true", help="list the gap, fetch nothing")
    args = ap.parse_args()

    tiers = set(args.tier.split(","))
    have = read_have(args.have)
    want = []
    for r in csv.DictReader(open(PLAN), delimiter="\t"):
        if r["db_tier"] in tiers and r["accession"] not in have:
            want.append(r["accession"])
    want = sorted(set(want))
    already = sum(1 for a in want if cds_path(a).exists()
                  and cds_path(a).stat().st_size > 0)
    print(f"tiers {sorted(tiers)}: {len(want):,} accessions after --have "
          f"({len(have):,} skipped as already on hand)")
    print(f"  {already:,} already downloaded locally, {len(want) - already:,} to fetch")
    if args.dry_run:
        for a in want:
            if not (cds_path(a).exists() and cds_path(a).stat().st_size > 0):
                print(f"  would fetch {a}")
        return

    results = []
    with ThreadPoolExecutor(max_workers=min(args.workers, 8)) as ex:
        futs = {ex.submit(fetch_one, a): a for a in want}
        done = 0
        for fut in as_completed(futs):
            acc, status, nbytes = fut.result()
            results.append((acc, status, nbytes))
            done += 1
            if status not in ("present", "fetched") or done % 50 == 0:
                print(f"  [{done}/{len(want)}] {acc} {status}")

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "status", "bytes"])
        for acc, status, nbytes in sorted(results):
            w.writerow([acc, status, nbytes])
    ok = sum(1 for _, s, _ in results if s in ("fetched", "present"))
    bad = [(a, s) for a, s, _ in results if s not in ("fetched", "present")]
    gb = sum(n for _, _, n in results) / 1e9
    print(f"\n  {ok:,}/{len(results):,} CDS present, {gb:.1f} GB")
    if bad:
        print(f"  {len(bad)} failed (re-run to retry): "
              + ", ".join(f"{a}:{s}" for a, s in bad[:8]))
    print(f"  -> {MANIFEST}")


if __name__ == "__main__":
    main()
