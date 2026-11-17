#!/usr/bin/env python3
"""One-off: move the fungal CDS out of research/02_classify into release/02_reference.

    python3 tools/migrate_cds_to_release.py --dry-run
    sbatch tools/migrate_cds.slurm

NOT part of the release pipeline, which is why it lives in tools/ and not in
release/02_reference/scripts/. A reviewer never runs this: they run fetch_cds.py, which
downloads the same CDS from NCBI. This exists once, because the bytes were fetched before the
release rewrite existed and re-downloading 92 GB to move it is absurd.

    from  research/02_classify/data/bronze/cds_db0/pathogen/<acc>/ncbi_dataset/data/<acc>/cds_from_genomic.fna
    to    release/02_reference/data/raw/cds/<acc>/cds_from_genomic.fna

## What it is actually doing

RENAME, not copy. Both trees are on the same Lustre filesystem, so each accession is a
metadata operation and the 92 GB never moves. That also makes it reversible: nothing is
deleted and no bytes are rewritten, so a wrong destination costs another rename.

**The layout FLATTENS.** The old tree kept the `datasets` package shape
(`<acc>/ncbi_dataset/data/<acc>/`), which is three directories carrying no information; the
release layout is `<acc>/cds_from_genomic.fna`, which is what release/02_reference's own
`fetch_cds.py` writes. Verified before writing this: an accession directory holds that one
file and nothing else, so flattening loses nothing.

**db0 and db1 are ONE tree.** `cds_db1/pathogen` is empty and db1's manifest resolves into
`cds_db0/pathogen`; all 3,533 db0 and 5,272 db1 accessions are there, in 5,302 directories.
So this one move covers both databases, and the 30 directories in neither manifest come
across too rather than being judged here.

**The built indexes are NOT touched.** `research/02_classify/data/silver/db/{db0,db1,db_v4}`
stay where they are: db1 is about to be superseded by db1' (= db1 + the host sink), which is a
fresh kraken2-build writing into the release tree, and `db1/library/` is 101 GB of tagged
copies that the rebuild regenerates. Delete the old indexes once db1' validates, not now.

## The manifest is the point

The CDS tree is gitignored (92 GB on a filesystem under a 21-day purge), so what is present
has to be recorded in a tracked file or it is not recorded at all. This writes
`release/02_reference/data/raw/cds_fetch_manifest.tsv` in exactly the three-column shape
`fetch_cds.py` writes, from a DISK SCAN rather than from the old manifests: the old ones
describe what a fetch intended, and the question a reviewer asks is what the build actually
read.

Also reported, because it is the next action: accessions in `build_plan.tsv` with no CDS on
disk (what fetch_cds.py still owes db1') and accessions on disk that no longer appear in the
plan.

Idempotent. An accession already at its destination is skipped, so a killed run resumes.
"""

import argparse
import csv
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "research/02_classify/data/bronze/cds_db0/pathogen"
DEST = REPO / "release/02_reference/data/raw/cds"
PLAN = REPO / "release/02_reference/data/curation/build_plan.tsv"
MANIFEST = REPO / "release/02_reference/data/raw/cds_fetch_manifest.tsv"
FNA = "cds_from_genomic.fna"


def rel(path):
    """Repo-relative for readability, absolute when the path is outside it (--src/--dest)."""
    try:
        return path.relative_to(REPO)
    except ValueError:
        return path


def source_file(acc_dir):
    """The one CDS file inside a `datasets` package directory, or None."""
    direct = acc_dir / FNA
    if direct.exists():
        return direct
    hits = sorted(acc_dir.glob(f"ncbi_dataset/data/*/{FNA}"))
    if len(hits) > 1:
        raise SystemExit(f"{acc_dir.name}: {len(hits)} CDS files in one package, expected 1")
    return hits[0] if hits else None


def prune(path, stop):
    """Remove empty directories upward from `path`, never past `stop`."""
    p = path
    while p != stop and p.is_dir():
        try:
            p.rmdir()
        except OSError:
            return
        p = p.parent


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--src", type=Path, default=SRC)
    ap.add_argument("--dest", type=Path, default=DEST)
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would move, rename nothing")
    ap.add_argument("--manifest-only", action="store_true",
                    help="rescan the destination and rewrite the manifest, move nothing")
    args = ap.parse_args()

    if not args.manifest_only:
        if not args.src.is_dir():
            raise SystemExit(f"source tree is not there: {args.src}")
        accs = sorted(d.name for d in args.src.iterdir() if d.is_dir())
        print(f"{len(accs):,} accession directories in {rel(args.src)}")

        moved = skipped = 0
        empty, mismatch = [], []
        for acc in accs:
            dst = args.dest / acc / FNA
            if dst.exists() and dst.stat().st_size:
                skipped += 1
                continue
            src = source_file(args.src / acc)
            if src is None or not src.stat().st_size:
                empty.append(acc)
                continue
            size = src.stat().st_size
            if args.dry_run:
                moved += 1
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.rename(src, dst)
            except OSError as e:
                raise SystemExit(f"{acc}: rename failed ({e}). Same filesystem? "
                                 "A cross-device move needs a copy, not this script")
            if dst.stat().st_size != size:
                mismatch.append(acc)
                continue
            prune(src.parent, args.src)
            moved += 1

        verb = "would move" if args.dry_run else "moved"
        print(f"  {verb} {moved:,}, {skipped:,} already in place")
        if empty:
            print(f"  {len(empty)} with NO CDS file, left alone: "
                  + ", ".join(empty[:8]) + (" ..." if len(empty) > 8 else ""))
        if mismatch:
            raise SystemExit(f"{len(mismatch)} changed size across the rename: {mismatch[:8]}")
        if args.dry_run:
            return

    # The manifest, from the disk and not from the old tables.
    on_disk = {}
    for d in sorted(args.dest.iterdir()) if args.dest.is_dir() else []:
        f = d / FNA
        if d.is_dir() and f.exists() and f.stat().st_size:
            on_disk[d.name] = f.stat().st_size
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "status", "bytes"])
        for acc, nbytes in sorted(on_disk.items()):
            w.writerow([acc, "present", nbytes])
    gb = sum(on_disk.values()) / 1e9
    print(f"\n  {len(on_disk):,} CDS present, {gb:.1f} GB")
    print(f"  -> {rel(MANIFEST)}")

    # What the move did and did not settle, which is the next action either way.
    if PLAN.exists():
        planned = {}
        with open(PLAN) as fh:
            for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                    delimiter="\t"):
                # The host tier is NOT in this tree. Host CDS come from Ensembl Plants,
                # land in raw/host_cds/ and are recorded in host_cds_manifest.tsv, so
                # counting their GenBank accessions here would report 21 permanent failures.
                if r.get("accession") and r.get("db_tier") != "host":
                    planned.setdefault(r.get("db_tier", "?"), set()).add(r["accession"])
        want = set().union(*planned.values()) if planned else set()
        missing = sorted(want - set(on_disk))
        orphan = sorted(set(on_disk) - want)
        for tier in sorted(planned):
            gap = len(planned[tier] - set(on_disk))
            print(f"  {tier:14s} {len(planned[tier]):>6,} planned, {gap:>5,} with no CDS")
        print(f"  {len(missing):,} planned accessions still to fetch, "
              f"{len(orphan):,} on disk outside the plan")
        if missing:
            print("    next: fetch_cds.py --tier <tier>  (sbatch, copy partition)")


if __name__ == "__main__":
    main()
