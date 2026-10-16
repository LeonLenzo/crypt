#!/usr/bin/env python3
"""build_db.py — build db_v4 with kraken2-build, from the CDS that fetch_cds.py put on disk.

Step 3 of 3 in 02_build (fetch -> saturate -> build). Fetches nothing. Chooses nothing. It
reads `cds_manifest.tsv`, tags each FASTA with its taxid, adds it to the Kraken2 library and
builds. Everything in the manifest goes in: there is no per-species ceiling and no BUSCO gate
(see fetch_cds.py for why both were rejected, and re-measured).

## Three things the 03_kraken build got wrong, which is why this is a rewrite

**1. It warned instead of failing.** `add_to_library` returned False on error and the caller
printed a WARNING and carried on. db_v3's own log closes with `Selected total: 2115 / Found on
disk: 2081 / Skipped: 32 / Missing: 2` and then builds anyway. Nothing downstream ever learned
that 34 selected assemblies are not in the database, and the shortfall survives only as a line
in a log. Here a missing file or a failed add is FATAL before `--build` is reached, because an
index that is quietly short is worse than no index: every coverage denominator computed from it
is wrong in an unrecorded direction.

**2. It tagged FASTA headers in place, destroying the download.** `tag_fasta_headers` wrote a
`.tagged.fna` then `unlink()`ed the original and renamed over it. That mutates bronze, which
`_layout.py` defines as immutable captures of somebody else's data, and it means a re-build
cannot verify what NCBI actually sent. Here the tagged copy is DERIVED: written to a staging
directory, handed to `--add-to-library` (which copies it into `<db>/library/` regardless), then
deleted. Bronze is never written.

**3. The tagger was idempotent on the PRESENCE of a tag, not its value.** It read the first 30
bytes, saw `>kraken:taxid|`, and returned unchanged. That is a silent-staleness bug and it is
live right now, because the two Acacia buckets disagree:

    personal  pawsey1168-llenzo-kraken-db  2026-09-01  UNTAGGED, pristine NCBI headers
    project   pawsey1168-kraken-db         2026-09-25  ALREADY TAGGED with db_v3 taxids

fetch_cds.py pulls from the PROJECT bucket, so ~1,750 of db_v4's FASTAs arrive carrying a tag
written for a different database. Most are right, which is exactly what makes it dangerous:
`ref_candidates.tsv` assigns taxids through holobase with species-collapse applied, and the
2026-10-15 route-1 fix changed which taxon some accessions belong to. So `tag_fasta()` here
parses the existing tag, compares it to the manifest, re-tags on any difference, and reports
the count. A re-tag is never silent.

Run on Setonix. `--build` needs a high-memory node: db_v3 peaked near 27 GB of hash for 4.8 B
minimizers and db_v4 holds roughly three times the sequence.

    sbatch 02_classify/slurm/build_db.slurm
    python 02_classify/02_build/build_db.py --db-dir /scratch/.../db_v4 --dry-run
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import urllib.request
from collections import Counter
from pathlib import Path

MODULE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import CDS_MANIFEST, CDS_PATHOGEN, DB, DB_INSPECT, LOGS, SILVER  # noqa: E402
from _util import _Tee, atomic_write, link_latest, make_log_dir  # noqa: E402

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAG_PREFIX = "kraken:taxid|"

# Kraken2 defaults, stated rather than inherited so a change is visible in the diff. db_v3 used
# these and the detection thresholds in 05_detect are calibrated against them: k=35, l=31.
# Changing either invalidates every coverage number in the chapter.
KMER_LEN = 35
MINIMIZER_LEN = 31


def load_manifest(path: Path) -> list:
    """The accessions to build from, with their taxids. Refuses an incomplete manifest.

    fetch_cds.py exits non-zero and writes cds_missing.tsv when anything failed, so a
    manifest row with source=MISSING means that exit was ignored. Do not build past it.
    """
    if not path.exists():
        sys.exit(f"no manifest at {path} — run fetch_cds.py first")
    rows = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
    if not rows:
        sys.exit(f"{path} is empty")
    missing = [r for r in rows if r["source"] == "MISSING" or not r["path"]]
    if missing:
        sys.exit(f"{len(missing):,} of {len(rows):,} manifest rows have no CDS on disk "
                 f"(e.g. {', '.join(r['accession'] for r in missing[:3])}).\n"
                 "fetch_cds.py exits non-zero for exactly this; re-run it rather than "
                 "building an index that is silently short.")
    bad = [r for r in rows if not (r["taxid"] or "").strip().isdigit()]
    if bad:
        sys.exit(f"{len(bad):,} manifest rows have no usable taxid, first {bad[0]['accession']}")
    return rows


MIN_FASTA_BYTES = 100_000
DERIVED_SUFFIXES = ("_clean.fna", "_combined.fna", ".tagged.fna")


def find_cds(acc_dir: Path):
    """Fallback locator, using the same rule as fetch_cds.have_cds().

    The manifest path is authoritative; this only covers a relocated tree or a manifest
    written by an earlier version. Largest non-derived `.fna` at or above the size floor, so
    an empty purge skeleton resolves to None rather than to a zero-byte file.
    """
    if not acc_dir.is_dir():
        return None
    cands = [f for f in acc_dir.glob("**/*.fna")
             if f.is_file() and not f.name.endswith(DERIVED_SUFFIXES)
             and f.stat().st_size >= MIN_FASTA_BYTES]
    return max(cands, key=lambda f: f.stat().st_size) if cands else None


def existing_tag(fasta: Path) -> int | None:
    """The taxid already in this FASTA's first header, or None if it is untagged.

    Reads one line, not 30 bytes: the prefix check alone is what let a db_v3 taxid survive
    into a db_v4 build unnoticed.
    """
    with open(fasta, "r", errors="replace") as fh:
        first = fh.readline()
    if not first.startswith(">" + TAG_PREFIX):
        return None
    rest = first[len(">" + TAG_PREFIX):]
    tid = rest.split("|", 1)[0]
    return int(tid) if tid.isdigit() else None


def tag_fasta(src: Path, taxid: int, dest: Path) -> None:
    """Write src to dest with every header tagged `>kraken:taxid|TAXID|...`.

    Tagging lets Kraken2 assign sequence to a taxon directly and skips the 15 GB
    accession2taxid lookup. Any tag already present is REPLACED, not preserved.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    pre = ">" + TAG_PREFIX
    with open(src, "r", errors="replace") as fin, atomic_write(dest) as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                if header.startswith(TAG_PREFIX):
                    # strip `kraken:taxid|<old>|` and keep the original header after it
                    parts = header.split("|", 2)
                    header = parts[2] if len(parts) == 3 else header
                fout.write(f"{pre}{taxid}|{header}\n")
            else:
                fout.write(line)


def add_to_library(fna: Path, db_dir: Path) -> str:
    """kraken2-build --add-to-library. Returns "" on success or the error text.

    `--no-masking` is deliberate: Kraken2's LCA resolves a k-mer shared between species to
    their common ancestor natively, so low-complexity masking buys nothing and removes real
    sequence. This is the db_v2-era decision and it still holds.
    """
    r = subprocess.run(["kraken2-build", "--add-to-library", str(fna),
                        "--db", str(db_dir), "--no-masking"],
                       capture_output=True, text=True)
    return "" if r.returncode == 0 else (r.stderr or r.stdout or "unknown").strip()[:200]


def download_taxonomy(db_dir: Path) -> None:
    tax = db_dir / "taxonomy"
    if (tax / "nodes.dmp").exists() and (tax / "names.dmp").exists():
        print("taxonomy already installed")
        return
    tax.mkdir(parents=True, exist_ok=True)
    tgz = db_dir / "taxdump.tar.gz"
    print(f"downloading taxonomy from {TAXDUMP_URL}", flush=True)
    urllib.request.urlretrieve(TAXDUMP_URL, tgz)
    subprocess.run(["tar", "-xzf", str(tgz), "-C", str(tax)], check=True)
    tgz.unlink(missing_ok=True)
    if not (tax / "nodes.dmp").exists():
        sys.exit("taxdump extracted but nodes.dmp is absent")
    print("taxonomy installed")


def build(db_dir: Path, threads: int) -> None:
    print(f"\nkraken2-build --build  (k={KMER_LEN}, l={MINIMIZER_LEN}, threads={threads})",
          flush=True)
    r = subprocess.run(["kraken2-build", "--build", "--db", str(db_dir),
                        "--kmer-len", str(KMER_LEN), "--minimizer-len", str(MINIMIZER_LEN),
                        "--threads", str(threads)], text=True)
    if r.returncode != 0:
        sys.exit(f"kraken2-build --build failed (exit {r.returncode})")


def inspect(db_dir: Path, out: Path, threads: int) -> None:
    """kraken2-inspect -> the minimizers-per-taxon table 05_detect divides by.

    This is not optional bookkeeping. Coverage is measured against the distinct minimizers a
    taxon holds in THIS database, so without a fresh inspect table every coverage number is
    computed against db_v3's denominators.
    """
    print(f"\nkraken2-inspect → {out}", flush=True)
    r = subprocess.run(["kraken2-inspect", "--db", str(db_dir), "--threads", str(threads)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"kraken2-inspect failed: {(r.stderr or '').strip()[:300]}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write(out) as fh:
        fh.write(r.stdout)
    total = 0
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) > 5 and f[3] == "" and f[5].strip() == "root":
            total = int(f[1])
    print(f"  root minimizers: {total:,}   (db_v3 held 4,800,752,701)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", default=None, help=f"default {CDS_MANIFEST}")
    ap.add_argument("--cds-dir", default=None, help=f"default {CDS_PATHOGEN}")
    ap.add_argument("--db-dir", default=None,
                    help=f"Kraken2 database directory (default {DB}); use scratch on Setonix")
    ap.add_argument("--stage-dir", default=None,
                    help="where tagged copies are written before being added (default "
                         "<db-dir>/.staging). Derived; deleted as it goes.")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--dry-run", action="store_true",
                    help="tag-check and report, add nothing and build nothing")
    ap.add_argument("--build-only", action="store_true",
                    help="skip the library pass; taxonomy + --build + inspect on an existing "
                         "library (for resuming a build that died in --build)")
    args = ap.parse_args()

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "build_db.log")
    link_latest(LOGS, log_dir)
    try:
        man_path = Path(args.manifest) if args.manifest else CDS_MANIFEST
        cds_dir = Path(args.cds_dir) if args.cds_dir else CDS_PATHOGEN
        db_dir = Path(args.db_dir) if args.db_dir else DB
        stage = Path(args.stage_dir) if args.stage_dir else db_dir / ".staging"
        db_dir.mkdir(parents=True, exist_ok=True)

        if not args.build_only:
            rows = load_manifest(man_path)
            print(f"manifest: {len(rows):,} accessions, "
                  f"{len({r['organism_name'] for r in rows}):,} species  ({man_path})")
            print(f"  taxids: {len({r['taxid'] for r in rows}):,} distinct")
            by_source = Counter(r["source"] for r in rows)
            print(f"  provenance: " + ", ".join(f"{k} {v:,}" for k, v in by_source.items()))

            tags = Counter()
            failures, added = [], 0
            for i, r in enumerate(rows, 1):
                acc, taxid = r["accession"], int(r["taxid"])
                # The manifest's `path` is relative to CDS_DIR and already contains the
                # accession directory, because fetch_cds.py writes
                # `f.relative_to(cds_dir)`. Joining cds_dir / acc / path doubles it.
                src = cds_dir / r["path"]
                if not src.is_file():
                    src = find_cds(cds_dir / acc)      # tree moved, or an older manifest
                if src is None or not src.is_file():
                    failures.append((acc, f"no CDS at {cds_dir / r['path']} "
                                          f"and none found under {cds_dir / acc}"))
                    continue
                old = existing_tag(src)
                tags["untagged" if old is None
                     else "tag matches" if old == taxid else "RE-TAGGED"] += 1
                if old is not None and old != taxid:
                    print(f"  re-tag {acc}: {old} → {taxid}  ({r['organism_name']})")
                if args.dry_run:
                    continue
                dest = stage / acc / "tagged.fna"
                tag_fasta(src, taxid, dest)
                err = add_to_library(dest, db_dir)
                shutil.rmtree(dest.parent, ignore_errors=True)
                if err:
                    failures.append((acc, err))
                else:
                    added += 1
                if i % 250 == 0 or i == len(rows):
                    print(f"  [{i:,}/{len(rows):,}]  added {added:,}", flush=True)

            print("\n── header tags ──────────────────────────────────────────")
            for k in ("untagged", "tag matches", "RE-TAGGED"):
                if tags.get(k):
                    print(f"  {k:12s} {tags[k]:6,d}")
            if tags.get("RE-TAGGED"):
                print("  RE-TAGGED rows carried a taxid from a previous database. The manifest "
                      "wins; this is the bug the 03_kraken tagger hid.")

            shutil.rmtree(stage, ignore_errors=True)

            if failures:
                path = SILVER / "build_failures.tsv"
                with atomic_write(path) as fh:
                    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                    w.writerow(["accession", "error"])
                    w.writerows(failures)
                sys.exit(f"\n{len(failures):,} of {len(rows):,} assemblies did not reach the "
                         f"library → {path}\nNOT building. db_v3 warned here and built anyway, "
                         "which is why nothing downstream knew it was 34 assemblies short.")

            if args.dry_run:
                print("\n--dry-run: nothing added, nothing built.")
                return
            print(f"\nlibrary complete: {added:,} of {len(rows):,}")

        download_taxonomy(db_dir)
        build(db_dir, args.threads)
        inspect(db_dir, DB_INSPECT, args.threads)

        hash_k2d = db_dir / "hash.k2d"
        if hash_k2d.exists():
            print(f"\ndb ready: {db_dir}")
            print(f"  hash.k2d {hash_k2d.stat().st_size/1e9:.1f} GB  "
                  f"(db_v3 was 27.4 GB; this is the RAM a classify job must hold)")
        print(f"\nclassify with: kraken2 --db {db_dir} --confidence 0.15 "
              "--minimum-hit-groups 3 --threads N ...")
    finally:
        log.close()


if __name__ == "__main__":
    main()
