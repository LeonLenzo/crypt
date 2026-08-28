#!/usr/bin/env python3
"""kraken_manifest.py — generate a manifest.tsv for a large, gitignored data directory.

Setonix scratch holds ~200GB of kraken pipeline data (CDS downloads, BUSCO lineage
caches, Kraken2 DBs, downloaded FASTQ reads) that lives under crypt/kraken/output/*/data/
but is too large to commit. This script walks a data directory and writes a small
manifest.tsv (relative path, size in bytes, mtime) that IS committed, so the repo shows
what exists on Setonix scratch without needing the data itself.

For Kraken2 DB files (*.k2d), also records an md5 checksum — the one artifact where
byte-level reproducibility actually matters (confirms a DB wasn't silently rebuilt
differently); skipped for everything else since checksumming ~200GB of FASTQ/CDS would
take too long for no real benefit.

Usage:
    python kraken/kraken_manifest.py --data-dir <path> [--out <path>/manifest.tsv]

Run once per data directory after any structural change (new DB build, new reads
downloaded, etc.) — not on every pipeline run.
"""
import argparse
import hashlib
import sys
from pathlib import Path


def _md5(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def build_manifest(data_dir: Path, out_path: Path) -> None:
    rows = []
    n_files = 0
    total_bytes = 0
    for p in sorted(data_dir.rglob("*")):
        if p.name == "manifest.tsv" or not p.is_file():
            continue
        rel = p.relative_to(data_dir)
        st = p.stat()
        checksum = _md5(p) if p.suffix == ".k2d" else ""
        rows.append((str(rel), st.st_size, int(st.st_mtime), checksum))
        n_files += 1
        total_bytes += st.st_size

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("path\tsize_bytes\tmtime_epoch\tmd5\n")
        for rel, size, mtime, checksum in rows:
            f.write(f"{rel}\t{size}\t{mtime}\t{checksum}\n")

    print(f"{out_path}: {n_files} files, {total_bytes / 1e9:.2f} GB")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path, help="Data dir to manifest")
    ap.add_argument("--out", type=Path, default=None,
                     help="Output manifest path (default: <data-dir>/manifest.tsv)")
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        sys.exit(f"Not a directory: {data_dir}")
    out_path = args.out or (data_dir / "manifest.tsv")
    build_manifest(data_dir, out_path)


if __name__ == "__main__":
    main()
