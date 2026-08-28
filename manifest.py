#!/usr/bin/env python3
"""
manifest.py — generate a manifest.tsv for a large, gitignored data directory.

Any module that keeps large regeneratable-but-expensive data on Setonix scratch
(gitignored, too big to commit) can use this to record what exists there — a
small manifest.tsv (relative path, size, mtime, +md5 for checksum-worthy files)
that IS committed, so the repo shows what's on remote scratch without holding
the data itself. First used by kraken/ (CDS downloads, BUSCO lineage caches,
Kraken2 DBs, downloaded FASTQ, ~200GB total) — reusable by any module.

Usage:
    python manifest.py --data-dir <path> [--out <path>/manifest.tsv]
                        [--checksum-suffix .k2d [--checksum-suffix .ext2 ...]]

Run once per data directory after any structural change (new DB build, new
reads downloaded, etc.) — not on every pipeline run. Requires python/3.11
module loaded on Setonix (same as every other script that imports _util).
"""
import argparse
import sys
from pathlib import Path

from _util import build_manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True, type=Path, help="Data dir to manifest")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output manifest path (default: <data-dir>/manifest.tsv)")
    ap.add_argument("--checksum-suffix", action="append", default=None,
                    help="File suffix worth md5-checksumming (repeatable; "
                         "default: .k2d — Kraken2 DB files)")
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        sys.exit(f"Not a directory: {data_dir}")
    out_path = args.out or (data_dir / "manifest.tsv")
    suffixes = tuple(args.checksum_suffix) if args.checksum_suffix else (".k2d",)
    build_manifest(data_dir, out_path, checksum_suffixes=suffixes)


if __name__ == "__main__":
    main()
