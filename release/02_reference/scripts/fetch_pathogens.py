#!/usr/bin/env python3
"""Download PHI-base, the source of the pathogen label.

    python release/02_reference/scripts/fetch_pathogens.py
    python release/02_reference/scripts/fetch_pathogens.py --force

Reads  data/inputs/scope_rule.tsv        the phibase pin (url, expected file, md5)
Writes data/raw/raw_phibase.csv          PHI-base v4-19, 17 MB, NOT tracked
       data/raw/raw_sources_manifest.tsv this source's provenance rows, TRACKED

Its own step, paired with screen_pathogens.py, because the pathogen label is a separate
concern from the reference scope: it never enters the Kraken2 database and is joined to a
detection afterwards. fetch_species.py fetches the scope sources (USDA, NCBI Taxonomy);
fetch_assemblies.py fetches what is buildable; this fetches what is pathogenic.

PHI-base publishes dated, immutable releases under releases/, so a specific version is pinned
by md5 (not releases/phi-base_current.csv, which rolls forward to the newest version). A
mismatch is a hard stop: re-pin inputs/scope_rule.tsv deliberately after looking at the diff.
"""

import argparse
import csv
import hashlib
import shutil
import urllib.request
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
RULE = HERE / "data/inputs/scope_rule.tsv"
RAW = HERE / "data/raw/raw_phibase.csv"
MANIFEST = HERE / "data/raw/raw_sources_manifest.tsv"
SOURCE = "phibase"
UA = {"User-Agent": "crypt/02_reference (leon.lenzo@curtin.edu.au)"}


def read_rule():
    out = {}
    with open(RULE) as fh:
        for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                  delimiter="\t"):
            out.setdefault(row["source"], {})[row["field"]] = row["value"]
    return out


def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    tmp.replace(dest)


def ensure(dest, url, want_md5, force):
    if dest.exists() and not force:
        have = md5_of(dest)
        if have == want_md5:
            print(f"  {dest.name} already present and verified")
            return dest.stat().st_size, have
        raise SystemExit(f"{dest} md5 is {have}, not the pinned {want_md5}. "
                         f"Delete it or pass --force.")
    print(f"  downloading -> {dest.name}")
    download(url, dest)
    have = md5_of(dest)
    if have != want_md5:
        dest.unlink()
        raise SystemExit(f"{dest.name} md5 {have} != pinned {want_md5}; file removed")
    print(f"  verified {dest.stat().st_size:,} bytes")
    return dest.stat().st_size, have


def census(path):
    """Rows, and distinct fungal/plant pathogen species is a job for screen_pathogens; here
    just the row count so the manifest records what was fetched."""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in csv.reader(fh)) - 1


def write_manifest(rows):
    keep = []
    if MANIFEST.exists():
        with open(MANIFEST) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["source"] != SOURCE:
                    keep.append(r)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source", "field", "value"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(keep)
        for k, v in rows:
            w.writerow({"source": SOURCE, "field": k, "value": v})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the local copy already verifies")
    args = ap.parse_args()

    pin = read_rule()[SOURCE]
    print(f"PHI-base: {pin['expected_file']}")
    size, md5 = ensure(RAW, pin["url"], pin["expected_md5"], args.force)
    n = census(RAW)
    print(f"  {n:,} interaction rows")
    write_manifest([("source_name", pin["source_name"]), ("licence", pin["licence"]),
                    ("url", pin["url"]), ("file", pin["expected_file"]),
                    ("bytes", size), ("md5", md5),
                    ("retrieved_on", date.today().isoformat()), ("n_rows", n)])
    print(f"  -> {MANIFEST}")
    print("Next: python release/02_reference/scripts/screen_pathogens.py")


if __name__ == "__main__":
    main()
