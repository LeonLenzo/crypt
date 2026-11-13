#!/usr/bin/env python3
"""Download the two public datasets the reference scope is built from.

    python release/02_reference/scripts/fetch_species.py
    python release/02_reference/scripts/fetch_species.py --force   # re-download

Reads  data/inputs/scope_rule.tsv      the pins: figshare article, URLs, expected md5s
Writes data/raw/raw_usda_fungus_host.csv   USDA Fungus-Host, 187 MB, NOT tracked
       data/raw/raw_taxdmp.zip             NCBI Taxonomy, 80 MB, NOT tracked
       data/raw/raw_sources_manifest.tsv   what was fetched and its checksum, TRACKED

The manifest is the tracked half, so a reviewer can verify their own copies against what
this pipeline ran on without the 267 MB entering history.

Both md5s are pinned in inputs/ and every fetch is verified against the pin. USDA is checked
twice, against figshare's `computed_md5` before downloading and against the bytes after.
NCBI publishes no checksum for its dated archive, so only the post-download check applies
there; the dated file is immutable, which is why the DATED url is pinned and not
`pub/taxonomy/taxdmp.zip`, which is rebuilt nightly.
"""

import argparse
import csv
import hashlib
import json
import shutil
import urllib.request
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
RULE = HERE / "data/inputs/scope_rule.tsv"
RAW = HERE / "data/raw"
USDA = RAW / "raw_usda_fungus_host.csv"
TAXDMP = RAW / "raw_taxdmp.zip"
MANIFEST = RAW / "raw_sources_manifest.tsv"
FIGSHARE = "https://api.figshare.com/v2/articles/{}"
UA = {"User-Agent": "crypt/02_reference (leon.lenzo@curtin.edu.au)"}


def read_rule():
    """{source: {field: value}} from inputs/scope_rule.tsv."""
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
    with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    tmp.replace(dest)


def ensure(dest, url, want_md5, force):
    """Download unless a verified copy is already here. Returns (bytes, md5)."""
    if dest.exists() and not force:
        have = md5_of(dest)
        if have == want_md5:
            print(f"  {dest.name} already present and verified")
            return dest.stat().st_size, have
        raise SystemExit(f"{dest} is present but its md5 is {have}, not {want_md5}. "
                         f"Delete it or pass --force.")
    print(f"  downloading -> {dest.name}")
    download(url, dest)
    have = md5_of(dest)
    if have != want_md5:
        dest.unlink()
        raise SystemExit(f"{dest.name} md5 {have} != pinned {want_md5}; file removed")
    print(f"  verified {dest.stat().st_size:,} bytes")
    return dest.stat().st_size, have


def fetch_usda(pin, force):
    req = urllib.request.Request(FIGSHARE.format(pin["figshare_article"]), headers=UA)
    with urllib.request.urlopen(req, timeout=60) as resp:
        meta = json.load(resp)
    if len(meta["files"]) != 1:
        raise SystemExit(f"expected 1 file on figshare {pin['figshare_article']}, found "
                         f"{[f['name'] for f in meta['files']]}")
    f = meta["files"][0]
    print(f"USDA: {meta['title']}")
    print(f"  {f['name']}  {f['size']:,} bytes  md5 {f['computed_md5']}")
    if f["name"] != pin["expected_file"]:
        raise SystemExit(f"file name changed: pinned {pin['expected_file']!r}, figshare "
                         f"serves {f['name']!r}. Re-pin inputs/scope_rule.tsv.")
    if f["computed_md5"] != pin["expected_md5"]:
        raise SystemExit(f"md5 changed: pinned {pin['expected_md5']}, figshare reports "
                         f"{f['computed_md5']}. The dataset was revised; look at what "
                         f"changed before re-pinning.")
    size, md5 = ensure(USDA, f["download_url"], pin["expected_md5"], force)
    return {"source_name": meta["title"], "source_doi": meta["doi"],
            "licence": meta["license"]["name"], "file": f["name"],
            "bytes": size, "md5": md5}


def fetch_taxonomy(pin, force):
    print(f"NCBI Taxonomy: {pin['expected_file']}")
    size, md5 = ensure(TAXDMP, pin["url"], pin["expected_md5"], force)
    return {"source_name": pin["source_name"], "source_doi": "", "licence": pin["licence"],
            "file": pin["expected_file"], "bytes": size, "md5": md5}


def census_usda():
    """Records and distinct names. csv.reader rather than a line count: the `citation`
    field holds quoted text and cannot be assumed newline-free."""
    n, fungi, hosts = 0, set(), set()
    with open(USDA, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            n += 1
            if row.get("sciName"):
                fungi.add(" ".join(row["sciName"].split()))
            if row.get("host"):
                hosts.add(" ".join(row["host"].split()))
    return {"n_records": n, "n_fungus_names": len(fungi), "n_host_names": len(hosts)}


def write_manifest(rows):
    """Long format (source, field, value) so each fetch script owns its own rows and the
    sources need not share a schema. Rows for these sources are replaced, others kept."""
    mine = set(rows)
    keep = []
    if MANIFEST.exists():
        with open(MANIFEST) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["source"] not in mine:
                    keep.append(r)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source", "field", "value"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(keep)
        for src, d in rows.items():
            for k, v in [("retrieved_on", date.today().isoformat())] + list(d.items()):
                w.writerow({"source": src, "field": k, "value": v})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the local copy already verifies")
    args = ap.parse_args()

    rule = read_rule()
    rows = {"usda": fetch_usda(rule["usda"], args.force),
            "ncbi_taxonomy": fetch_taxonomy(rule["ncbi_taxonomy"], args.force)}
    rows["usda"].update(census_usda())
    print(f"  {rows['usda']['n_records']:,} records, "
          f"{rows['usda']['n_fungus_names']:,} distinct fungus names, "
          f"{rows['usda']['n_host_names']:,} distinct host names")

    write_manifest(rows)
    print(f"\n  -> {MANIFEST}")
    print("Next: python release/02_reference/scripts/screen_species.py")


if __name__ == "__main__":
    main()
