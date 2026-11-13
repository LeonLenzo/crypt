#!/usr/bin/env python3
"""Download NCBI's GenBank fungal assembly summary: what exists to build the database from.

    python release/02_reference/scripts/fetch_assemblies.py
    python release/02_reference/scripts/fetch_assemblies.py --force

Reads  data/inputs/scope_rule.tsv           the url, and an md5 pin if one has been set
Writes data/raw/raw_assembly_summary.txt    26,826 rows, 10 MB, NOT tracked
       data/raw/raw_sources_manifest.tsv    this source's provenance rows, TRACKED

One file instead of 4,138 `datasets summary genome taxon <genus>` calls, and it answers the
question holobase could not: `annotation_provider`, `annotation_date` and
`protein_coding_gene_count` say whether an assembly carries a gene model, which decides
whether it can contribute CDS at all.

GenBank only. RefSeq holds re-annotated copies of assemblies already in this file, and
taking both is how db_v4 ended up building 432 genomes twice under paired GCA/GCF
accessions. `gbrs_paired_asm` names the RefSeq twin if one is ever wanted.

UNPINNED BY DESIGN. USDA is frozen and the taxdump has dated archives, but this file is
rebuilt daily and NCBI publishes no archive of it, so there is no stable artefact to pin to.
The md5 and the server's `Last-Modified` are recorded in the manifest so the snapshot the
pipeline ran on is identifiable after the fact. Setting `expected_md5` in scope_rule.tsv
turns the pin back on and a mismatch then becomes a hard stop.
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
RAW = HERE / "data/raw/raw_assembly_summary.txt"
MANIFEST = HERE / "data/raw/raw_sources_manifest.tsv"
SOURCE = "ncbi_assembly"
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
    """Returns the server's Last-Modified, which is the only version label this file has."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as fh:
        last_mod = resp.headers.get("Last-Modified", "")
        shutil.copyfileobj(resp, fh, length=1 << 20)
    tmp.replace(dest)
    return last_mod


def census(path):
    """Rows, and how many carry a gene model. `na` is NCBI's own null in this file."""
    n = annotated = 0
    levels = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        hdr = None
        for line in fh:
            if line.startswith("#assembly_accession"):
                hdr = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or hdr is None:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            n += 1
            levels[row.get("assembly_level", "")] = \
                levels.get(row.get("assembly_level", ""), 0) + 1
            if row.get("annotation_provider", "").strip() not in ("", "na"):
                annotated += 1
    return n, annotated, levels


def write_manifest(rows):
    """Long format (source, field, value) so each fetch script owns its own rows and the
    three sources need not share a schema. This source's rows are replaced, others kept."""
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
                    help="re-download even if a local copy is present")
    args = ap.parse_args()

    pin = read_rule()[SOURCE]
    want = (pin.get("expected_md5") or "").strip()
    last_mod = ""

    if RAW.exists() and not args.force:
        have = md5_of(RAW)
        if want and have != want:
            raise SystemExit(f"{RAW} md5 is {have}, not the pinned {want}. Delete it, "
                             f"pass --force, or re-pin inputs/scope_rule.tsv.")
        print(f"  {RAW.name} already present (md5 {have}); pass --force to refresh")
    else:
        print(f"  downloading {pin['url']}")
        last_mod = download(pin["url"], RAW)
        have = md5_of(RAW)
        if want and have != want:
            RAW.unlink()
            raise SystemExit(f"md5 {have} != pinned {want}; file removed. NCBI rebuilds "
                             f"this file daily, so a stale pin is the likely cause.")
        print(f"  {RAW.stat().st_size:,} bytes, md5 {have}")
        print(f"  Last-Modified: {last_mod or '(not reported)'}")

    n, annotated, levels = census(RAW)
    print(f"  {n:,} assemblies, {annotated:,} with a gene model "
          f"({100 * annotated / n:.1f}%)")
    print("  " + ", ".join(f"{k} {v:,}" for k, v in
                           sorted(levels.items(), key=lambda kv: -kv[1])))
    if not want:
        print(f"\n  UNPINNED. To freeze this snapshot, set expected_md5 in "
              f"inputs/scope_rule.tsv to:\n    {have}")

    write_manifest([("source_name", pin["source_name"]), ("licence", pin["licence"]),
                    ("url", pin["url"]), ("file", pin["expected_file"]),
                    ("bytes", RAW.stat().st_size), ("md5", have),
                    ("last_modified", last_mod), ("pinned", "yes" if want else "no"),
                    ("retrieved_on", date.today().isoformat()),
                    ("n_assemblies", n), ("n_annotated", annotated)])
    print(f"\n  -> {MANIFEST}")
    print("Next: screen_assemblies.py, joining these against the in-scope genera")


if __name__ == "__main__":
    main()
