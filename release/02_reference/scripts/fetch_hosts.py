#!/usr/bin/env python3
"""Download Ensembl Plants' species table: what host CDS exists to build the sink from.

    python release/02_reference/scripts/fetch_hosts.py
    python release/02_reference/scripts/fetch_hosts.py --force

Reads  data/inputs/scope_rule.tsv                   source=ensembl_plants, the pinned release
Writes data/raw/raw_ensembl_plants_species.txt      NOT tracked
       data/raw/raw_sources_manifest.tsv            this source's provenance rows, TRACKED

THE HOSTS ARE NOT IN THE SCOPE RULE, AND THAT IS WHY THIS SCRIPT EXISTS SEPARATELY. The
reference scope is "every fungus USDA ever recorded on a plant"; a wheat genome satisfies no
clause of it. When 02_reference was written around that sentence the host sink fell out of the
database silently, and db0/db1 were built with `Fungi 100%` and no plant CDS at all.

WHY THE DATABASE NEEDS PLANTS IN IT, which is easy to argue backwards:

1. THE SINK. Without host CDS a plant read is absorbed by the nearest fungal CDS and
   reported as a fungal hit. These are plant RNA-seq libraries, so most reads are wheat, and
   that is a false-positive mechanism aimed straight at the chapter's claim.
2. THE HOST CALL IS A RESULT. 04_assign ranks reads across hosts to decide which plant a
   library came from, and agrees with the LLM-curated host 91.2% of the time. With no host
   CDS there is no host call to make.

Neither is contradicted by the measured finding that HOST-READ REMOVAL is unnecessary (PST
read counts ratio 1.00x raw versus HISAT2-split). That tested deleting host reads BEFORE
classification, which changes the percentage and not what Kraken2 finds. Putting host CDS IN
the database is the opposite operation: it gives a plant read a correct home rather than a
fungal one.

## ENSEMBL PLANTS, NOT NCBI, AND THE DIFFERENCE IS DECISIVE

NCBI GenBank holds the wheat pangenome as ASSEMBLIES and not as gene models. Of its twelve
*T. aestivum* entries carrying an `annotation_provider`, **eight have
protein_coding_gene_count = 0** -- the NIAB cultivars (Alchemy, Xi19, Robigus, Soissons,
Brompton, Rialto, Claire) are assembled and not annotated -- leaving only Chinese Spring,
annotated twice. Taking "all available wheat CDS" from NCBI is therefore two annotations of
one cultivar plus three fragments.

**Ensembl Plants publishes 19 *T. aestivum* cultivars WITH CDS**, which is the breadth a sink
actually wants, and it is a better provenance story as well: Ensembl ships immutable numbered
releases, so release-116 is a stable artefact a reviewer can re-fetch byte-identically. The
fungal `assembly_summary.txt` is rebuilt daily and cannot be pinned at all.

This script fetches only the SPECIES TABLE. `screen_hosts.py` decides which hosts and which
cultivars, and `fetch_host_cds.py` pulls the FASTA, mirroring
fetch_assemblies.py -> screen_assemblies.py -> fetch_cds.py on the fungal side.
"""

import argparse
import csv
import hashlib
import shutil
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
RULE = HERE / "data/inputs/scope_rule.tsv"
RAW = HERE / "data/raw/raw_ensembl_plants_species.txt"
MANIFEST = HERE / "data/raw/raw_sources_manifest.tsv"
SOURCE = "ensembl_plants"
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
        last_mod = resp.headers.get("Last-Modified", "")
        shutil.copyfileobj(resp, fh, length=1 << 20)
    tmp.replace(dest)
    return last_mod


def rows_of(path):
    """The table's header line starts with '#name', so it is read positionally."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        hdr = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [c.lstrip("#") for c in f]
                continue
            if line.strip():
                yield dict(zip(hdr, f))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true",
                    help="re-download even if a local copy is present")
    args = ap.parse_args()

    rule = read_rule()
    pin = rule[SOURCE]
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
            raise SystemExit(f"md5 {have} != pinned {want}; file removed. Release "
                             f"{pin['release']} is immutable, so a mismatch means the pin "
                             f"names a different release than the url.")
        print(f"  {RAW.stat().st_size:,} bytes, md5 {have}")

    rows = list(rows_of(RAW))
    by_tid = Counter(r["taxonomy_id"] for r in rows)
    print(f"  Ensembl Plants release {pin['release']}: {len(rows)} genomes, "
          f"{len(by_tid)} distinct taxids")

    # Report what is on offer for the hosts the rule names, as a reachability check. A host
    # declared in the rule with no Ensembl genome has no CDS and screen_hosts.py would
    # otherwise drop it without saying why.
    host = rule["host"]
    for label, tid in (("analysed", host.get("analysed")), ("outgroup", host.get("outgroup"))):
        if not tid:
            continue
        hits = [r for r in rows if r["taxonomy_id"] == tid]
        print(f"    {label:<9} taxid {tid:<8} {len(hits)} Ensembl genomes")
        if not hits:
            print(f"      note: no genome under this taxid. Ensembl may file it under a "
                  f"SUBSPECIES node (it files barley under 112509, not 4513), which "
                  f"screen_hosts.py resolves by walking the taxonomy.")

    if not want:
        print(f"\n  UNPINNED. Release {pin['release']} is immutable, so freezing this table "
              f"is cheap: set expected_md5 in inputs/scope_rule.tsv to\n    {have}")

    write = [("source_name", pin["source_name"]), ("licence", pin["licence"]),
             ("release", pin["release"]), ("url", pin["url"]),
             ("file", pin["expected_file"]), ("bytes", RAW.stat().st_size), ("md5", have),
             ("last_modified", last_mod), ("pinned", "yes" if want else "no"),
             ("retrieved_on", date.today().isoformat()), ("n_genomes", len(rows)),
             ("n_taxids", len(by_tid))]
    keep = []
    if MANIFEST.exists():
        with open(MANIFEST) as fh:
            keep = [r for r in csv.DictReader(fh, delimiter="\t") if r["source"] != SOURCE]
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source", "field", "value"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(keep)
        for k, v in write:
            w.writerow({"source": SOURCE, "field": k, "value": v})
    print(f"\n  -> {MANIFEST}")
    print("Next: screen_hosts.py, choosing the hosts and their cultivars")


if __name__ == "__main__":
    main()
