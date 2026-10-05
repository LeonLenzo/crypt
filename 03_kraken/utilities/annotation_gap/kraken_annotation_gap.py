#!/usr/bin/env python3
"""kraken_annotation_gap.py — how many assemblies in our pathogen genera have no CDS, and why.

NCBI's `cds_from_genomic.fna` is an extraction of annotation already in the record, not
something NCBI computes. GenBank (GCA) carries whatever the submitter deposited and annotation
is optional; RefSeq (GCF) runs NCBI's own eukaryotic pipeline selectively. So assembly quality
and annotation availability are independent, and a chromosome-level assembly can have no CDS
(Puccinia polysora GCA_025617555.3, N50 54 Mb, unannotated).

Because kraken_search.py requires scaffold_plus AND has_annotation, every unannotated assembly
was ineligible for the database regardless of quality. This measures the size of that exclusion
per genus, which bounds what gene prediction could add.

Counting trap this script exists to avoid: taking the first word of `organism_name` and querying
NCBI per name inflates the total badly. Synonyms resolve to the same clade (`Fusarium` and
`[Neocosmospora]` both return 2,238) and the reference table holds some order-level names
(`Pleosporales`, `Xylariales`) whose counts overlap dozens of real genera. Each name is resolved
through the taxonomy API, only rank GENUS is kept, and tax_ids are deduplicated.

A utility, not a pipeline step: it informs database sizing, it does not shape a result.

    python 03_kraken/utilities/annotation_gap/kraken_annotation_gap.py
"""
import collections
import csv
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from _paths import REF_CANDIDATES  # noqa: E402

OUT = Path(__file__).resolve().parent / "data/annotation_gap.tsv"
TAXON = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/"
GENOME = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/"
PAUSE = 0.34   # NCBI allows ~3 requests/second unauthenticated


def _get(url, tries=3):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "crypt/annotation_gap"})
            with urllib.request.urlopen(req, timeout=45) as h:
                return json.load(h)
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def resolve(name):
    """-> (taxid, RANK) for a taxon name, or (None, None)."""
    d = _get(TAXON + urllib.request.quote(name))
    for node in (d or {}).get("taxonomy_nodes", []):
        t = node.get("taxonomy") or {}
        if t.get("tax_id"):
            return t["tax_id"], (t.get("rank") or "").upper()
    return None, None


def count(taxid, annotated=None):
    q = f"{taxid}/dataset_report?page_size=1&filters.assembly_version=current"
    if annotated is not None:
        q += f"&filters.has_annotation={'true' if annotated else 'false'}"
    d = _get(GENOME + q)
    return None if d is None else (d.get("total_count") or 0)


def main():
    cand, spp = collections.Counter(), collections.defaultdict(set)
    for r in csv.DictReader(open(REF_CANDIDATES), delimiter="\t"):
        parts = r["organism_name"].split()
        if len(parts) >= 2:
            g = parts[0].strip("[]")
            cand[g] += 1
            spp[g].add(" ".join(parts[:2]))

    print(f"{len(cand)} candidate genus names; resolving ranks", flush=True)
    rows, by_taxid, dropped = [], {}, collections.Counter()
    for i, (g, nacc) in enumerate(sorted(cand.items()), 1):
        taxid, rank = resolve(g)
        time.sleep(PAUSE)
        if taxid is None:
            dropped["unresolved"] += 1
            continue
        if rank != "GENUS":
            dropped[f"not a genus ({rank or 'unknown'})"] += 1
            continue
        if taxid in by_taxid:                      # synonym of a genus already counted
            dropped["synonym"] += 1
            by_taxid[taxid]["db_accessions"] += nacc
            continue
        tot, ann = count(taxid), count(taxid, True)
        time.sleep(PAUSE)
        if tot is None or ann is None:
            dropped["query failed"] += 1
            continue
        row = {"genus": g, "taxid": taxid, "db_accessions": nacc, "db_species": len(spp[g]),
               "ncbi_total": tot, "ncbi_annotated": ann, "ncbi_unannotated": tot - ann}
        rows.append(row)
        by_taxid[taxid] = row
        if i % 40 == 0:
            print(f"  {i}/{len(cand)}", flush=True)

    rows.sort(key=lambda r: -r["ncbi_unannotated"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    with open(OUT, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    total = lambda k: sum(r[k] for r in rows)
    print(f"\n{len(rows)} true genera (dropped: {dict(dropped)})")
    print(f"  accessions selected into the database : {total('db_accessions'):>7,}")
    print(f"  assemblies at NCBI in those genera    : {total('ncbi_total'):>7,}")
    print(f"    annotated, CDS available            : {total('ncbi_annotated'):>7,}")
    print(f"    unannotated, the prediction gain    : {total('ncbi_unannotated'):>7,}")
    print(f"\n  written: {OUT}")


if __name__ == "__main__":
    main()
