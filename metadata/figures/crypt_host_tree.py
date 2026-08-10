#!/usr/bin/env python3
"""
crypt_host_tree.py

Build NCBI taxonomy tree for plant host species confirmed in crypt MAL mode.
Each tip = one plant species; n_confirmed = gate-passed runs that detected that host.

Run as: python3 figure/host_tree/crypt_host_tree.py   (from crypt/)
Must use system python3 (ete3 + sqlite3 incompatibility with miniconda).

Inputs:
  output/02_filter_runs/data/runs.tsv

Outputs:
  figure/host_tree/crypt_host_tree.nwk
  figure/host_tree/crypt_host_tree_meta.tsv
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from ete3 import NCBITaxa

RUNS_TSV           = Path("stat/output/data/runs.tsv")
OUT_NWK             = Path("metadata/output/figures/host_tree/crypt_host_tree.nwk")
OUT_META            = Path("metadata/output/figures/host_tree/crypt_host_tree_meta.tsv")
VIRIDIPLANTAE_TAXID = 33090

FAMILY_TAXIDS: list[tuple[int, str]] = [
    ( 4479, "Poaceae"),
    ( 3803, "Fabaceae"),
    ( 4070, "Solanaceae"),
    ( 4210, "Asteraceae"),
    ( 3745, "Rosaceae"),
    ( 3700, "Brassicaceae"),
    (23513, "Rutaceae"),
    ( 3602, "Vitaceae"),
    ( 3629, "Malvaceae"),
    ( 4637, "Musaceae"),
    ( 3650, "Cucurbitaceae"),
    ( 3977, "Euphorbiaceae"),
    ( 4668, "Amaryllidaceae"),
    ( 3623, "Actinidiaceae"),
    ( 4710, "Arecaceae"),
    ( 3647, "Caricaceae"),
]


def safe_label(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", name.strip())


def is_binomial(name: str) -> bool:
    words = name.strip().split()
    return len(words) >= 2 and words[0][0].isupper() and words[1][0].islower()


def assign_family(taxid: int, lineage_cache: dict, ncbi: NCBITaxa) -> str:
    if taxid not in lineage_cache:
        lineage_cache[taxid] = set(ncbi.get_lineage(taxid))
    lineage = lineage_cache[taxid]
    for fam_taxid, fam_name in FAMILY_TAXIDS:
        if fam_taxid in lineage:
            return fam_name
    return ""


def main() -> None:
    host_data: dict = defaultdict(lambda: {"single": 0, "multi": 0})
    with open(RUNS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("biosample_representative", "").strip() != "True":
                continue
            host = row.get("host", "").strip()
            flag = row.get("co_infection_flag", "")
            if host and is_binomial(host):
                if flag == "single":
                    host_data[host]["single"] += 1
                else:
                    host_data[host]["multi"] += 1

    host_counts = {sp: v["single"] + v["multi"] for sp, v in host_data.items()}
    print(f"Candidate host species (binomial): {len(host_counts):,}")
    print(f"Confirmed runs with binomial host:  {sum(host_counts.values()):,}")

    ncbi = NCBITaxa()
    name2taxid = ncbi.get_name_translator(list(host_counts.keys()))

    lineage_cache: dict[int, set] = {}
    resolved, skipped = [], []
    for sp, count in host_counts.items():
        if sp not in name2taxid:
            skipped.append((sp, "unresolved"))
            continue
        taxid = name2taxid[sp][0]
        lineage = set(ncbi.get_lineage(taxid))
        lineage_cache[taxid] = lineage
        if VIRIDIPLANTAE_TAXID not in lineage:
            skipped.append((sp, "not Viridiplantae"))
            continue
        resolved.append({
            "species":     sp,
            "taxid":       taxid,
            "n_single":    host_data[sp]["single"],
            "n_multi":     host_data[sp]["multi"],
            "n_confirmed": count,
        })

    print(f"Resolved + Viridiplantae: {len(resolved):,}  |  skipped: {len(skipped):,}")
    if skipped[:5]:
        print(f"  First skipped: {skipped[:5]}")

    fam_counts: Counter = Counter()
    for r in resolved:
        fam = assign_family(r["taxid"], lineage_cache, ncbi)
        r["family"] = fam
        fam_counts[fam or "(other)"] += 1

    print("Family breakdown:")
    for fam, n in sorted(fam_counts.items(), key=lambda x: -x[1]):
        print(f"  {fam:<20} {n:>5}")

    taxid2row: dict[int, dict] = {}
    for r in resolved:
        tid = r["taxid"]
        if tid not in taxid2row or r["n_confirmed"] > taxid2row[tid]["n_confirmed"]:
            taxid2row[tid] = r
        else:
            taxid2row[tid]["n_single"]    += r["n_single"]
            taxid2row[tid]["n_multi"]     += r["n_multi"]
            taxid2row[tid]["n_confirmed"] += r["n_confirmed"]

    taxids = list(taxid2row.keys())
    print(f"Unique taxids: {len(taxids):,}")

    print("Building NCBI taxonomy tree …")
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)

    for leaf in tree.get_leaves():
        tid = int(leaf.name)
        if tid in taxid2row:
            leaf.name = safe_label(taxid2row[tid]["species"])

    for node in tree.traverse():
        if not node.is_leaf():
            node.name = ""
        node.dist = 1.0

    OUT_NWK.parent.mkdir(parents=True, exist_ok=True)
    tree.write(format=1, outfile=str(OUT_NWK))
    print(f"Written: {OUT_NWK}  ({len(tree.get_leaves())} tips)")

    fields = ["label", "species", "n_single", "n_multi", "n_confirmed", "family"]
    with open(OUT_META, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in resolved:
            w.writerow({
                "label":       safe_label(r["species"]),
                "species":     r["species"],
                "n_single":    r["n_single"],
                "n_multi":     r["n_multi"],
                "n_confirmed": r["n_confirmed"],
                "family":      r["family"],
            })
    print(f"Written: {OUT_META}")


if __name__ == "__main__":
    main()
