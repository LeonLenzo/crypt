#!/usr/bin/env python3
"""
crypt_host_tree.py

NCBI taxonomy tree of plant host species across the 2,719-sample field/aerial
validation cohort (metadata/output/meta_classify/data/samples.tsv, same
selection as kraken/run/kraken_run_select.py's default scope). Each tip = one
host species; n_biosamples = BioSamples in the cohort resolved to that host.

Rebuilt 2026-09-01 from the retired metadata/legacy/figures/crypt_host_tree.py
(which used stat/stat_filter's raw, sometimes-generic STAT-inferred `host`
column, MAL mode only, and split single- vs multi-pathogen bars). Now uses
meta_classify.py's LLM-extracted + disambiguated `llm_host_resolved_taxid` —
a materially more accurate host call, already deterministically resolved to a
taxid (_util.resolve_taxon_name, not re-resolved here) — and the single-vs-
multi split is dropped: this is a host census, not a co-infection-rate figure
(see metadata/figures/sample_funnel_v3.py for that).

Must use system python3, not miniconda Python (ete3/sqlite3 incompatibility).

Run as: python3 metadata/figures/crypt_host_tree.py   (from crypt/)

Inputs:
  metadata/output/meta_classify/data/samples.tsv

Outputs:
  metadata/output/figures/host_tree/crypt_host_tree.nwk
  metadata/output/figures/host_tree/crypt_host_tree_meta.tsv
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from ete3 import NCBITaxa

SAMPLES_TSV         = Path("metadata/output/meta_classify/data/samples.tsv")
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


def ti(r: dict) -> str:
    t = r.get("llm_tissue", "")
    return t if t in ("aerial", "non-aerial") else "unclear"


def se(r: dict) -> str:
    s = r.get("llm_study_setting", "")
    return s if s in ("field", "greenhouse", "growth_chamber",
                       "detached_leaf_assay", "in_vitro") else "unclear"


def assign_family(taxid: int, lineage_cache: dict, ncbi: NCBITaxa) -> str:
    if taxid not in lineage_cache:
        lineage_cache[taxid] = set(ncbi.get_lineage(taxid))
    lineage = lineage_cache[taxid]
    for fam_taxid, fam_name in FAMILY_TAXIDS:
        if fam_taxid in lineage:
            return fam_name
    return ""


def main() -> None:
    with open(SAMPLES_TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    targets = [r for r in rows if ti(r) == "aerial" and se(r) == "field"]
    print(f"Target BioSamples (field, aerial): {len(targets):,}")

    # llm_host_resolved_taxid already covers every single-candidate BioProject
    # (meta_classify.py sets host_resolved = named_hosts[0] there directly) as
    # well as every confidently-disambiguated multi-host BioSample. What's left
    # blank is genuinely ambiguous (unresolved multi-host, no bs_host signal) —
    # excluded here rather than guessed, same convention as the rest of the
    # pipeline (kraken_run_split.py resolves those via read alignment instead).
    taxid_counts: Counter = Counter()
    taxid_to_name: dict[int, str] = {}
    n_resolved = n_unresolved = 0
    for r in targets:
        tid_s = (r.get("llm_host_resolved_taxid", "") or "").strip()
        if not tid_s:
            n_unresolved += 1
            continue
        tid = int(tid_s)
        taxid_counts[tid] += 1
        if tid not in taxid_to_name:
            taxid_to_name[tid] = r.get("llm_host_resolved", "").strip()
        n_resolved += 1

    print(f"Resolved to a host taxid: {n_resolved:,}  |  "
          f"unresolved (excluded, genuinely ambiguous): {n_unresolved:,}")
    print(f"Distinct host taxids: {len(taxid_counts):,}")

    ncbi = NCBITaxa()
    lineage_cache: dict[int, set] = {}
    resolved, skipped = [], []
    for tid, count in taxid_counts.items():
        lineage = set(ncbi.get_lineage(tid))
        lineage_cache[tid] = lineage
        if VIRIDIPLANTAE_TAXID not in lineage:
            skipped.append((taxid_to_name.get(tid, str(tid)), "not Viridiplantae"))
            continue
        resolved.append({
            "taxid":        tid,
            "species":      taxid_to_name.get(tid, str(tid)),
            "n_biosamples": count,
        })

    print(f"Viridiplantae-confirmed: {len(resolved):,}  |  skipped: {len(skipped):,}")
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

    taxids = [r["taxid"] for r in resolved]
    print(f"Unique taxids for tree: {len(taxids):,}")

    print("Building NCBI taxonomy tree …")
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)

    taxid2row = {r["taxid"]: r for r in resolved}
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

    fields = ["label", "species", "n_biosamples", "family"]
    with open(OUT_META, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in resolved:
            w.writerow({
                "label":        safe_label(r["species"]),
                "species":      r["species"],
                "n_biosamples": r["n_biosamples"],
                "family":       r["family"],
            })
    print(f"Written: {OUT_META}")


if __name__ == "__main__":
    main()
