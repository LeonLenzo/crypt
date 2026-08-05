#!/usr/bin/env python3
"""
Select runs for the kallisto MAL gate pilot.

Three organisms × three tiers:
  zero   — euk_pct = 0 (zero-kmer, negative control)
  low    — 0 < euk_pct < 1% (Dikarya-stall, main test group)
  high   — euk_pct >= 1%, library_detected=True (positive controls)

Stratified by BioProject so runs come from different studies.

Run from crypt/: python kraken/pilot/kallisto/select_runs.py
Output: kraken/pilot/kallisto/pilot_runs.tsv
"""
import csv
import random
from collections import defaultdict
from pathlib import Path

RUNS_TSV = Path("output/02_filter_runs/data/runs.tsv")
OUT_TSV  = Path("kraken/pilot/kallisto/pilot_runs.tsv")

ORGANISMS = {
    "pst": "Puccinia striiformis f. sp. tritici",
    "pgt": "Puccinia graminis f. sp. tritici",
    "por": "Pyricularia oryzae",
}

N_PER_TIER = 5
RANDOM_SEED = 42

random.seed(RANDOM_SEED)

# ── Load MAL biosample_rep runs ───────────────────────────────────────────────

# org_key → tier → list of (BioProject, Run, euk_pct, analyzed, library_detected)
buckets = {k: {"zero": [], "low": [], "high": []} for k in ORGANISMS}

with open(RUNS_TSV) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["mode"] != "mal" or r["biosample_representative"] != "True":
            continue
        lib_org = r["library_organism"].strip()
        org_key = None
        for k, name in ORGANISMS.items():
            if lib_org == name:
                org_key = k
                break
        if org_key is None:
            continue
        euk_pct = (float(r["fungi_pct"] or 0)
                   + float(r["oomycete_pct"] or 0)
                   + float(r["nematode_pct"] or 0))
        row = (r["BioProject"], r["Run"], round(euk_pct, 4),
               int(float(r["analyzed"])), r["library_detected"])
        if euk_pct == 0:
            buckets[org_key]["zero"].append(row)
        elif euk_pct < 1.0:
            buckets[org_key]["low"].append(row)
        elif r["library_detected"] == "True":
            buckets[org_key]["high"].append(row)

# ── Stratified sampling: prefer different BioProjects ─────────────────────────

def stratified_sample(rows, n):
    """Pick n rows, preferring one per BioProject."""
    by_bp = defaultdict(list)
    for row in rows:
        by_bp[row[0]].append(row)
    selected = []
    # round-robin across BioProjects
    pools = list(by_bp.values())
    random.shuffle(pools)
    for pool in pools:
        random.shuffle(pool)
    i = 0
    while len(selected) < n and i < max(len(p) for p in pools):
        for pool in pools:
            if i < len(pool) and len(selected) < n:
                selected.append(pool[i])
        i += 1
    return selected

# ── Write output ──────────────────────────────────────────────────────────────

print(f"{'Organism':<40} {'Tier':<6} {'n_avail':>8} {'n_sel':>6}")
print("-" * 65)

selected_rows = []
for org_key, tiers in buckets.items():
    org_name = ORGANISMS[org_key]
    for tier, rows in tiers.items():
        sel = stratified_sample(rows, N_PER_TIER)
        print(f"{org_name:<40} {tier:<6} {len(rows):>8} {len(sel):>6}")
        for bp, run, euk_pct, analyzed, lib_det in sel:
            selected_rows.append({
                "org_key": org_key,
                "organism": org_name,
                "tier": tier,
                "BioProject": bp,
                "Run": run,
                "euk_pct": euk_pct,
                "analyzed": analyzed,
                "library_detected": lib_det,
            })

print()
print(f"Total runs selected: {len(selected_rows)}")

with open(OUT_TSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["org_key", "organism", "tier", "BioProject",
                                       "Run", "euk_pct", "analyzed", "library_detected"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(selected_rows)

print(f"Written: {OUT_TSV}")
