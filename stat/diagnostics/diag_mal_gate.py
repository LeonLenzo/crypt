#!/usr/bin/env python3
"""
Diagnostic 1: Absolute read count check on low-euk MAL group.

Tests whether fungal signal in low-euk MAL runs (euk_pct < 1%) scales with
sequencing depth (real signal) or is constant in absolute reads (noise floor).

  constant count → noise / background contamination
  constant pct   → biological signal, just low resolution

Run from crypt/: python scripts/diag_mal_gate.py
"""
import csv
import math
import statistics
from pathlib import Path

RUNS_TSV = Path("output/02_filter_runs/data/runs.tsv")

# ── Load MAL biosample_rep runs ───────────────────────────────────────────────

rows = []
with open(RUNS_TSV) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["mode"] != "mal" or r["biosample_representative"] != "True":
            continue
        analyzed = float(r["analyzed"] or 0)
        if analyzed == 0:
            continue
        fungi     = float(r["fungi_pct"]    or 0)
        oomycete  = float(r["oomycete_pct"] or 0)
        nematode  = float(r["nematode_pct"] or 0)
        euk_pct   = fungi + oomycete + nematode
        euk_count = euk_pct * analyzed / 100
        rows.append({
            "analyzed":  analyzed,
            "euk_pct":   euk_pct,
            "euk_count": euk_count,
            "fungi_pct": fungi,
            "library_organism": r["library_organism"],
            "library_detected": r["library_detected"],
        })

# ── Split into tiers ──────────────────────────────────────────────────────────

zero = [r for r in rows if r["euk_pct"] == 0]
low  = [r for r in rows if 0 < r["euk_pct"] < 1]
high = [r for r in rows if r["euk_pct"] >= 1]

print(f"MAL biosample_rep: {len(rows):,} total")
print(f"  zero euk (pct=0):  {len(zero):,}  ({100*len(zero)/len(rows):.1f}%)")
print(f"  low  euk (0<x<1):  {len(low):,}  ({100*len(low)/len(rows):.1f}%)  ← Dikarya-stall candidates")
print(f"  high euk (>=1%):   {len(high):,}  ({100*len(high)/len(rows):.1f}%)")

# ── Diagnostic: abs count vs depth for low-euk group ─────────────────────────

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n)
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n)
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (n * sx * sy)

def percentile(vals, p):
    vals_s = sorted(vals)
    idx = (p / 100) * (len(vals_s) - 1)
    lo, hi = int(idx), min(int(idx) + 1, len(vals_s) - 1)
    return vals_s[lo] + (idx - lo) * (vals_s[hi] - vals_s[lo])

print()
print("── LOW-EUK GROUP (0 < euk_pct < 1%) ─────────────────────────────────")

analyzed_low  = [r["analyzed"]  for r in low]
euk_cnt_low   = [r["euk_count"] for r in low]
euk_pct_low   = [r["euk_pct"]  for r in low]

r_cnt_depth = pearson(analyzed_low, euk_cnt_low)
r_pct_depth = pearson(analyzed_low, euk_pct_low)

print(f"  Pearson r (abs_count vs depth):  {r_cnt_depth:+.3f}")
print(f"  Pearson r (euk_pct  vs depth):   {r_pct_depth:+.3f}")
print()
print("  Interpretation:")
if r_cnt_depth > 0.8:
    print("    → abs_count SCALES with depth → signal is proportional (real biology)")
elif r_cnt_depth < 0.3 and r_pct_depth < 0:
    print("    → abs_count roughly CONSTANT; pct negatively correlated with depth → NOISE FLOOR")
else:
    print(f"    → intermediate signal (r_cnt={r_cnt_depth:.3f}) — mixed or weak correlation")

print()
print("  Absolute euk count distribution (low-euk group):")
for p in [5, 25, 50, 75, 95]:
    print(f"    p{p:2d}: {percentile(euk_cnt_low, p):>10,.0f} reads")

print()
print("  Analyzed depth distribution (low-euk group):")
for p in [5, 25, 50, 75, 95]:
    print(f"    p{p:2d}: {percentile(analyzed_low, p):>10,.0f} reads")

# ── Compare low vs high abs counts ───────────────────────────────────────────

print()
print("── COMPARISON: median absolute euk counts ────────────────────────────")
euk_cnt_high = [r["euk_count"] for r in high]
print(f"  LOW  group median abs count:  {statistics.median(euk_cnt_low):>12,.0f} reads")
print(f"  HIGH group median abs count:  {statistics.median(euk_cnt_high):>12,.0f} reads")
print(f"  Ratio HIGH/LOW:               {statistics.median(euk_cnt_high)/statistics.median(euk_cnt_low):>12.1f}×")

# ── Depth-binned analysis: is pct stable across depth bins? ──────────────────

print()
print("── DEPTH-BINNED euk_pct IN LOW-EUK GROUP (noise test) ──────────────")
bins = [(0, 5e6), (5e6, 20e6), (20e6, 60e6), (60e6, float("inf"))]
labels = ["<5M", "5-20M", "20-60M", ">60M"]
for (lo_d, hi_d), label in zip(bins, labels):
    grp = [r for r in low if lo_d <= r["analyzed"] < hi_d]
    if not grp:
        continue
    pcts = [r["euk_pct"] for r in grp]
    cnts = [r["euk_count"] for r in grp]
    print(f"  {label:>8}  n={len(grp):>4}  "
          f"median_pct={statistics.median(pcts):.4f}%  "
          f"median_abs={statistics.median(cnts):>9,.0f} reads")

print()
print("  → If median_pct FALLS as depth increases: noise (constant count)")
print("  → If median_pct STABLE across bins: real signal (constant pct)")

# ── Library detection check ───────────────────────────────────────────────────

print()
print("── LIBRARY_DETECTED IN LOW-EUK GROUP ───────────────────────────────")
det     = sum(1 for r in low if r["library_detected"] == "True")
not_det = sum(1 for r in low if r["library_detected"] != "True")
print(f"  library_detected=True:   {det:>4}  ({100*det/len(low):.1f}%)")
print(f"  library_detected=False:  {not_det:>4}  ({100*not_det/len(low):.1f}%)")
print()
print("  (library_detected=False + no euk species = submitted as pathogen but")
print("   pathogen absent from library — host-only experiment)")
