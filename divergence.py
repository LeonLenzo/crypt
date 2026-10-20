#!/usr/bin/env python3
"""Divergence calibration for the rusts: is within-species separable from between-species?

Pairwise CDS divergence among Puccinia assemblies, from the cached FracMinHash sketches
(k=31, scaled=1000). Jaccard -> Mash distance -> approximate ANI. If within-species pairs
cluster tight and high-ANI while between-species pairs sit clearly lower, an alignment/allelic
method can separate them; if the distributions overlap, it cannot and we stop.
"""
import csv, collections, itertools, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(".")
SK = ROOT / "superseded/03_kraken/utilities/saturation/data/sketches"
REF = ROOT / "superseded/03_kraken/01_search/data/ref_candidates_v3.tsv"
K = 31

def sp2(n):
    p = n.split(); return " ".join(p[:2]) if len(p) >= 2 else n

# Puccinia accessions with a sketch
acc = collections.defaultdict(list)
for r in csv.DictReader(open(REF), delimiter="\t"):
    s = sp2(r["organism_name"])
    if s.startswith("Puccinia") and (SK / f"{r['accession']}.npy").exists():
        acc[s].append(r["accession"])
print("Puccinia assemblies with sketches:")
for s in sorted(acc): print(f"  {s:28s} {len(acc[s])}")

sketches = {a: set(np.load(SK / f"{a}.npy").tolist())
            for alist in acc.values() for a in alist}

def ani(a, b):
    A, B = sketches[a], sketches[b]
    if not A or not B: return None
    j = len(A & B) / len(A | B)
    if j == 0: return 0.0
    d = -math.log(2 * j / (1 + j)) / K      # Mash distance
    return (1 - d) * 100

# within-species and between-species ANI
within = collections.defaultdict(list)
for s, al in acc.items():
    for a, b in itertools.combinations(al, 2):
        v = ani(a, b)
        if v is not None: within[s].append(v)

between = collections.defaultdict(list)
specs = sorted(acc)
for i in range(len(specs)):
    for j in range(i + 1, len(specs)):
        pair = tuple(sorted((specs[i], specs[j])))
        for a in acc[specs[i]]:
            for b in acc[specs[j]]:
                v = ani(a, b)
                if v is not None: between[pair].append(v)

def stat(vs):
    vs = sorted(vs)
    return (f"n={len(vs):3d} min={vs[0]:.2f} med={vs[len(vs)//2]:.2f} max={vs[-1]:.2f}"
            if vs else "n=0")

print("\n=== WITHIN-species CDS ANI (conspecific divergence) ===")
for s in sorted(within):
    print(f"  {s:28s} {stat(within[s])}")
allw = [v for vs in within.values() for v in vs]
print(f"  {'ALL within':28s} {stat(allw)}")

print("\n=== BETWEEN-species CDS ANI (congeneric divergence) ===")
for pair in sorted(between):
    print(f"  {pair[0].split()[1]:12s} vs {pair[1].split()[1]:14s} {stat(between[pair])}")
allb = [v for vs in between.values() for v in vs]
print(f"  {'ALL between':28s} {stat(allb)}")

if allw and allb:
    wmin, bmax = min(allw), max(allb)
    print(f"\n=== GAP CHECK ===")
    print(f"  lowest within-species ANI  : {wmin:.2f}")
    print(f"  highest between-species ANI: {bmax:.2f}")
    if wmin > bmax:
        print(f"  CLEAN GAP of {wmin-bmax:.2f} ANI points: within and between do not overlap.")
        print(f"  -> species are separable by divergence magnitude.")
    else:
        print(f"  OVERLAP of {bmax-wmin:.2f} points: some between-species pairs are closer than")
        print(f"     the loosest within-species pair. Separation is marginal for those.")
