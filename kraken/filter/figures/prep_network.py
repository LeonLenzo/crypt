#!/usr/bin/env python3
"""Co-occurrence ASSOCIATION edges/nodes per host from the filtered detections.

A pathogen is present in a biosample if any of its runs scores >= THRESH. The universe is
every biosample of that host (from the host calls), so absence counts too. For each pair we
test whether they co-occur MORE than their individual prevalences predict (Fisher's exact,
one-sided) rather than just counting raw co-occurrences: two common pathogens land in the
same sample by chance, and a raw-count edge would link everything to everything. Edge weight
is fold-enrichment over expected; keep the pair's p for significance filtering in R.
"""
import csv, collections, itertools, math
from pathlib import Path
from scipy.stats import fisher_exact
ROOT = Path(__file__).resolve().parents[3]
THRESH = 0.7       # real-shape score: not an artefact
COVFLOOR = 0.01    # genome k-mer coverage in the sample (~50k distinct k-mers): DNA present, not a trace
csv.field_size_limit(10 ** 7)

klass = {r["label"]: r["class"] for r in
         csv.DictReader(open(ROOT / "kraken/filter/data/patho_lineage.tsv"), delimiter="\t")}

# universe: every biosample seen per host (denominator for the null), from the host calls
universe = collections.defaultdict(set)
for r in csv.DictReader(open(ROOT / "kraken/assign/host/data/host_calls.tsv"), delimiter="\t"):
    bs = (r.get("biosample") or "").strip(); h = (r.get("host") or "").strip()
    if bs and h:
        universe[h].add(bs)

present = collections.defaultdict(lambda: collections.defaultdict(set))   # host -> bs -> {species}
for r in csv.DictReader(open(ROOT / "kraken/filter/data/_classified.tsv"), delimiter="\t"):
    try:
        s = float(r["score"]); cov = float(r["coverage"])
    except ValueError:
        continue
    bs = r["biosample"].strip()
    if s >= THRESH and cov >= COVFLOOR and bs:
        present[r["host"]][bs].add(r["species"])

def bh(ps):                                       # Benjamini-Hochberg adjusted p-values
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); adj = [0.0]*m; prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        k = m - rank + 1
        prev = min(prev, ps[i]*m/k); adj[i] = prev
    return adj

nodes, edges, summary = [], [], []
for host, N in ((h, len(universe[h])) for h in present):
    bsmap = present[host]
    prev = collections.Counter()
    for sps in bsmap.values():
        for sp in sps:
            prev[sp] += 1
    obs = collections.Counter()
    for sps in bsmap.values():
        for a, b in itertools.combinations(sorted(sps), 2):
            obs[(a, b)] += 1
    for sp, n in prev.items():
        nodes.append([host, sp, n, klass.get(sp, "other")])
    rows, ps = [], []
    for (a, b), o in obs.items():
        if o < 3:                                 # min support
            continue
        both = o; aonly = prev[a]-o; bonly = prev[b]-o; neither = N-aonly-bonly-both
        if neither < 0:
            continue
        exp = prev[a]*prev[b]/N
        p = fisher_exact([[both, aonly], [bonly, neither]], alternative="greater")[1]
        l2 = math.log2((both+0.5)/(exp+0.5))
        rows.append([host, a, b, both, round(exp, 1), round(l2, 2), p]); ps.append(p)
    for row, padj in zip(rows, bh(ps) if ps else []):
        edges.append(row + [round(padj, 4)])
    summary.append([host, N, sum(1 for s in bsmap.values() if len(s) >= 2), len(rows)])

def dump(path, header, rows):
    with open(ROOT / path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(header); w.writerows(rows)

dump("kraken/filter/data/network_nodes.tsv", ["host", "species", "n_biosamples", "class"], nodes)
dump("kraken/filter/data/network_edges.tsv",
     ["host", "source", "target", "obs", "expected", "log2_oe", "p", "p_adj"], edges)
dump("kraken/filter/data/network_summary.tsv", ["host", "bs_total", "bs_coinfected", "n_tested"], summary)
print(f"  {len(present)} hosts, {len(nodes)} host-nodes, {len(edges)} tested pairs (thresh {THRESH})")
