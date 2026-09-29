#!/usr/bin/env python3
"""Per-species points + per-host tightness stats for the individual complexity images."""
import csv, sys, math, statistics, collections
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "kraken" / "assign"))
from kraken_report_analysis import load_reference, load_run_to_biosample, parse_report  # noqa
csv.field_size_limit(10 ** 7)
pathogens, hosts = load_reference()
# read-based host call per run (kraken_host_call.py), which resolves all 3,220 runs;
# the LLM host label left ~450 unresolved
host_of = {}
for r in csv.DictReader(open(ROOT / "kraken/assign/host/data/host_calls.tsv"), delimiter="\t"):
    host_of[r["run"].strip()] = (r.get("host") or "").strip()

pts = collections.defaultdict(lambda: collections.defaultdict(list))
for p in sorted((ROOT / "kraken/assign/data/reports").glob("*.txt")):
    h = host_of.get(p.stem, "") or "unresolved"
    for rec in parse_report(p):
        if rec["rank"] != "S" or rec["distinct"] is None:
            continue
        tid = str(rec["taxid"])
        if tid in hosts or tid not in pathogens or rec["reads_clade"] < 100 or rec["distinct"] < 1:
            continue
        pts[rec["name"]][h].append((rec["reads_clade"], rec["distinct"]))

def fit(ps):
    xs = [math.log10(r) for r, _ in ps]; ys = [math.log10(d) for _, d in ps]
    mx = statistics.mean(xs); my = statistics.mean(ys)
    sxx = sum((x-mx)**2 for x in xs); syy = sum((y-my)**2 for y in ys)
    sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    if sxx == 0 or syy == 0: return None, None
    return sxy/sxx, (sxy**2)/(sxx*syy)

# keep species with >=1 fittable host cell (n>=15, span>=0.7)
keep = set()
stats = []
for sp, hs in pts.items():
    for h, ps in hs.items():
        if len(ps) < 15: continue
        span = max(math.log10(r) for r,_ in ps) - min(math.log10(r) for r,_ in ps)
        if span < 0.7: continue
        s, r2 = fit(ps)
        stats.append((sp, h, len(ps), round(s,3), round(r2,3)))
        keep.add(sp)

with open(ROOT/"complexity_by_host/_points.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["species","host","reads","distinct"])
    for sp in keep:
        for h, ps in pts[sp].items():
            for r, d in ps: w.writerow([sp, h, r, d])
with open(ROOT/"complexity_by_host/_stats.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["species","host","n","slope","r2"])
    for row in stats: w.writerow(row)
print(f"  {len(keep)} species, {len(stats)} fittable host cells")
