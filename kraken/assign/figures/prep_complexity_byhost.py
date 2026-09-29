#!/usr/bin/env python3
"""
prep_complexity_byhost.py — reads vs distinct minimizers, per species, coloured by host.

For every species that appears in the per-host slope heatmap (slopes_byhost.tsv), emit
its raw (reads, distinct, host) detections across the 8 major hosts, plus the fitted
per-host slopes. The scatter then shows WITHIN each species panel how the host clouds
separate: e.g. Bipolaris maydis maize points ride the diagonal, wheat points sit on the
floor, so the per-host slope difference from the heatmap is visible as diverging clouds.

Output: kraken/assign/figures/complexity_byhost.tsv
        kraken/assign/figures/complexity_byhost_slopes.tsv  (per host x species slope)
"""
import csv, sys, math, statistics, collections
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "kraken" / "assign"))
from kraken_report_analysis import load_reference, load_run_to_biosample, parse_report  # noqa

csv.field_size_limit(10 ** 7)
pathogens, hosts = load_reference()
run2bs = load_run_to_biosample()
llm = {}
for r in csv.DictReader(open(ROOT / "metadata/classify/data/samples.tsv"), delimiter="\t"):
    llm[r["BioSample"].strip()] = (r.get("llm_host_resolved") or "").strip()

MAJOR = ["Triticum aestivum", "Hordeum vulgare", "Zea mays", "Oryza sativa",
         "Glycine max", "Solanum lycopersicum", "Vitis vinifera", "Brassica napus"]
# species to show = those in the heatmap
keep_sp = set()
for r in csv.DictReader(open(ROOT / "kraken/assign/figures/slopes_byhost.tsv"), delimiter="\t"):
    keep_sp.add(r["species"])

pts = collections.defaultdict(lambda: collections.defaultdict(list))   # sp -> host -> [(r,d)]
for p in sorted((ROOT / "kraken/assign/data/reports").glob("*.txt")):
    h = llm.get(run2bs.get(p.stem, ""), "")
    if h not in MAJOR:
        continue
    for rec in parse_report(p):
        if rec["rank"] != "S" or rec["distinct"] is None:
            continue
        tid = str(rec["taxid"])
        if tid in hosts or tid not in pathogens or rec["reads_clade"] < 100 or rec["distinct"] < 1:
            continue
        if rec["name"] in keep_sp:
            pts[rec["name"]][h].append((rec["reads_clade"], rec["distinct"]))

def fit(ps):
    xs = [math.log10(r) for r, _ in ps]; ys = [math.log10(d) for _, d in ps]
    mx = statistics.mean(xs); my = statistics.mean(ys)
    den = sum((x - mx) ** 2 for x in xs)
    return None if den == 0 else sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den

with open(ROOT / "kraken/assign/figures/complexity_byhost.tsv", "w", newline="") as fh, \
     open(ROOT / "kraken/assign/figures/complexity_byhost_slopes.tsv", "w", newline="") as sh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["species", "host", "reads", "distinct"])
    ws = csv.writer(sh, delimiter="\t"); ws.writerow(["species", "host", "n", "slope"])
    for sp in keep_sp:
        for h in MAJOR:
            ps = pts[sp].get(h, [])
            for r, d in ps:
                w.writerow([sp, h, r, d])
            if len(ps) >= 15:
                span = max(math.log10(r) for r, _ in ps) - min(math.log10(r) for r, _ in ps)
                s = fit(ps) if span >= 0.7 else None
                if s is not None:
                    ws.writerow([sp, h, len(ps), round(s, 2)])
print(f"  {len(keep_sp)} species written")
