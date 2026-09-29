#!/usr/bin/env python3
"""
prep_slopes.py — per-species and per-(host x species) accumulation slopes.

slope = fit of log(distinct minimizers) ~ log(reads) over a taxon's detections. A real
detection accumulates k-mers (slope ~0.7-0.9); a pile-up or cross-map plateaus
(slope <0.5). Two outputs feed the slope figures:

  slopes_cohort.tsv  one slope per pathogen species over the whole cohort
  slopes_byhost.tsv  slope per (host, species) for species seen on >=2 major hosts,
                     which is what shows the filter MUST be per-host (Bipolaris maydis
                     is real on maize, cross-maps on wheat)
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

def fit(ps):
    xs = [x for x, _ in ps]; ys = [y for _, y in ps]
    mx = statistics.mean(xs); my = statistics.mean(ys)
    den = sum((x - mx) ** 2 for x in xs)
    return None if den == 0 else sum((x - mx) * (y - my) for x, y in ps) / den

cohort = collections.defaultdict(list)
byhost = collections.defaultdict(lambda: collections.defaultdict(list))
for p in sorted((ROOT / "kraken/assign/data/reports").glob("*.txt")):
    h = llm.get(run2bs.get(p.stem, ""), "")
    for rec in parse_report(p):
        if rec["rank"] != "S" or rec["distinct"] is None:
            continue
        tid = str(rec["taxid"])
        if tid in hosts or tid not in pathogens or rec["reads_clade"] < 100 or rec["distinct"] < 1:
            continue
        pt = (math.log10(rec["reads_clade"]), math.log10(rec["distinct"]))
        cohort[rec["name"]].append(pt)
        if h:
            byhost[rec["name"]][h].append(pt)

MIN_N, MIN_SPAN = 30, 1.0
with open(ROOT / "kraken/assign/figures/slopes_cohort.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["species", "n", "slope"])
    for sp, ps in cohort.items():
        span = max(x for x, _ in ps) - min(x for x, _ in ps)
        s = fit(ps) if len(ps) >= MIN_N and span >= MIN_SPAN else None
        if s is not None:
            w.writerow([sp, len(ps), round(s, 3)])

# per-host: keep species present on >=2 major hosts with enough points each
MAJOR = ["Triticum aestivum", "Zea mays", "Hordeum vulgare", "Oryza sativa",
         "Glycine max", "Solanum lycopersicum", "Vitis vinifera", "Brassica napus"]
MIN_HN, MIN_HSPAN = 15, 0.7
rows = []
for sp, hs in byhost.items():
    for h in MAJOR:
        ps = hs.get(h, [])
        if len(ps) < MIN_HN:
            continue
        span = max(x for x, _ in ps) - min(x for x, _ in ps)
        s = fit(ps) if span >= MIN_HSPAN else None
        if s is not None:
            rows.append((sp, h, len(ps), round(s, 3)))
keep = {sp for sp in {r[0] for r in rows}
        if sum(1 for r in rows if r[0] == sp) >= 2}
with open(ROOT / "kraken/assign/figures/slopes_byhost.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["species", "host", "n", "slope"])
    for r in rows:
        if r[0] in keep:
            w.writerow(r)
print(f"  slopes_cohort.tsv : {sum(1 for _ in open(ROOT/'kraken/assign/figures/slopes_cohort.tsv'))-1} species")
print(f"  slopes_byhost.tsv : {len(keep)} species across major hosts")
