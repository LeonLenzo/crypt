#!/usr/bin/env python3
"""
prep_detection_panels.py — data for detection_panels.R.

One row per species-rank detection clearing a 100-read floor, tagged by whether the
declaring study named that pathogen. Splitting on the study's own declaration rather
than on species is what makes the figure species-agnostic: it asks what a known-real
infection looks like, and then how everything else compares against that yardstick.

Reuses the loaders in kraken_report_analysis.py rather than re-deriving the joins,
so the declared test stays ancestry-aware in both places.

Run from crypt/:
  python kraken/assign/figures/prep_detection_panels.py
Output: kraken/assign/figures/detection_panels.tsv
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "kraken" / "assign"))
from kraken_report_analysis import (          # noqa: E402
    load_reference, load_inspect, load_run_to_biosample, load_samples, parse_report,
    expand_declared)

MIN_READS = 100
CRITERION = 0.01

pathogens, hosts = load_reference()
db_min = load_inspect(ROOT / "kraken/build/data/db_v3_inspect.tsv")
run2bs = load_run_to_biosample()
samples = load_samples()
if not db_min:
    sys.exit("no db_v3_inspect.tsv — the k-mer fraction cannot be computed")

# host comes along so the panels can be split by it: the cohort is 73 species, not
# wheat, and mixing them lets one host's pathogens read as another's noise
import csv as _csv                                                        # noqa: E402
host_of = {}
with open(ROOT / "metadata/classify/data/samples.tsv", newline="") as fh:
    for row in _csv.DictReader(fh, delimiter="\t"):
        host_of[row["BioSample"].strip()] = (
            (row.get("llm_host_resolved") or row.get("stat_host") or "").strip() or "unresolved")

rows = []
for p in sorted((ROOT / "kraken/assign/data/reports").glob("*.txt")):
    run = p.stem
    bs = run2bs.get(run, "")
    host = host_of.get(bs, "unresolved")
    recs = list(parse_report(p))
    # symmetric: a study declaring a forma specialis still counts when Kraken2 calls
    # the species. Without this every barley study reads as undeclared.
    declared = expand_declared(recs, samples.get(bs, {}).get("declared", set()))
    for rec in recs:
        if rec["rank"] != "S" or rec["distinct"] is None:
            continue
        tid = str(rec["taxid"])
        # Hosts stay out: a wheat call in a wheat study is the library, not a finding.
        # Non-pathogens stay IN, as their own class. They are mostly the congeners
        # db_v3 added as LCA competitors, and the analysis has been blind to them
        # because PHI-base supplies only 665 pathogen taxids against db_v3's 1,043.
        if tid in hosts:
            continue
        if rec["reads_clade"] < MIN_READS:
            continue
        total = db_min.get(rec["taxid"])
        if not total:
            continue
        frac = rec["distinct"] / total
        if frac <= 0:
            continue
        lineage = {tid} | {str(t) for _, t, _ in rec["lineage"]}
        if lineage & declared:
            cls = "declared"
        elif tid in pathogens:
            cls = "secondary pathogen"
        else:
            cls = "secondary non-pathogen"
        rows.append({
            "run": run, "host": host, "taxon": rec["name"], "reads": rec["reads_clade"],
            "kmer_frac": round(frac, 8),
            "class": cls,
            "passes": "yes" if frac >= CRITERION else "no",
        })

out = ROOT / "kraken/assign/figures/detection_panels.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

import collections                              # noqa: E402
c = collections.Counter((r["class"], r["passes"]) for r in rows)
for d in ("declared", "secondary pathogen", "secondary non-pathogen"):
    a, b = c[(d, "yes")], c[(d, "no")]
    print(f"  {d:<24} n={a+b:>6}  above criterion={a:>5} ({100*a/(a+b):.1f}%)")
print(f"  wrote {out}")
