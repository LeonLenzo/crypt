#!/usr/bin/env python3
"""
prep_complexity.py — reads vs distinct minimizers per wheat detection, for the
library-complexity view of the criterion.

The k-mer-fraction criterion (distinct / taxon's total in db_v3) removes pile-up
artefacts but is biased by genome size and misses that real detections form a diagonal
(distinct grows with reads) while artefacts form a flat floor (distinct pinned, reads
pile up). The discriminating quantity is distinct minimizers PER READ, which needs no
denominator and so also escapes the Alternaria genome-size bias.

Output: kraken/assign/figures/complexity.tsv  (Triticum aestivum)
"""
import csv, sys
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

# taxa that the read-level alignment / absence-of-competitor work flagged as
# wheat-context artefacts; used only to colour the demonstration, not to threshold
ART = {"Phakopsora pachyrhizi", "Melampsora laricis-populina", "Cronartium quercuum",
       "Austropuccinia psidii", "Cercospora zeae-maydis", "Puccinia sorghi"}

rows = []
for p in sorted((ROOT / "kraken/assign/data/reports").glob("*.txt")):
    if llm.get(run2bs.get(p.stem, "")) != "Triticum aestivum":
        continue
    for rec in parse_report(p):
        if rec["rank"] != "S" or rec["distinct"] is None:
            continue
        tid = str(rec["taxid"])
        if tid in hosts or tid not in pathogens:
            continue
        if rec["reads_clade"] < 100 or rec["distinct"] < 1:
            continue
        rows.append({
            "reads": rec["reads_clade"], "distinct": rec["distinct"],
            "ratio": round(rec["distinct"] / rec["reads_clade"], 6),
            "kind": "artefact taxa (wheat context)" if rec["name"] in ART
                    else "other wheat pathogens",
        })

out = ROOT / "kraken/assign/figures/complexity.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"  wrote {out}  ({len(rows)} detections)")
