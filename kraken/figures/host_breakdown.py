#!/usr/bin/env python3
"""
host_breakdown.py — horizontal bar chart of top plant hosts in runs.tsv.
Shows biosample-representative rows only; broad clade terms excluded.

Run from crypt/:
    python kraken/figures/host_breakdown.py
"""

import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
from pathlib import Path

RUNS_TSV = Path("stat/output/stat_filter/data/runs.tsv")
OUT_DIR  = Path("kraken/output/figures/host_breakdown")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_N = 35

BROAD = {
    "Viridiplantae","Mesangiospermae","BOP clade","Triticinae","IRL clade",
    "PACMAD clade","eudicotyledons","Triticodae","NPAAA clade","Eukaryota",
    "Embryophyta","Tracheophyta","core eudicotyledons","Oryza","Aegilops",
}

# Clade groupings for colour
CLADE_COLOURS = {
    "Cereals":    "#e6a817",
    "Legumes":    "#4a9e5c",
    "Solanaceae": "#c0392b",
    "Brassicas":  "#8e44ad",
    "Other":      "#7f8c8d",
}

def assign_clade(name: str) -> str:
    n = name.lower()
    if any(x in n for x in ["triticum","aegilops","hordeum","oryza","zea","brachypodium",
                              "sorghum","avena","secale","musa","panicum","setaria"]):
        return "Cereals"
    if any(x in n for x in ["glycine","phaseolus","cicer","vigna","medicago","pisum",
                              "cajanus","lens","arachis","vicia"]):
        return "Legumes"
    if any(x in n for x in ["solanum","lycopersicum","tuberosum","capsicum","nicotiana",
                              "petunia","tomato"]):
        return "Solanaceae"
    if any(x in n for x in ["arabidopsis","brassica","raphanus","sinapis"]):
        return "Brassicas"
    return "Other"


rows = list(csv.DictReader(open(RUNS_TSV), delimiter='\t'))
rep  = [r for r in rows if r.get('biosample_representative') == 'True']

counts = Counter(r['host'] for r in rep if r['host'] not in BROAD)
top    = counts.most_common(TOP_N)
names  = [h for h, _ in top]
vals   = [n for _, n in top]
clades = [assign_clade(h) for h in names]
colors = [CLADE_COLOURS[c] for c in clades]

# Reverse so largest is at top
names, vals, colors, clades = names[::-1], vals[::-1], colors[::-1], clades[::-1]

fig, ax = plt.subplots(figsize=(8, 10))
bars = ax.barh(range(len(names)), vals, color=colors, height=0.7, edgecolor="white", linewidth=0.4)

# Value labels
for i, v in enumerate(vals):
    ax.text(v + 8, i, f"{v:,}", va="center", ha="left", fontsize=7.5, color="#2a2a2a")

ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel("Biosample-representative runs", fontsize=10)
ax.set_title(f"Top {TOP_N} plant hosts in STAT-screened RNA-seq runs", fontsize=11, fontweight="bold")
ax.set_xlim(0, max(vals) * 1.18)
ax.spines[["top","right","left"]].set_visible(False)
ax.tick_params(left=False)
ax.xaxis.grid(True, color="#e0e0e0", linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

legend_patches = [mpatches.Patch(color=v, label=k) for k, v in CLADE_COLOURS.items()]
ax.legend(handles=legend_patches, loc="lower right", fontsize=8, framealpha=0.85,
          title="Clade", title_fontsize=8)

n_total = sum(counts.values())
ax.text(0.99, 0.01, f"n = {n_total:,} total named-host BioSamples",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="#555")

plt.tight_layout()
for ext in ("pdf", "png"):
    p = OUT_DIR / f"host_breakdown.{ext}"
    fig.savefig(p, dpi=180, bbox_inches="tight")
    print(f"Written: {p}")
