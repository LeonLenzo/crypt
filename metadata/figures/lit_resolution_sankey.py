#!/usr/bin/env python3
"""
Literature resolution pipeline Sankey.

Shows how 1,287 BioProjects were resolved through each strategy and what
text data is available downstream.

Run from crypt/: python metadata/figures/lit_resolution_sankey.py
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

try:
    import plotly.graph_objects as go
except ImportError:
    sys.exit("pip install plotly kaleido")

CACHE_PATH = Path("metadata/output/fetch_lit/data/lit_cache.json")
RUNS_TSV   = Path("stat/output/filter_runs/data/runs.tsv")
OUT_DIR    = Path("metadata/output/figures/sankey")
OUT_HTML   = OUT_DIR / "lit_resolution_sankey.html"
OUT_PNG    = OUT_DIR / "lit_resolution_sankey.png"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────

with open(CACHE_PATH) as f:
    cache = json.load(f)

bps_in_runs: set[str] = set()
with open(RUNS_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        bps_in_runs.add(row["BioProject"])

# ── Classify each BP into (resolution_path, text_coverage) ───────────────────

flows: list[tuple[str, str]] = []
for bp in bps_in_runs:
    e   = cache.get(bp, {})
    src = e.get("pmid_source", "")
    pmid  = bool(e.get("primary_pmid"))
    doi   = bool(e.get("primary_doi"))
    abstr = bool(e.get("abstract"))
    meth  = bool(e.get("methods_text"))

    if "PMC" in src:
        path = "PMC full-text"
    elif src == "BioProject XML":
        path = "BioProject XML"
    elif "Serper" in src:
        path = "Web search (Serper)"
    elif not pmid and not doi:
        path = "Unresolved"
    else:
        path = "Other / legacy"

    if meth:
        text = "Methods + abstract"
    elif abstr and pmid:
        text = "Abstract (PMID)"
    elif abstr and doi:
        text = "Abstract (DOI only)"
    elif pmid or doi:
        text = "Title only"
    else:
        text = "No publication data"

    flows.append((path, text))

counts: dict[tuple[str, str], int] = defaultdict(int)
for pair in flows:
    counts[pair] += 1

# ── Node and link definitions ─────────────────────────────────────────────────

LEFT_NODES = [
    "PMC full-text",
    "BioProject XML",
    "Web search (Serper)",
    "Other / legacy",
    "Unresolved",
]

RIGHT_NODES = [
    "Methods + abstract",
    "Abstract (PMID)",
    "Abstract (DOI only)",
    "Title only",
    "No publication data",
]

ALL_NODES  = LEFT_NODES + RIGHT_NODES
node_index = {n: i for i, n in enumerate(ALL_NODES)}

LEFT_COLORS = {
    "PMC full-text":        "#27ae60",   # green
    "BioProject XML":       "#2980b9",   # blue
    "Web search (Serper)":  "#e67e22",   # orange
    "Other / legacy":       "#95a5a6",   # grey
    "Unresolved":           "#c0392b",   # red
}

RIGHT_COLORS = {
    "Methods + abstract":   "#1a5e36",   # dark green
    "Abstract (PMID)":      "#27ae60",   # green
    "Abstract (DOI only)":  "#16a085",   # teal
    "Title only":           "#f39c12",   # amber
    "No publication data":  "#7f8c8d",   # dark grey
}

node_colors = [LEFT_COLORS[n] for n in LEFT_NODES] + \
              [RIGHT_COLORS[n] for n in RIGHT_NODES]

# Build links
link_sources, link_targets, link_values, link_colors = [], [], [], []
for (src_name, tgt_name), val in counts.items():
    link_sources.append(node_index[src_name])
    link_targets.append(node_index[tgt_name])
    link_values.append(val)
    # link colour = left node colour, semi-transparent
    base = LEFT_COLORS[src_name].lstrip("#")
    r, g, b = int(base[0:2], 16), int(base[2:4], 16), int(base[4:6], 16)
    link_colors.append(f"rgba({r},{g},{b},0.38)")

# Node totals for labels
left_totals  = defaultdict(int)
right_totals = defaultdict(int)
for (src, tgt), v in counts.items():
    left_totals[src]  += v
    right_totals[tgt] += v

node_labels = []
for n in LEFT_NODES:
    node_labels.append(f"<b>{n}</b><br>{left_totals[n]:,} BPs")
for n in RIGHT_NODES:
    node_labels.append(f"<b>{n}</b><br>{right_totals[n]:,} BPs")

# ── Figure ────────────────────────────────────────────────────────────────────

fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(
        pad=18,
        thickness=22,
        line=dict(color="white", width=0.5),
        label=node_labels,
        color=node_colors,
    ),
    link=dict(
        source=link_sources,
        target=link_targets,
        value=link_values,
        color=link_colors,
    ),
))

total_bps = len(bps_in_runs)
resolved  = total_bps - left_totals["Unresolved"]

fig.update_layout(
    title=dict(
        text=(f"Literature resolution pipeline — {total_bps:,} BioProjects<br>"
              f"<sup>Resolved: {resolved:,} ({100*resolved/total_bps:.1f}%)   "
              f"Unresolved: {left_totals['Unresolved']:,} "
              f"({100*left_totals['Unresolved']/total_bps:.1f}%)</sup>"),
        font=dict(size=16),
        x=0.01,
    ),
    font=dict(size=12, family="Arial"),
    paper_bgcolor="white",
    width=900,
    height=560,
    margin=dict(l=10, r=10, t=80, b=10),
)

fig.write_html(OUT_HTML)
print(f"Written {OUT_HTML}")

try:
    fig.write_image(OUT_PNG, scale=2)
    print(f"Written {OUT_PNG}")
except Exception as e:
    print(f"PNG skipped ({e})")
