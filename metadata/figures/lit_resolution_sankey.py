#!/usr/bin/env python3
"""
Literature resolution pipeline Sankey.

3-column flow: Origin → Resolution strategy → Text coverage.
Shows how 1,287 BioProjects were resolved and what text data is available.

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

total_bps = len(bps_in_runs)

# ── Classify each BP ─────────────────────────────────────────────────────────

flows: list[tuple[str, str]] = []
for bp in bps_in_runs:
    e   = cache.get(bp, {})
    src = e.get("pmid_source", "")
    pmid  = bool(e.get("primary_pmid"))
    doi   = bool(e.get("primary_doi"))
    abstr = bool(e.get("abstract"))
    meth  = bool(e.get("methods_text"))

    if "PMC" in src:               path = "PMC full-text"
    elif src == "BioProject XML":  path = "BioProject XML"
    elif "Serper" in src:          path = "Web search (Serper)"
    elif src == "manual_doi":      path = "Manual DOI lookup"
    elif not pmid and not doi:     path = "Unresolved"
    else:                          path = "Other / legacy"

    if meth:              text = "Methods + abstract"
    elif abstr and pmid:  text = "Abstract (PMID)"
    elif abstr and doi:   text = "Abstract (DOI only)"
    elif pmid or doi:     text = "Title only"
    else:                 text = "No publication data"

    flows.append((path, text))

counts: dict[tuple[str, str], int] = defaultdict(int)
for pair in flows:
    counts[pair] += 1

# ── Node definitions ──────────────────────────────────────────────────────────

ORIGIN = "All BioProjects"

STRATEGY_NODES = [
    "PMC full-text",
    "BioProject XML",
    "Web search (Serper)",
    "Manual DOI lookup",
    "Other / legacy",
    "Unresolved",
]

COVERAGE_NODES = [
    "Methods + abstract",
    "Abstract (PMID)",
    "Abstract (DOI only)",
    "Title only",
    "No publication data",
]

ALL_NODES  = [ORIGIN] + STRATEGY_NODES + COVERAGE_NODES
node_index = {n: i for i, n in enumerate(ALL_NODES)}

STRATEGY_COLORS = {
    "PMC full-text":        "#27ae60",
    "BioProject XML":       "#2980b9",
    "Web search (Serper)":  "#e67e22",
    "Manual DOI lookup":    "#8e44ad",
    "Other / legacy":       "#95a5a6",
    "Unresolved":           "#c0392b",
}

COVERAGE_COLORS = {
    "Methods + abstract":   "#1a5e36",
    "Abstract (PMID)":      "#27ae60",
    "Abstract (DOI only)":  "#16a085",
    "Title only":           "#f39c12",
    "No publication data":  "#7f8c8d",
}

ORIGIN_COLOR = "#2c3e50"

node_colors = (
    [ORIGIN_COLOR]
    + [STRATEGY_COLORS[n] for n in STRATEGY_NODES]
    + [COVERAGE_COLORS[n] for n in COVERAGE_NODES]
)

# ── Build links ───────────────────────────────────────────────────────────────

link_sources, link_targets, link_values, link_colors = [], [], [], []

def add_link(src_name: str, tgt_name: str, value: int, color_hex: str) -> None:
    if value <= 0:
        return
    link_sources.append(node_index[src_name])
    link_targets.append(node_index[tgt_name])
    link_values.append(value)
    r, g, b = int(color_hex[1:3], 16), int(color_hex[3:5], 16), int(color_hex[5:7], 16)
    link_colors.append(f"rgba({r},{g},{b},0.38)")

# Layer 1: Origin → Strategies
left_totals: dict[str, int] = defaultdict(int)
for (src, _), v in counts.items():
    left_totals[src] += v

for strat in STRATEGY_NODES:
    add_link(ORIGIN, strat, left_totals[strat], STRATEGY_COLORS[strat])

# Layer 2: Strategies → Coverage
for (src_name, tgt_name), val in counts.items():
    add_link(src_name, tgt_name, val, STRATEGY_COLORS[src_name])

# ── Node labels ───────────────────────────────────────────────────────────────

right_totals: dict[str, int] = defaultdict(int)
for (_, tgt), v in counts.items():
    right_totals[tgt] += v

node_labels = [f"<b>{ORIGIN}</b><br>{total_bps:,} BPs"]
for n in STRATEGY_NODES:
    node_labels.append(f"<b>{n}</b><br>{left_totals[n]:,} BPs")
for n in COVERAGE_NODES:
    node_labels.append(f"<b>{n}</b><br>{right_totals[n]:,} BPs")

# ── Figure ────────────────────────────────────────────────────────────────────

resolved  = total_bps - left_totals["Unresolved"]

fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(
        pad=16,
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
    width=1000,
    height=580,
    margin=dict(l=10, r=10, t=85, b=10),
)

fig.write_html(OUT_HTML)
print(f"Written {OUT_HTML}")

try:
    fig.write_image(OUT_PNG, scale=2)
    print(f"Written {OUT_PNG}")
except Exception as e:
    print(f"PNG skipped ({e})")
