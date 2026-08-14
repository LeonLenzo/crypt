#!/usr/bin/env python3
"""
Literature resolution pipeline Sankey — temporal funnel.

Left-to-right funnel: each strategy stage receives the unresolved remainder
from the previous stage. Resolved BPs exit sideways into text coverage buckets.

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

# ── Classify each BP by resolution stage and text coverage ───────────────────

STAGES = ["fetch_lit.py", "Web search (Serper)", "Manual DOI lookup", "Unresolved"]

COVERAGE = [
    "Methods + abstract",
    "Abstract (PMID)",
    "Abstract (DOI only)",
    "Title only",
    "No publication data",
]

# stage_flows[stage][coverage] = count
stage_flows: dict[str, dict[str, int]] = {s: defaultdict(int) for s in STAGES}

for bp in bps_in_runs:
    e   = cache.get(bp, {})
    src = e.get("pmid_source", "")
    pmid  = bool(e.get("primary_pmid"))
    doi   = bool(e.get("primary_doi"))
    abstr = bool(e.get("abstract"))
    meth  = bool(e.get("methods_text"))

    if "Serper" in src:
        stage = "Web search (Serper)"
    elif src == "manual_doi":
        stage = "Manual DOI lookup"
    elif not pmid and not doi:
        stage = "Unresolved"
    else:
        stage = "fetch_lit.py"

    if meth:              cov = "Methods + abstract"
    elif abstr and pmid:  cov = "Abstract (PMID)"
    elif abstr and doi:   cov = "Abstract (DOI only)"
    elif pmid or doi:     cov = "Title only"
    else:                 cov = "No publication data"

    stage_flows[stage][cov] += 1

# How many BPs enter each stage (cumulative remainder)
stage_total = {s: sum(stage_flows[s].values()) for s in STAGES}
stage_cumulative_in = {}
remaining = total_bps
for s in STAGES:
    stage_cumulative_in[s] = remaining
    remaining -= stage_total[s]

# ── Node definitions ──────────────────────────────────────────────────────────
#
# Nodes:
#   0   : Origin ("All BioProjects")
#   1–3 : Strategy stages (fetch_lit, Serper, Manual)
#   4   : "Unresolved" exit node (stage 4)
#   5–9 : Text coverage output nodes
#
# Layout (x positions, 0=left, 1=right):
#   x=0.01  Origin
#   x=0.25  fetch_lit.py
#   x=0.50  Serper
#   x=0.75  Manual DOI
#   x=0.92  Unresolved  (exit)
#   x=0.99  Text coverage outputs (right edge)

STAGE_COLORS = {
    "fetch_lit.py":          "#2980b9",
    "Web search (Serper)":   "#e67e22",
    "Manual DOI lookup":     "#8e44ad",
    "Unresolved":            "#c0392b",
}

COVERAGE_COLORS = {
    "Methods + abstract":   "#1a5e36",
    "Abstract (PMID)":      "#27ae60",
    "Abstract (DOI only)":  "#16a085",
    "Title only":           "#f39c12",
    "No publication data":  "#7f8c8d",
}

ORIGIN_COLOR = "#2c3e50"

# Node index mapping
ORIGIN_IDX   = 0
STAGE_IDX    = {s: i + 1 for i, s in enumerate(STAGES[:3])}   # 1,2,3
UNRESOLV_IDX = 4
COV_IDX      = {c: i + 5 for i, c in enumerate(COVERAGE)}     # 5–9

n_nodes = 10
node_labels = [""] * n_nodes
node_colors = [""] * n_nodes
node_x      = [0.0] * n_nodes
node_y      = [0.5] * n_nodes

# Origin
node_labels[ORIGIN_IDX] = f"<b>All BioProjects</b><br>{total_bps:,} BPs"
node_colors[ORIGIN_IDX] = ORIGIN_COLOR
node_x[ORIGIN_IDX]      = 0.01
node_y[ORIGIN_IDX]      = 0.5

# Stage nodes
stage_x = {"fetch_lit.py": 0.25, "Web search (Serper)": 0.50, "Manual DOI lookup": 0.75}
for s, idx in STAGE_IDX.items():
    n_in = stage_cumulative_in[s]
    node_labels[idx] = f"<b>{s}</b><br>{n_in:,} BPs in"
    node_colors[idx] = STAGE_COLORS[s]
    node_x[idx]      = stage_x[s]
    node_y[idx]      = 0.5

# Unresolved exit
n_unresolved = stage_total["Unresolved"]
node_labels[UNRESOLV_IDX] = f"<b>Unresolved</b><br>{n_unresolved:,} BPs"
node_colors[UNRESOLV_IDX] = STAGE_COLORS["Unresolved"]
node_x[UNRESOLV_IDX]      = 0.92
node_y[UNRESOLV_IDX]      = 0.88

# Coverage output nodes (stacked top→bottom by size)
cov_totals = defaultdict(int)
for s in STAGES:
    for cov, v in stage_flows[s].items():
        cov_totals[cov] += v

cov_order = sorted(COVERAGE, key=lambda c: -cov_totals[c])
# y positions: evenly spaced with Methods+abstract at top
y_positions = [0.05, 0.22, 0.40, 0.55, 0.68]
for rank, cov in enumerate(cov_order):
    idx = COV_IDX[cov]
    node_labels[idx] = f"<b>{cov}</b><br>{cov_totals[cov]:,} BPs"
    node_colors[idx] = COVERAGE_COLORS[cov]
    node_x[idx]      = 0.99
    node_y[idx]      = y_positions[rank]

# ── Build links ───────────────────────────────────────────────────────────────

link_sources, link_targets, link_values, link_colors_list = [], [], [], []

def rgba(hex_color: str, alpha: float = 0.38) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"

def add_link(src: int, tgt: int, val: int, color: str) -> None:
    if val <= 0:
        return
    link_sources.append(src)
    link_targets.append(tgt)
    link_values.append(val)
    link_colors_list.append(rgba(color))

# Origin → Stage 1
add_link(ORIGIN_IDX, STAGE_IDX["fetch_lit.py"], total_bps, ORIGIN_COLOR)

# Each stage: resolved exits to coverage nodes, remainder passes to next stage
remaining_stages = ["fetch_lit.py", "Web search (Serper)", "Manual DOI lookup"]
next_stage = {
    "fetch_lit.py":         "Web search (Serper)",
    "Web search (Serper)":  "Manual DOI lookup",
    "Manual DOI lookup":    None,
}

for s in remaining_stages:
    s_idx   = STAGE_IDX[s]
    s_color = STAGE_COLORS[s]
    # Resolved: flow to coverage nodes
    for cov, count in stage_flows[s].items():
        add_link(s_idx, COV_IDX[cov], count, s_color)
    # Pass remainder forward
    nxt = next_stage[s]
    n_remain = stage_cumulative_in[s] - stage_total[s]
    if nxt and n_remain > 0:
        add_link(s_idx, STAGE_IDX[nxt], n_remain, s_color)

# Final stage: Manual DOI → Unresolved exit
add_link(STAGE_IDX["Manual DOI lookup"], UNRESOLV_IDX,
         stage_flows["Unresolved"]["No publication data"],
         STAGE_COLORS["Unresolved"])

# ── Figure ────────────────────────────────────────────────────────────────────

resolved = total_bps - n_unresolved

fig = go.Figure(go.Sankey(
    arrangement="fixed",
    node=dict(
        pad=14,
        thickness=20,
        line=dict(color="white", width=0.5),
        label=node_labels,
        color=node_colors,
        x=node_x,
        y=node_y,
    ),
    link=dict(
        source=link_sources,
        target=link_targets,
        value=link_values,
        color=link_colors_list,
    ),
))

fig.update_layout(
    title=dict(
        text=(f"Literature resolution pipeline — {total_bps:,} BioProjects<br>"
              f"<sup>Resolved: {resolved:,} ({100*resolved/total_bps:.1f}%)   "
              f"Unresolved: {n_unresolved:,} "
              f"({100*n_unresolved/total_bps:.1f}%)</sup>"),
        font=dict(size=16),
        x=0.01,
    ),
    font=dict(size=12, family="Arial"),
    paper_bgcolor="white",
    width=1050,
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
