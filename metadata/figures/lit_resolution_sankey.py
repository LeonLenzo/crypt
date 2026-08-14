#!/usr/bin/env python3
"""
Literature resolution pipeline Sankey — waterfall funnel.

A horizontal spine (unresolved remainder) flows left-to-right through each
strategy stage. At each stage, resolved BPs drop out directly below into
coverage nodes. The spine gets thinner as BPs are resolved, creating the
funnel effect. No link crossings.

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

# ── Load and classify ─────────────────────────────────────────────────────────

with open(CACHE_PATH) as f:
    cache = json.load(f)

bps_in_runs: set[str] = set()
with open(RUNS_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        bps_in_runs.add(row["BioProject"])

total_bps = len(bps_in_runs)

STAGES = ["fetch_lit.py", "Web search (Serper)", "Manual DOI lookup"]
COVERAGE_KEYS = [
    "Methods + abstract",
    "Abstract (PMID)",
    "Abstract (DOI only)",
    "Title only",
]

stage_flows: dict[str, dict[str, int]] = {s: defaultdict(int) for s in STAGES}
unresolved = 0

for bp in bps_in_runs:
    e   = cache.get(bp, {})
    src = e.get("pmid_source", "")
    pmid  = bool(e.get("primary_pmid"))
    doi   = bool(e.get("primary_doi"))
    abstr = bool(e.get("abstract"))
    meth  = bool(e.get("methods_text"))

    if "Serper" in src:           stage = "Web search (Serper)"
    elif src == "manual_doi":     stage = "Manual DOI lookup"
    elif not pmid and not doi:    stage = None   # unresolved
    else:                         stage = "fetch_lit.py"

    if meth:              cov = "Methods + abstract"
    elif abstr and pmid:  cov = "Abstract (PMID)"
    elif abstr and doi:   cov = "Abstract (DOI only)"
    elif pmid or doi:     cov = "Title only"
    else:                 cov = None

    if stage is None:
        unresolved += 1
    else:
        stage_flows[stage][cov] += 1

stage_total   = {s: sum(v for v in stage_flows[s].values()) for s in STAGES}
# BPs entering each stage = total minus all previously resolved
remaining = total_bps
stage_in: dict[str, int] = {}
for s in STAGES:
    stage_in[s] = remaining
    remaining  -= stage_total[s]
# remaining == unresolved

# ── Node layout ───────────────────────────────────────────────────────────────
#
# Spine (y ≈ 0.04, flows left → right, gets thinner as BPs drop off):
#   Origin → fetch_lit → Serper → Manual → Unresolved
#
# Drop nodes (x just right of their stage, y below spine):
#   fetch_lit drops  → x=0.37
#   Serper drops     → x=0.65
#   Manual drops     → x=0.87

SPINE_Y  = 0.04

STAGE_X = {
    "fetch_lit.py":         0.20,
    "Web search (Serper)":  0.50,
    "Manual DOI lookup":    0.73,
}
DROP_X = {
    "fetch_lit.py":         0.37,
    "Web search (Serper)":  0.65,
    "Manual DOI lookup":    0.87,
}

# Drop y-positions (per stage, from top to bottom by coverage key order)
# Spacing chosen so nodes don't overlap given pad=14
DROP_Y: dict[str, dict[str, float]] = {
    "fetch_lit.py": {
        "Methods + abstract": 0.33,
        "Abstract (PMID)":    0.73,
        "Abstract (DOI only)":0.88,
        "Title only":         0.96,
    },
    "Web search (Serper)": {
        "Methods + abstract": 0.33,
        "Abstract (PMID)":    0.50,
        "Abstract (DOI only)":0.68,
        "Title only":         0.83,
    },
    "Manual DOI lookup": {
        "Methods + abstract": 0.33,
        "Abstract (PMID)":    0.44,
        "Abstract (DOI only)":0.59,
        "Title only":         0.72,
    },
}

STAGE_COLORS = {
    "fetch_lit.py":         "#2980b9",
    "Web search (Serper)":  "#e67e22",
    "Manual DOI lookup":    "#8e44ad",
}
COVERAGE_COLORS = {
    "Methods + abstract":   "#1a5e36",
    "Abstract (PMID)":      "#27ae60",
    "Abstract (DOI only)":  "#16a085",
    "Title only":           "#f39c12",
}
ORIGIN_COLOR    = "#2c3e50"
UNRESOLV_COLOR  = "#7f8c8d"

def rgba(hex_color: str, alpha: float = 0.40) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"

# ── Build nodes ───────────────────────────────────────────────────────────────

labels, colors, xs, ys = [], [], [], []

def add_node(label: str, color: str, x: float, y: float) -> int:
    idx = len(labels)
    labels.append(label)
    colors.append(color)
    xs.append(x)
    ys.append(y)
    return idx

ORIGIN_IDX  = add_node(f"<b>All BioProjects</b><br>{total_bps:,} BPs",
                        ORIGIN_COLOR, 0.01, SPINE_Y)

stage_idx: dict[str, int] = {}
for s in STAGES:
    n_in = stage_in[s]
    stage_idx[s] = add_node(
        f"<b>{s}</b><br>{n_in:,} BPs in",
        STAGE_COLORS[s], STAGE_X[s], SPINE_Y,
    )

UNRESOLV_IDX = add_node(
    f"<b>Unresolved</b><br>{unresolved:,} BPs",
    UNRESOLV_COLOR, 0.92, SPINE_Y,
)

# Drop nodes: only emit if count > 0
drop_idx: dict[str, dict[str, int]] = {s: {} for s in STAGES}
for s in STAGES:
    for cov in COVERAGE_KEYS:
        count = stage_flows[s].get(cov, 0)
        if count > 0:
            drop_idx[s][cov] = add_node(
                f"<b>{cov}</b><br>{count:,} BPs",
                COVERAGE_COLORS[cov],
                DROP_X[s],
                DROP_Y[s][cov],
            )

# ── Build links ───────────────────────────────────────────────────────────────

srcs, tgts, vals, link_colors = [], [], [], []

def add_link(src: int, tgt: int, val: int, color_hex: str) -> None:
    if val <= 0:
        return
    srcs.append(src)
    tgts.append(tgt)
    vals.append(val)
    link_colors.append(rgba(color_hex))

# Origin → first stage
add_link(ORIGIN_IDX, stage_idx["fetch_lit.py"], total_bps, ORIGIN_COLOR)

# Stage → drops + pass-forward to next stage
spine_pairs = [
    ("fetch_lit.py",        "Web search (Serper)"),
    ("Web search (Serper)", "Manual DOI lookup"),
    ("Manual DOI lookup",   None),
]

for s, next_s in spine_pairs:
    s_color = STAGE_COLORS[s]
    s_idx   = stage_idx[s]
    # Drops
    for cov in COVERAGE_KEYS:
        count = stage_flows[s].get(cov, 0)
        if count > 0:
            add_link(s_idx, drop_idx[s][cov], count, s_color)
    # Pass remainder to next stage or to Unresolved
    n_pass = stage_in[s] - stage_total[s]
    if next_s:
        add_link(s_idx, stage_idx[next_s], n_pass, s_color)
    else:
        add_link(s_idx, UNRESOLV_IDX, n_pass, UNRESOLV_COLOR)

# ── Figure ────────────────────────────────────────────────────────────────────

resolved = total_bps - unresolved

fig = go.Figure(go.Sankey(
    arrangement="fixed",
    node=dict(
        pad=14,
        thickness=20,
        line=dict(color="white", width=0.4),
        label=labels,
        color=colors,
        x=xs,
        y=ys,
    ),
    link=dict(
        source=srcs,
        target=tgts,
        value=vals,
        color=link_colors,
    ),
))

fig.update_layout(
    title=dict(
        text=(f"Literature resolution pipeline — {total_bps:,} BioProjects<br>"
              f"<sup>Resolved: {resolved:,} ({100*resolved/total_bps:.1f}%)   "
              f"Unresolved: {unresolved:,} ({100*unresolved/total_bps:.1f}%)</sup>"),
        font=dict(size=16),
        x=0.01,
    ),
    font=dict(size=12, family="Arial"),
    paper_bgcolor="white",
    width=1080,
    height=600,
    margin=dict(l=10, r=10, t=85, b=10),
)

fig.write_html(OUT_HTML)
print(f"Written {OUT_HTML}")

try:
    fig.write_image(OUT_PNG, scale=2)
    print(f"Written {OUT_PNG}")
except Exception as e:
    print(f"PNG skipped ({e})")
