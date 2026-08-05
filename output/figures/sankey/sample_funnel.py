#!/usr/bin/env python3
"""
figure/sankey/sample_funnel.py
5-column funnel Sankey: All BioSamples → Evidence → Setting → Stress → Co-infection
Unclear nodes terminate at each stage — only classifiable BioSamples proceed.
Run from crypt/: python figure/sankey/sample_funnel.py
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

try:
    import plotly.graph_objects as go
except ImportError:
    sys.exit("pip install plotly kaleido")

OUT_HTML = Path("figure/sankey/sample_funnel.html")
OUT_PNG  = Path("figure/sankey/sample_funnel.png")

# ── Load data ──────────────────────────────────────────────────────────────────

llm = {}
with open("output/05_llm_classify/data/bioproject_llm.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        llm[r["BioProject"]] = r

rows = []
with open("output/04_filter_kw/data/biosample_kw.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)

# ── Classification ─────────────────────────────────────────────────────────────

bp_has_xml = {}
for r in rows:
    bp = r["BioProject"]
    if bp not in bp_has_xml:
        bp_has_xml[bp] = False
    if r.get("geo_loc_name") or r.get("tissue") or r.get("collection_date"):
        bp_has_xml[bp] = True

def evidence(r):
    if r.get("primary_pmid", "").strip():    return "publication"
    if bp_has_xml.get(r["BioProject"], False): return "metadata"
    return "unclear"

def setting(r):
    s = llm.get(r["BioProject"], {}).get("llm_study_setting", "")
    if s == "field": return "field"
    if s == "lab":   return "controlled"
    return "unclear"

def stress(r):
    t = llm.get(r["BioProject"], {}).get("llm_treatment", "unclear")
    if t in ("single", "coinf_experiment", "surveillance"): return "biotic"
    if t == "abiotic_stress": return "abiotic"
    if t == "host_study":     return "none"
    return "unclear"

# ── Cell counts ────────────────────────────────────────────────────────────────

# Active (classifiable) categories — "unclear" terminates at each stage
EV_ACTIVE   = ["publication", "metadata"]
SET_ACTIVE  = ["field", "controlled"]
STR_ACTIVE  = ["biotic", "abiotic", "none"]

cells = defaultdict(lambda: {"n": 0, "coinf": 0})
for r in rows:
    e, s, t = evidence(r), setting(r), stress(r)
    cells[(e, s, t)]["n"] += 1
    if r["co_infection_flag"] != "single":
        cells[(e, s, t)]["coinf"] += 1

n_total = len(rows)
pct = lambda n, d: f"{100*n/d:.1f}%" if d else "0%"

# Marginal counts (including unclear at each level)
ev_n    = {e: sum(cells[(e,s,t)]["n"] for s in SET_ACTIVE+["unclear"] for t in STR_ACTIVE+["unclear"])
           for e in EV_ACTIVE + ["unclear"]}
set_n   = {s: sum(cells[(e,s,t)]["n"] for e in EV_ACTIVE for t in STR_ACTIVE+["unclear"])
           for s in SET_ACTIVE + ["unclear"]}
# Counts reaching stress column (active ev + active setting)
str_n   = {(s,t): sum(cells[(e,s,t)]["n"] for e in EV_ACTIVE)
           for s in SET_ACTIVE for t in STR_ACTIVE + ["unclear"]}

n_ev_unclear  = ev_n["unclear"]
n_set_unclear = sum(
    sum(cells[(e,"unclear",t)]["n"] for t in STR_ACTIVE+["unclear"])
    for e in EV_ACTIVE
)
n_str_unclear = {s: str_n[(s,"unclear")] for s in SET_ACTIVE}

# Print funnel summary
print(f"\nEvidence:   {ev_n['publication']:,} publication  |  {ev_n['metadata']:,} metadata  |  {n_ev_unclear:,} no-evidence (terminal)")
print(f"Setting:    {set_n['field']:,} field  |  {set_n['controlled']:,} controlled  |  {n_set_unclear:,} unclear (terminal)")
for s in SET_ACTIVE:
    active = sum(str_n[(s,t)] for t in STR_ACTIVE)
    print(f"  {s}: {active:,} continue  |  {n_str_unclear[s]:,} unclear design (terminal)")
print()
print(f"{'Evidence':<14} {'Setting':<12} {'Stress':<10} {'n':>6}  {'coinf%':>7}")
print("-" * 58)
for e in EV_ACTIVE:
    for s in SET_ACTIVE:
        for t in STR_ACTIVE:
            c = cells[(e,s,t)]
            if c["n"] == 0: continue
            print(f"{e:<14} {s:<12} {t:<10} {c['n']:>6}  {100*c['coinf']/c['n']:>6.1f}%")

# ── Colors ─────────────────────────────────────────────────────────────────────

# 3-family palette: greens (field), blues (controlled), greys (terminal/unclear)
# Co-infection uses a single bold orange-red to stand out from both families.

TERMINAL_COL = "#c8d0d4"   # neutral grey — all dead-end nodes

# Evidence (neutral blues — pre-classification, no setting yet)
EV_COLS = {"publication": "#2471a3", "metadata": "#7fb3d3"}

# Setting
SET_COLS = {"field": "#1e8449", "controlled": "#2e4057"}

# Stress: lighter shades of the parent setting colour
STRESS_COLS = {
    ("field",      "biotic"):  "#1a6b3c",   # dark green
    ("field",      "abiotic"): "#52be80",   # medium green
    ("field",      "none"):    "#a9dfbf",   # light green
    ("controlled", "biotic"):  "#1f3f5b",   # dark blue
    ("controlled", "abiotic"): "#5b8db8",   # medium blue
    ("controlled", "none"):    "#a9c4d8",   # light blue
}

COINF_COL = "#e67e22"   # bold orange — co-infection signal, distinct from both families
GREY_NC   = "#c8d0d4"   # same as terminal — no co-infection is unremarkable

def rgba(c, a=0.42):
    h = c.lstrip("#")
    r_, g_, b_ = int(h[0:2],16), int(h[2:4],16), int(h[4:6],16)
    return f"rgba({r_},{g_},{b_},{a})"

# ── Build nodes and links ──────────────────────────────────────────────────────

node_labels, node_colors, node_x, node_y = [], [], [], []
link_src, link_tgt, link_val, link_col, link_lbl = [], [], [], [], []

def add_node(label, color, x, y=0.5):
    idx = len(node_labels)
    node_labels.append(label)
    node_colors.append(color)
    node_x.append(x)
    node_y.append(y)
    return idx

def add_link(src, tgt, val, col, lbl):
    if val <= 0: return
    link_src.append(src); link_tgt.append(tgt)
    link_val.append(val); link_col.append(col); link_lbl.append(lbl)

# Col 0
idx_all = add_node(f"All BioSamples<br><b>{n_total:,}</b>", "#1a5276", 0.01, 0.44)

# Col 1: Evidence — unclear terminates here
ev_idx = {}
for e in EV_ACTIVE:
    n = ev_n[e]
    lbl = "With publication" if e == "publication" else "Metadata only"
    idx = add_node(f"{lbl}<br><b>{n:,}</b> ({pct(n,n_total)})", EV_COLS[e], 0.22)
    ev_idx[e] = idx
    add_link(idx_all, idx, n, rgba(EV_COLS[e]), lbl)

idx_ev_unk = add_node(
    f"No evidence<br><b>{n_ev_unclear:,}</b> ({pct(n_ev_unclear,n_total)})<br><i>not classified</i>",
    TERMINAL_COL, 0.22)
add_link(idx_all, idx_ev_unk, n_ev_unclear, rgba(TERMINAL_COL, 0.3), "No evidence — terminal")

# Col 2: Setting — unclear terminates here
set_idx = {}
for s in SET_ACTIVE:
    n = set_n[s]
    lbl = "Field" if s == "field" else "Controlled"
    idx = add_node(f"{lbl}<br><b>{n:,}</b> ({pct(n,n_total)})", SET_COLS[s], 0.44)
    set_idx[s] = idx

idx_set_unk = add_node(
    f"Setting unclear<br><b>{n_set_unclear:,}</b> ({pct(n_set_unclear,n_total)})<br><i>not classified</i>",
    TERMINAL_COL, 0.44)

for e in EV_ACTIVE:
    for s in SET_ACTIVE:
        n = sum(cells[(e,s,t)]["n"] for t in STR_ACTIVE + ["unclear"])
        add_link(ev_idx[e], set_idx[s], n, rgba(EV_COLS[e]), f"{e} → {s}: {n:,}")
    n_unk = sum(cells[(e,"unclear",t)]["n"] for t in STR_ACTIVE+["unclear"])
    add_link(ev_idx[e], idx_set_unk, n_unk, rgba(TERMINAL_COL, 0.3), f"{e} → unclear: {n_unk:,}")

# Col 3: Stress — unclear terminates here (separate terminal per setting)
stress_idx    = {}
str_unk_idx   = {}
STR_LABELS    = {"biotic": "Biotic", "abiotic": "Abiotic", "none": "No stress"}

for s in SET_ACTIVE:
    for t in STR_ACTIVE:
        n = str_n[(s,t)]
        if n == 0: continue
        col = STRESS_COLS[(s,t)]
        set_lbl = "Field" if s == "field" else "Ctrl"
        idx = add_node(f"{STR_LABELS[t]}<br><b>{n:,}</b>", col, 0.66)
        stress_idx[(s,t)] = idx
        add_link(set_idx[s], idx, n, rgba(col), f"{s}/{t}: {n:,}")

    n_unk = n_str_unclear[s]
    if n_unk > 0:
        set_lbl = "Field" if s == "field" else "Controlled"
        idx_unk = add_node(
            f"Design unclear<br><b>{n_unk:,}</b><br><i>not classified</i>",
            TERMINAL_COL, 0.66)
        str_unk_idx[s] = idx_unk
        add_link(set_idx[s], idx_unk, n_unk, rgba(TERMINAL_COL, 0.3),
                 f"{s}/unclear design: {n_unk:,}")

# Col 4: Co-infection nodes per setting (nc + coinf × field/controlled = 4 nodes)
coinf_nc_idx = {}
coinf_c_idx  = {}

for s in SET_ACTIVE:
    active_n = sum(str_n[(s,t)] for t in STR_ACTIVE)
    n_c  = sum(sum(cells[(e,s,t)]["coinf"] for e in EV_ACTIVE) for t in STR_ACTIVE)
    n_nc = active_n - n_c
    lbl  = "Field" if s == "field" else "Controlled"
    coinf_nc_idx[s] = add_node(
        f"No co-infection<br><b>{n_nc:,}</b> ({pct(n_nc, active_n)})", GREY_NC, 0.92)
    coinf_c_idx[s]  = add_node(
        f"Co-infected<br><b>{n_c:,}</b> ({pct(n_c, active_n)})", COINF_COL, 0.92)

for s in SET_ACTIVE:
    for t in STR_ACTIVE:
        if (s,t) not in stress_idx: continue
        tidx  = stress_idx[(s,t)]
        n_c   = sum(cells[(e,s,t)]["coinf"] for e in EV_ACTIVE)
        n_tot = str_n[(s,t)]
        n_nc  = n_tot - n_c
        col   = STRESS_COLS[(s,t)]
        add_link(tidx, coinf_nc_idx[s], n_nc, rgba(GREY_NC, 0.25),
                 f"{s}/{t} no coinf: {n_nc:,}")
        add_link(tidx, coinf_c_idx[s],  n_c,  rgba(col, 0.65),
                 f"{s}/{t} co-infected: {n_c:,} ({pct(n_c,n_tot)})")

# ── Y positions ────────────────────────────────────────────────────────────────

PAD = 0.01

# At the setting column, the two active settings + unclear terminal fill the space.
# Field sits at top, controlled in middle, unclear terminal at bottom.
n_field      = set_n["field"]
n_controlled = set_n["controlled"]
n_set_total  = n_field + n_controlled + n_set_unclear  # = active ev total

def proportional_bands(sizes, pad=PAD):
    total = sum(sizes)
    bands, y = [], pad
    for sz in sizes:
        h = (sz / total) * (1 - 2*pad)
        bands.append((y, y+h))
        y += h
    return bands

# Evidence bands (proportional to publication, metadata, unclear)
ev_sizes  = [ev_n["publication"], ev_n["metadata"], n_ev_unclear]
ev_bands  = proportional_bands(ev_sizes)
node_y[ev_idx["publication"]] = sum(ev_bands[0]) / 2
node_y[ev_idx["metadata"]]    = sum(ev_bands[1]) / 2
node_y[idx_ev_unk]            = sum(ev_bands[2]) / 2
node_y[idx_all]               = 0.44

# Setting bands
set_sizes = [n_field, n_controlled, n_set_unclear]
set_bands = proportional_bands(set_sizes)
node_y[set_idx["field"]]      = sum(set_bands[0]) / 2
node_y[set_idx["controlled"]] = sum(set_bands[1]) / 2
node_y[idx_set_unk]           = sum(set_bands[2]) / 2

# Stress bands: within field band and controlled band; unclear appended at bottom of each
field_top, field_bot       = set_bands[0]
controlled_top, ctrl_bot   = set_bands[1]

def place_stress(s, band_top, band_bot):
    sizes = [str_n[(s,t)] for t in STR_ACTIVE if str_n[(s,t)] > 0] + [n_str_unclear[s]]
    keys  = [t for t in STR_ACTIVE if str_n[(s,t)] > 0] + ["unclear"]
    total = sum(sizes)
    y = band_top
    for key, sz in zip(keys, sizes):
        h = (sz / total) * (band_bot - band_top)
        mid = max(PAD, min(1-PAD, y + h/2))
        if key == "unclear":
            if s in str_unk_idx:
                node_y[str_unk_idx[s]] = mid
        else:
            if (s,key) in stress_idx:
                node_y[stress_idx[(s,key)]] = mid
        y += h

place_stress("field",      field_top, field_bot)
place_stress("controlled", controlled_top, ctrl_bot)

# Co-infection nodes: nc top, coinf bottom within each setting's band
def place_coinf(s, band_top, band_bot):
    active_n = sum(str_n[(s,t)] for t in STR_ACTIVE)
    n_c  = sum(sum(cells[(e,s,t)]["coinf"] for e in EV_ACTIVE) for t in STR_ACTIVE)
    n_nc = active_n - n_c
    span = band_bot - band_top
    node_y[coinf_nc_idx[s]] = band_top + (n_nc / active_n) * span / 2
    node_y[coinf_c_idx[s]]  = band_top + (n_nc / active_n) * span + (n_c / active_n) * span / 2

place_coinf("field",      field_top, field_bot)
place_coinf("controlled", controlled_top, ctrl_bot)

# ── Figure ─────────────────────────────────────────────────────────────────────

n_classified = sum(str_n[(s,t)] for s in SET_ACTIVE for t in STR_ACTIVE)
n_coinf_tot  = sum(sum(cells[(e,s,t)]["coinf"] for e in EV_ACTIVE)
                   for s in SET_ACTIVE for t in STR_ACTIVE)
upstream = (
    "Upstream: <b>608,368</b> STAT-screened runs "
    "→ <b>10,995</b> gate-pass → <b>9,002</b> BioSamples (eukaryotic pathogens only).  "
    f"<b>{n_classified:,}</b> fully classified; "
    f"<b>{n_coinf_tot:,}</b> co-infected ({pct(n_coinf_tot, n_classified)})"
)

fig = go.Figure(go.Sankey(
    arrangement = "snap",
    node = dict(
        label         = node_labels,
        color         = node_colors,
        x             = node_x,
        y             = node_y,
        pad           = 28,
        thickness     = 18,
        hovertemplate = "%{label}<extra></extra>",
    ),
    link = dict(
        source        = link_src,
        target        = link_tgt,
        value         = link_val,
        color         = link_col,
        label         = link_lbl,
        hovertemplate = "%{label}<extra></extra>",
    ),
))

fig.update_layout(
    title = dict(
        text = (
            "<b>Co-infection prevalence funnel — cryptic co-infection in plant RNA-seq</b><br>"
            f"<sup>{upstream}</sup>"
        ),
        font    = dict(size=14),
        x=0.0, xanchor="left",
    ),
    font          = dict(size=11, family="sans-serif"),
    width         = 2000,
    height        = 1100,
    margin        = dict(l=20, r=20, t=110, b=120),
    paper_bgcolor = "white",
)

fig.write_html(str(OUT_HTML))
print(f"\nWritten: {OUT_HTML}")
try:
    fig.write_image(str(OUT_PNG), scale=2)
    print(f"Written: {OUT_PNG}")
except Exception as e:
    print(f"PNG export skipped: {e}")
