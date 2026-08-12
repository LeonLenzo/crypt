#!/usr/bin/env python3
"""
6-column funnel Sankey:
  All BioSamples → Evidence → Tissue → Setting → Stress → Pathogens

Run from crypt/: python metadata/figures/sample_funnel.py
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import plotly.graph_objects as go
except ImportError:
    sys.exit("pip install plotly kaleido")

OUT_HTML = Path("metadata/output/figures/sankey/sample_funnel.html")
OUT_SVG  = Path("metadata/output/figures/sankey/sample_funnel.svg")
OUT_PNG  = Path("metadata/output/figures/sankey/sample_funnel.png")

# ── Y position overrides ───────────────────────────────────────────────────────
# After opening the HTML and dragging nodes to preferred positions, extract y
# values from the browser console and populate here to lock them permanently.
# Run in console after dragging:
#   JSON.stringify(Plotly.d3.select('.sankey').datum().nodes.map(n=>n.y0))
# Then set e.g. Y_OVERRIDES = {3: 0.62, 7: 0.15, ...}
Y_OVERRIDES: dict = {
    0:  0.387,  # All BioSamples
    1:  0.257,  # Publication
    2:  0.675,  # Metadata
    3:  0.836,  # None (no evidence)
    4:  0.195,  # Aerial
    5:  0.471,  # Non-aerial
    6:  0.670,  # Unclear (tissue)
    7:  0.043,  # Field
    8:  0.290,  # Controlled
    9:  0.473,  # Unclear (setting)
    10: 0.041,  # Field Biotic
    11: 0.099,  # Field Abiotic
    12: 0.136,  # Field None
    13: 0.281,  # Ctrl Biotic
    14: 0.422,  # Ctrl Abiotic
    15: 0.463,  # Ctrl None
    16: 0.503,  # Ctrl Unclear
    17: 0.115,  # Field Single
    18: 0.020,  # Field Multiple
    19: 0.366,  # Ctrl Single
    20: 0.185,  # Ctrl Multiple
}

# ── Load data ──────────────────────────────────────────────────────────────────

llm = {}
with open("metadata/output/llm_classify/data/bioproject_llm.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        llm[r["BioProject"]] = r

rows = []
with open("metadata/output/filter_kw/data/biosample_kw.tsv") as f:
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
    if r.get("primary_pmid", "").strip():      return "publication"
    if bp_has_xml.get(r["BioProject"], False): return "metadata"
    return "unclear"

def tissue_cat(r):
    cat = r.get("tissue_category", "")
    if cat in ("leaf", "aerial_other"):              return "aerial"
    if cat in ("root", "seed_fruit", "whole_plant"): return "non_aerial"
    return "unknown"

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

EV_ACTIVE  = ["publication", "metadata"]
TI_ACTIVE  = ["aerial"]
SET_ACTIVE = ["field", "controlled"]
STR_ACTIVE = ["biotic", "abiotic", "none"]

cells = defaultdict(lambda: {"n": 0, "coinf": 0})
for r in rows:
    e, ti, s, t = evidence(r), tissue_cat(r), setting(r), stress(r)
    cells[(e, ti, s, t)]["n"] += 1
    if r.get("co_infection_flag", "none") != "single":
        cells[(e, ti, s, t)]["coinf"] += 1

n_total = len(rows)
pct = lambda n, d: f"{100*n/d:.1f}%" if d else "0%"

# Evidence marginals (sum over ALL tissue categories)
ev_n = {
    e: sum(cells[(e, ti, s, t)]["n"]
           for ti in TI_ACTIVE + ["non_aerial", "unknown"]
           for s  in SET_ACTIVE + ["unclear"]
           for t  in STR_ACTIVE + ["unclear"])
    for e in EV_ACTIVE + ["unclear"]
}
n_ev_unclear = ev_n["unclear"]

# Tissue marginals (active evidence only)
ti_n = {
    ti: sum(cells[(e, ti, s, t)]["n"]
            for e in EV_ACTIVE
            for s in SET_ACTIVE + ["unclear"]
            for t in STR_ACTIVE + ["unclear"])
    for ti in TI_ACTIVE + ["non_aerial", "unknown"]
}
n_ti_nonaer  = ti_n["non_aerial"]
n_ti_unknown = ti_n["unknown"]

# Setting marginals (active evidence + aerial tissue)
set_n = {
    s: sum(cells[(e, ti, s, t)]["n"]
           for e  in EV_ACTIVE
           for ti in TI_ACTIVE
           for t  in STR_ACTIVE + ["unclear"])
    for s in SET_ACTIVE + ["unclear"]
}
n_set_unclear = set_n["unclear"]

# Stress marginals
str_n = {
    (s, t): sum(cells[(e, ti, s, t)]["n"]
                for e  in EV_ACTIVE
                for ti in TI_ACTIVE)
    for s in SET_ACTIVE for t in STR_ACTIVE + ["unclear"]
}
n_str_unclear = {s: str_n[(s, "unclear")] for s in SET_ACTIVE}

# Print summary
print(f"\nEvidence:  {ev_n['publication']:,} publication  |  "
      f"{ev_n['metadata']:,} metadata  |  {n_ev_unclear:,} none (terminal)")
print(f"Tissue:    {ti_n['aerial']:,} aerial  |  "
      f"{n_ti_nonaer:,} non-aerial (terminal)  |  {n_ti_unknown:,} unclear (terminal)")
print(f"Setting:   {set_n['field']:,} field  |  "
      f"{set_n['controlled']:,} controlled  |  {n_set_unclear:,} unclear (terminal)")
for s in SET_ACTIVE:
    n_c  = sum(cells[(e, ti, s, t)]["coinf"] for e in EV_ACTIVE for ti in TI_ACTIVE for t in STR_ACTIVE)
    n_nc = sum(str_n[(s, t)] for t in STR_ACTIVE) - n_c
    print(f"  {s}: {n_c:,} co-infected  |  {n_nc:,} no co-infection")

# ── Colors ─────────────────────────────────────────────────────────────────────

TERMINAL_COL = "#c8d0d4"
EV_COLS      = {"publication": "#2471a3", "metadata": "#7fb3d3"}
TI_COLS      = {"aerial": "#1e8449", "non_aerial": "#c8d0d4", "unknown": "#7f8c8d"}
SET_COLS     = {"field": "#1e8449", "controlled": "#2e4057"}
STRESS_COLS  = {
    ("field",      "biotic"):  "#1a6b3c",
    ("field",      "abiotic"): "#52be80",
    ("field",      "none"):    "#a9dfbf",
    ("controlled", "biotic"):  "#1f3f5b",
    ("controlled", "abiotic"): "#5b8db8",
    ("controlled", "none"):    "#a9c4d8",
}
COINF_COL = "#e67e22"   # orange — co-infected
GREY_NC   = "#c8d0d4"   # grey   — no co-infection

def rgba(c, a=0.42):
    h = c.lstrip("#")
    r_, g_, b_ = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
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
    if val <= 0:
        return
    link_src.append(src); link_tgt.append(tgt)
    link_val.append(val); link_col.append(col); link_lbl.append(lbl)

# Col 0: All BioSamples
idx_all = add_node(f"All BioSamples  {n_total:,}", "#1a5276", 0.01, 0.44)

# Col 1: Evidence
ev_idx = {}
for e in EV_ACTIVE:
    n   = ev_n[e]
    lbl = "Publication" if e == "publication" else "Metadata"
    idx = add_node(f"{lbl}  {n:,}", EV_COLS[e], 0.18)
    ev_idx[e] = idx
    add_link(idx_all, idx, n, rgba(EV_COLS[e]), lbl)

idx_ev_unk = add_node(f"None  {n_ev_unclear:,}", TERMINAL_COL, 0.18)
add_link(idx_all, idx_ev_unk, n_ev_unclear, rgba(TERMINAL_COL, 0.3), "No evidence")

# Col 2: Tissue
ti_idx = {}
for ti in TI_ACTIVE:
    n   = ti_n[ti]
    idx = add_node(f"Aerial  {n:,}", TI_COLS[ti], 0.34)
    ti_idx[ti] = idx

idx_ti_nonaer  = add_node(f"Non-aerial  {n_ti_nonaer:,}",  TI_COLS["non_aerial"], 0.34)
idx_ti_unknown = add_node(f"Unclear  {n_ti_unknown:,}",    TI_COLS["unknown"],    0.34)

for e in EV_ACTIVE:
    for ti in TI_ACTIVE:
        n = sum(cells[(e, ti, s, t)]["n"]
                for s in SET_ACTIVE + ["unclear"]
                for t in STR_ACTIVE + ["unclear"])
        add_link(ev_idx[e], ti_idx[ti], n, rgba(EV_COLS[e]), f"{e}→aerial: {n:,}")
    n_na = sum(cells[(e, "non_aerial", s, t)]["n"]
               for s in SET_ACTIVE + ["unclear"]
               for t in STR_ACTIVE + ["unclear"])
    add_link(ev_idx[e], idx_ti_nonaer, n_na, rgba(TERMINAL_COL, 0.3),
             f"{e}→non-aerial: {n_na:,}")
    n_unk = sum(cells[(e, "unknown", s, t)]["n"]
                for s in SET_ACTIVE + ["unclear"]
                for t in STR_ACTIVE + ["unclear"])
    add_link(ev_idx[e], idx_ti_unknown, n_unk, rgba(TI_COLS["unknown"], 0.3),
             f"{e}→unclear: {n_unk:,}")

# Col 3: Setting
set_idx = {}
for s in SET_ACTIVE:
    n   = set_n[s]
    lbl = "Field" if s == "field" else "Controlled"
    idx = add_node(f"{lbl}  {n:,}", SET_COLS[s], 0.50)
    set_idx[s] = idx

n_ti_active = sum(ti_n[ti] for ti in TI_ACTIVE)
idx_set_unk = add_node(f"Unclear  {n_set_unclear:,}", TERMINAL_COL, 0.50)

for ti in TI_ACTIVE:
    for s in SET_ACTIVE:
        n = sum(cells[(e, ti, s, t)]["n"]
                for e in EV_ACTIVE
                for t in STR_ACTIVE + ["unclear"])
        add_link(ti_idx[ti], set_idx[s], n, rgba(TI_COLS[ti]), f"aerial→{s}: {n:,}")
    n_unclear = sum(cells[(e, ti, "unclear", t)]["n"]
                    for e in EV_ACTIVE
                    for t in STR_ACTIVE + ["unclear"])
    add_link(ti_idx[ti], idx_set_unk, n_unclear, rgba(TERMINAL_COL, 0.3),
             f"aerial→unclear: {n_unclear:,}")

# Col 4: Stress
stress_idx  = {}
str_unk_idx = {}
STR_LABELS  = {"biotic": "Biotic", "abiotic": "Abiotic", "none": "None"}

for s in SET_ACTIVE:
    for t in STR_ACTIVE:
        n = str_n[(s, t)]
        if n == 0:
            continue
        col = STRESS_COLS[(s, t)]
        idx = add_node(f"{STR_LABELS[t]}  {n:,}", col, 0.67)
        stress_idx[(s, t)] = idx
        add_link(set_idx[s], idx, n, rgba(col), f"{s}/{t}: {n:,}")

    n_unk = n_str_unclear[s]
    if n_unk > 0:
        idx_unk = add_node(f"Unclear  {n_unk:,}", TERMINAL_COL, 0.67)
        str_unk_idx[s] = idx_unk
        add_link(set_idx[s], idx_unk, n_unk, rgba(TERMINAL_COL, 0.3),
                 f"{s}/unclear: {n_unk:,}")

# Col 5: Pathogens — blank labels; right-side text via annotations
coinf_nc_idx = {}
coinf_c_idx  = {}
col5_annot   = {}   # node_idx → label text

for s in SET_ACTIVE:
    active_n = sum(str_n[(s, t)] for t in STR_ACTIVE)
    n_c  = sum(cells[(e, ti, s, t)]["coinf"] for e in EV_ACTIVE for ti in TI_ACTIVE for t in STR_ACTIVE)
    n_nc = active_n - n_c
    coinf_nc_idx[s] = add_node("", GREY_NC,   0.84)
    coinf_c_idx[s]  = add_node("", COINF_COL, 0.84)
    col5_annot[coinf_nc_idx[s]] = f"Single  {n_nc:,}"
    col5_annot[coinf_c_idx[s]]  = f"Multiple  {n_c:,}"

for s in SET_ACTIVE:
    for t in STR_ACTIVE:
        if (s, t) not in stress_idx:
            continue
        tidx = stress_idx[(s, t)]
        n_c   = sum(cells[(e, ti, s, t)]["coinf"] for e in EV_ACTIVE for ti in TI_ACTIVE)
        n_tot = str_n[(s, t)]
        n_nc  = n_tot - n_c
        col   = STRESS_COLS[(s, t)]
        add_link(tidx, coinf_nc_idx[s], n_nc, rgba(GREY_NC,   0.25), f"{s}/{t} no coinf: {n_nc:,}")
        add_link(tidx, coinf_c_idx[s],  n_c,  rgba(COINF_COL, 0.65), f"{s}/{t} co-infected: {n_c:,}")

# ── Y positions ────────────────────────────────────────────────────────────────

PAD   = 0.01
POWER = 0.55

def proportional_bands(sizes, pad=PAD):
    adjusted = [max(s, 1) ** POWER for s in sizes]
    total    = sum(adjusted)
    bands, y = [], pad
    for adj in adjusted:
        h = (adj / total) * (1 - 2 * pad)
        bands.append((y, y + h))
        y += h
    return bands

# Col 1: Evidence
ev_sizes = [ev_n["publication"], ev_n["metadata"], n_ev_unclear]
ev_bands = proportional_bands(ev_sizes)
node_y[ev_idx["publication"]] = sum(ev_bands[0]) / 2
node_y[ev_idx["metadata"]]    = sum(ev_bands[1]) / 2
node_y[idx_ev_unk]            = sum(ev_bands[2]) / 2
node_y[idx_all]               = 0.44

# Col 2: Tissue
ti_sizes = [ti_n["aerial"], n_ti_nonaer, n_ti_unknown]
ti_bands = proportional_bands(ti_sizes)
node_y[ti_idx["aerial"]] = sum(ti_bands[0]) / 2
node_y[idx_ti_nonaer]    = sum(ti_bands[1]) / 2
node_y[idx_ti_unknown]   = sum(ti_bands[2]) / 2

# Col 3: Setting
set_sizes = [set_n["field"], set_n["controlled"], n_set_unclear]
set_bands = proportional_bands(set_sizes)
node_y[set_idx["field"]]      = sum(set_bands[0]) / 2
node_y[set_idx["controlled"]] = sum(set_bands[1]) / 2
node_y[idx_set_unk]           = sum(set_bands[2]) / 2

field_top, field_bot = set_bands[0]
ctrl_top,  ctrl_bot  = set_bands[1]

# Col 4: Stress (packed within field / controlled bands)
def place_stress(s, band_top, band_bot):
    sizes = [str_n[(s, t)] for t in STR_ACTIVE if str_n[(s, t)] > 0] + [n_str_unclear[s]]
    keys  = [t for t in STR_ACTIVE if str_n[(s, t)] > 0] + ["unclear"]
    total = sum(sizes)
    y     = band_top
    for key, sz in zip(keys, sizes):
        h   = (sz / total) * (band_bot - band_top)
        mid = max(PAD, min(1 - PAD, y + h / 2))
        if key == "unclear":
            if s in str_unk_idx:
                node_y[str_unk_idx[s]] = mid
        else:
            if (s, key) in stress_idx:
                node_y[stress_idx[(s, key)]] = mid
        y += h

place_stress("field",      field_top, field_bot)
place_stress("controlled", ctrl_top,  ctrl_bot)

# Col 5: Co-infection (2 nodes per setting, packed within setting bands)
def place_coinf(s, band_top, band_bot):
    active_n = sum(str_n[(s, t)] for t in STR_ACTIVE)
    n_c  = sum(cells[(e, ti, s, t)]["coinf"] for e in EV_ACTIVE for ti in TI_ACTIVE for t in STR_ACTIVE)
    n_nc = active_n - n_c
    span = band_bot - band_top
    node_y[coinf_nc_idx[s]] = band_top + (n_nc / active_n) * span / 2
    node_y[coinf_c_idx[s]]  = band_top + (n_nc / active_n) * span + (n_c / active_n) * span / 2

place_coinf("field",      field_top, field_bot)
place_coinf("controlled", ctrl_top,  ctrl_bot)

# Apply overrides last
for idx, y in Y_OVERRIDES.items():
    node_y[idx] = y

# Print node index table so Y_OVERRIDES can be populated after drag-and-drop
print("\nNode index → label (for Y_OVERRIDES):")
for i, (lbl, y) in enumerate(zip(node_labels, node_y)):
    clean = lbl.replace("<br>", " ").replace("<b>", "").replace("</b>", "")
    print(f"  [{i:2d}] y={y:.3f}  {clean}")

# ── Figure layout constants ────────────────────────────────────────────────────

W, H         = 1700, 1000
M_L, M_R     = 20, 120
M_T, M_B     = 110, 100

def s2px(sx):
    """Sankey x (0–1 within plot area) → paper x fraction (0–1 across full figure)."""
    return M_L / W + sx * (W - M_L - M_R) / W

def s2py(sy):
    """Sankey y (0=top, 1=bottom) → paper y fraction (0=bottom, 1=top)."""
    return 1.0 - (M_T / H + sy * (H - M_T - M_B) / H)

# ── Build figure ───────────────────────────────────────────────────────────────

n_classified = sum(str_n[(s, t)] for s in SET_ACTIVE for t in STR_ACTIVE)
n_coinf_tot  = sum(
    sum(cells[(e, ti, s, t)]["coinf"] for e in EV_ACTIVE for ti in TI_ACTIVE)
    for s in SET_ACTIVE for t in STR_ACTIVE
)
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
        pad           = 24,
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
        text    = "<b>Co-infection prevalence in field-collected plant RNA-seq</b>",
        font    = dict(size=20),
        x=0.0, xanchor="left",
    ),
    font          = dict(size=15, family="sans-serif"),
    width         = W,
    height        = H,
    margin        = dict(l=M_L, r=M_R, t=M_T, b=M_B),
    paper_bgcolor = "white",
)

# ── Column header annotations (bottom axis) ────────────────────────────────────

COLUMN_TITLES   = ["BioSamples", "Evidence", "Tissue", "Setting", "Stress", "Pathogens"]
COLUMN_X_SANKEY = [0.01,          0.18,       0.34,    0.50,      0.67,     0.84]

for title, sx in zip(COLUMN_TITLES, COLUMN_X_SANKEY):
    fig.add_annotation(
        x         = s2px(sx),
        y         = 0.03,       # in bottom margin (below plot area)
        xref      = "paper",
        yref      = "paper",
        text      = f"<b>{title}</b>",
        showarrow = False,
        xanchor   = "center",
        yanchor   = "middle",
        font      = dict(size=16),
    )

# ── SVG post-processing for col5 right-side labels ────────────────────────────
# Plotly hard-codes label flip for right-half nodes; we inject <text> elements
# at the correct absolute pixel positions after export instead.

def _svg_text(html, x, y, font_size=15, line_height=18):
    """Convert 'Line1<br><b>Line2</b>' to an SVG <text> element string."""
    parts = re.split(r'<br\s*/?>', html, flags=re.I)
    y0    = y - (len(parts) - 1) * line_height / 2   # centre block vertically
    tspans = []
    for i, part in enumerate(parts):
        bold   = bool(re.search(r'<b>', part, re.I))
        text   = re.sub(r'<[^>]+>', '', part).strip()
        weight = ' font-weight="bold"' if bold else ''
        tspans.append(
            f'<tspan x="{x:.1f}" y="{y0 + i * line_height:.1f}"{weight}>{text}</tspan>'
        )
    return (
        f'<text text-anchor="start" dominant-baseline="central" '
        f'font-family="sans-serif" font-size="{font_size}" fill="#2a2a2a">'
        + ''.join(tspans) + '</text>'
    )

def fix_col5_labels(svg_path):
    plot_w  = W - M_L - M_R
    plot_h  = H - M_T - M_B
    label_x = M_L + 0.84 * plot_w + 18 + 5   # node right edge + gap

    inserts = [
        _svg_text(html, label_x, M_T + node_y[nidx] * plot_h)
        for nidx, html in col5_annot.items()
    ]
    content = svg_path.read_text()
    svg_path.write_text(content.replace('</svg>', '\n'.join(inserts) + '\n</svg>'))

# ── Export ─────────────────────────────────────────────────────────────────────

fig.write_html(str(OUT_HTML))
print(f"\nWritten: {OUT_HTML}")

fig.write_image(str(OUT_SVG))
fix_col5_labels(OUT_SVG)
print(f"Written: {OUT_SVG}")

try:
    import cairosvg
    cairosvg.svg2png(url=str(OUT_SVG.resolve()), write_to=str(OUT_PNG), scale=2)
    print(f"Written: {OUT_PNG}")
except ImportError:
    print("PNG: pip install cairosvg  (SVG is ready)")
except Exception as e:
    print(f"PNG export skipped: {e}")
