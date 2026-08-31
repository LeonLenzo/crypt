#!/usr/bin/env python3
"""
metadata/figures/sample_funnel_v2.py — Sample funnel Sankey.

Reads: metadata/output/classify_metadata/data/samples.tsv  (single source)
Flow:  All BioSamples → Evidence → Tissue → Setting → Stress → Pathogens

Run from crypt/:  python metadata/figures/sample_funnel_v2.py
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

SAMPLES_TSV = Path("metadata/output/classify_metadata/data/samples.tsv")
OUT_HTML    = Path("metadata/output/figures/sankey/sample_funnel_v2.html")
OUT_SVG     = Path("metadata/output/figures/sankey/sample_funnel_v2.svg")
OUT_PNG     = Path("metadata/output/figures/sankey/sample_funnel_v2.png")

# ── Load & classify rows ───────────────────────────────────────────────────────

rows = list(csv.DictReader(open(SAMPLES_TSV), delimiter="\t"))

def ev(r):
    return "doi" if r.get("doi", "").strip() else "none"

def ti(r):
    t = r.get("llm_tissue", "")
    return t if t in ("aerial", "non-aerial") else "unclear"

def se(r):
    s = r.get("llm_study_setting", "")
    return s if s in ("field", "lab") else "unclear"

def st(r):
    s = r.get("llm_stress", "")
    return s if s in ("biotic", "abiotic", "none") else "unclear"

def coinf(r):
    return "multiple" if r.get("co_infection_flag", "") in ("multi_species", "multi_kingdom") else "single"

# cells[(ev, ti, se, st)] = {"n": int, "multiple": int}
cells = defaultdict(lambda: {"n": 0, "multiple": 0})
for r in rows:
    key = (ev(r), ti(r), se(r), st(r))
    cells[key]["n"] += 1
    if coinf(r) == "multiple":
        cells[key]["multiple"] += 1

ALL_EV  = ["doi", "none"]
ALL_TI  = ["aerial", "non-aerial", "unclear"]
ALL_SE  = ["field", "lab", "unclear"]
ALL_ST  = ["biotic", "abiotic", "none", "unclear"]

EV_ACT  = ["doi"]          # flows forward
TI_ACT  = ["aerial"]       # flows forward
SE_ACT  = ["field", "lab"] # flows forward
ST_ACT  = ["biotic"]       # flows forward (abiotic/none terminal)

def C(ev_l, ti_l, se_l, st_l, fld="n"):
    return sum(cells[(e,t,s,p)][fld]
               for e in ev_l for t in ti_l for s in se_l for p in st_l)

n_total   = len(rows)
ev_n      = {e: C([e], ALL_TI, ALL_SE, ALL_ST) for e in ALL_EV}
ti_n      = {t: C(EV_ACT, [t], ALL_SE, ALL_ST) for t in ALL_TI}
se_n      = {s: C(EV_ACT, TI_ACT, [s], ALL_ST) for s in ALL_SE}
str_n     = {(s,p): C(EV_ACT, TI_ACT, [s], [p]) for s in SE_ACT for p in ALL_ST}
pa_n      = {(s,c): C(EV_ACT, TI_ACT, [s], ST_ACT, fld=c)
             for s in SE_ACT for c in ("n","multiple")}

print(f"Total BioSamples: {n_total:,}")
print(f"Evidence — doi: {ev_n['doi']:,}  none: {ev_n['none']:,}")
print(f"Tissue   — aerial: {ti_n['aerial']:,}  non-aerial: {ti_n['non-aerial']:,}  unclear: {ti_n['unclear']:,}")
print(f"Setting  — field: {se_n['field']:,}  lab: {se_n['lab']:,}  unclear: {se_n['unclear']:,}")
for s in SE_ACT:
    bio = str_n[(s,"biotic")]
    abio = str_n[(s,"abiotic")]
    none = str_n[(s,"none")]
    unk  = str_n[(s,"unclear")]
    nc   = pa_n[(s,"multiple")]
    print(f"  {s}: biotic={bio:,}  abiotic={abio:,}  none={none:,}  unclear={unk:,}  → coinf={nc:,}/{bio:,}")

# ── Colors ─────────────────────────────────────────────────────────────────────

TERMINAL  = "#c8d0d4"
EV_COL    = "#2471a3"
TI_COL    = "#1e8449"
SET_COLS  = {"field": "#1e8449", "lab": "#2e4057"}
STR_COLS  = {
    ("field", "biotic"):  "#1a6b3c",
    ("field", "abiotic"): "#52be80",
    ("field", "none"):    "#a9dfbf",
    ("lab",   "biotic"):  "#1f3f5b",
    ("lab",   "abiotic"): "#5b8db8",
    ("lab",   "none"):    "#a9c4d8",
}
COINF_COL = "#e67e22"
GREY_NC   = "#c8d0d4"

def rgba(c, a=0.42):
    h = c.lstrip("#")
    return f"rgba({int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)},{a})"

# ── Build Sankey nodes & links ─────────────────────────────────────────────────

node_labels, node_colors, node_x, node_y = [], [], [], []
link_src, link_tgt, link_val, link_col, link_lbl = [], [], [], [], []

def add_node(label, color, x, y=0.5):
    i = len(node_labels)
    node_labels.append(label); node_colors.append(color)
    node_x.append(x);          node_y.append(y)
    return i

def add_link(src, tgt, val, col, lbl=""):
    if val <= 0: return
    link_src.append(src); link_tgt.append(tgt)
    link_val.append(val); link_col.append(col); link_lbl.append(lbl)

# Col 0
idx_all = add_node(f"All BioSamples  {n_total:,}", "#1a5276", 0.01, 0.44)

# Col 1: Evidence
idx_doi  = add_node(f"DOI  {ev_n['doi']:,}",         EV_COL,   0.18)
idx_none = add_node(f"No DOI  {ev_n['none']:,}",      TERMINAL, 0.18)
add_link(idx_all, idx_doi,  ev_n["doi"],  rgba(EV_COL,   0.5), "DOI")
add_link(idx_all, idx_none, ev_n["none"], rgba(TERMINAL, 0.3), "No DOI")

# Col 2: Tissue
idx_aerial  = add_node(f"Aerial  {ti_n['aerial']:,}",         TI_COL,   0.34)
idx_nonaer  = add_node(f"Non-aerial  {ti_n['non-aerial']:,}", TERMINAL, 0.34)
idx_ti_unk  = add_node(f"Unclear  {ti_n['unclear']:,}",       TERMINAL, 0.34)

n_doi_aerial  = C(["doi"], ["aerial"],    ALL_SE, ALL_ST)
n_doi_nonaer  = C(["doi"], ["non-aerial"],ALL_SE, ALL_ST)
n_doi_tiunk   = C(["doi"], ["unclear"],   ALL_SE, ALL_ST)
add_link(idx_doi, idx_aerial,  n_doi_aerial,  rgba(EV_COL, 0.5), f"aerial: {n_doi_aerial:,}")
add_link(idx_doi, idx_nonaer,  n_doi_nonaer,  rgba(TERMINAL,0.3), f"non-aerial: {n_doi_nonaer:,}")
add_link(idx_doi, idx_ti_unk,  n_doi_tiunk,   rgba(TERMINAL,0.3), f"unclear: {n_doi_tiunk:,}")

# Col 3: Setting
idx_field   = add_node(f"Field  {se_n['field']:,}", SET_COLS["field"], 0.50)
idx_lab     = add_node(f"Lab  {se_n['lab']:,}",     SET_COLS["lab"],   0.50)
idx_set_unk = add_node(f"Unclear  {se_n['unclear']:,}", TERMINAL,      0.50)

for s, idx in [("field", idx_field), ("lab", idx_lab)]:
    n = C(EV_ACT, TI_ACT, [s], ALL_ST)
    add_link(idx_aerial, idx, n, rgba(TI_COL, 0.5), f"aerial→{s}: {n:,}")
n_set_unk = C(EV_ACT, TI_ACT, ["unclear"], ALL_ST)
add_link(idx_aerial, idx_set_unk, n_set_unk, rgba(TERMINAL,0.3), f"unclear: {n_set_unk:,}")

# Col 4: Stress
str_idx     = {}
str_unk_idx = {}
STR_LABELS  = {"biotic": "Biotic", "abiotic": "Abiotic", "none": "None"}

for s, set_idx in [("field", idx_field), ("lab", idx_lab)]:
    for p in ["biotic", "abiotic", "none"]:
        n = str_n[(s, p)]
        if n == 0: continue
        col = STR_COLS[(s, p)]
        idx = add_node(f"{STR_LABELS[p]}  {n:,}", col, 0.67)
        str_idx[(s, p)] = idx
        add_link(set_idx, idx, n, rgba(col), f"{s}/{p}")
    n_unk = str_n[(s, "unclear")]
    if n_unk > 0:
        idx_unk = add_node(f"Unclear  {n_unk:,}", TERMINAL, 0.67)
        str_unk_idx[s] = idx_unk
        add_link(set_idx, idx_unk, n_unk, rgba(TERMINAL,0.3), f"{s}/unclear")

# Col 5: Pathogens (biotic only → single vs multiple)
pa_nc_idx  = {}
pa_c_idx   = {}
col5_annot = {}

for s in SE_ACT:
    n_c  = pa_n[(s, "multiple")]
    n_nc = pa_n[(s, "n")] - n_c
    if n_c + n_nc == 0: continue
    pa_nc_idx[s] = add_node("", GREY_NC,   0.84)
    pa_c_idx[s]  = add_node("", COINF_COL, 0.84)
    col5_annot[pa_nc_idx[s]] = f"Single  {n_nc:,}"
    col5_annot[pa_c_idx[s]]  = f"Multiple  {n_c:,}"

for s in SE_ACT:
    if (s, "biotic") not in str_idx: continue
    n_c  = pa_n[(s, "multiple")]
    n_nc = pa_n[(s, "n")] - n_c
    add_link(str_idx[(s,"biotic")], pa_nc_idx[s], n_nc, rgba(GREY_NC,   0.25), f"{s}/single")
    add_link(str_idx[(s,"biotic")], pa_c_idx[s],  n_c,  rgba(COINF_COL, 0.65), f"{s}/co-infected")

# ── Y positioning ──────────────────────────────────────────────────────────────

PAD   = 0.01
POWER = 0.55

def place_nodes(pairs: list[tuple[int, int]], pad=PAD):
    """Place nodes at proportional-sqrt y midpoints; return band (top, bot) list."""
    adj   = [max(s, 1) ** POWER for _, s in pairs]
    total = sum(adj)
    bands, y = [], pad
    for (idx, _), a in zip(pairs, adj):
        h   = (a / total) * (1 - 2 * pad)
        mid = y + h / 2
        node_y[idx] = max(pad, min(1 - pad, mid))
        bands.append((y, y + h))
        y += h
    return bands

# Col 1
ev_bands = place_nodes([(idx_doi, ev_n["doi"]), (idx_none, ev_n["none"])])

# Col 2
ti_bands = place_nodes([
    (idx_aerial, ti_n["aerial"]),
    (idx_nonaer, ti_n["non-aerial"]),
    (idx_ti_unk, ti_n["unclear"]),
])

# Col 3
set_bands = place_nodes([
    (idx_field,   se_n["field"]),
    (idx_lab,     se_n["lab"]),
    (idx_set_unk, se_n["unclear"]),
])
field_band = set_bands[0]
lab_band   = set_bands[1]

# Col 4: stress within each setting band
def place_stress_in_band(s, top, bot):
    pairs = []
    for p in ["biotic", "abiotic", "none"]:
        if (s, p) in str_idx:
            pairs.append((str_idx[(s, p)], str_n[(s, p)]))
    if s in str_unk_idx:
        pairs.append((str_unk_idx[s], str_n[(s, "unclear")]))
    if not pairs: return
    adj   = [max(sz, 1) ** POWER for _, sz in pairs]
    total = sum(adj)
    y = top
    for (idx, _), a in zip(pairs, adj):
        h   = (a / total) * (bot - top)
        mid = max(PAD, min(1 - PAD, y + h / 2))
        node_y[idx] = mid
        y += h

place_stress_in_band("field", field_band[0], field_band[1])
place_stress_in_band("lab",   lab_band[0],   lab_band[1])

# Col 5: co-infection within each setting band
def place_coinf_in_band(s, top, bot):
    if s not in pa_nc_idx: return
    n_c  = pa_n[(s, "multiple")]
    n_nc = pa_n[(s, "n")] - n_c
    tot  = n_c + n_nc
    if tot == 0: return
    span = bot - top
    node_y[pa_nc_idx[s]] = top + (n_nc / tot) * span / 2
    node_y[pa_c_idx[s]]  = top + (n_nc / tot) * span + (n_c / tot) * span / 2

place_coinf_in_band("field", field_band[0], field_band[1])
place_coinf_in_band("lab",   lab_band[0],   lab_band[1])

node_y[idx_all] = 0.44

print("\nNode index → y (for manual overrides if needed):")
for i, (lbl, y) in enumerate(zip(node_labels, node_y)):
    print(f"  [{i:2d}] y={y:.3f}  {lbl}")

# ── Figure ─────────────────────────────────────────────────────────────────────

W, H     = 1700, 1000
M_L, M_R = 20, 120
M_T, M_B = 110, 100

def s2px(sx): return M_L / W + sx * (W - M_L - M_R) / W

n_field_coinf = pa_n[("field","multiple")]
n_lab_coinf   = pa_n[("lab","multiple")]
n_field_bio   = pa_n[("field","n")]
n_lab_bio     = pa_n[("lab","n")]
pct = lambda n, d: f"{100*n/d:.1f}%" if d else "0%"
upstream = (
    "Upstream: <b>608,368</b> STAT-screened runs "
    "→ <b>10,995</b> gate-pass → <b>9,002</b> BioSamples (eukaryotic pathogens only).  "
    f"Co-infection: field <b>{n_field_coinf:,}/{n_field_bio:,}</b> ({pct(n_field_coinf,n_field_bio)})"
    f"  lab <b>{n_lab_coinf:,}/{n_lab_bio:,}</b> ({pct(n_lab_coinf,n_lab_bio)})"
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
        text    = "<b>Co-infection prevalence in field-collected plant RNA-seq</b><br>"
                  f"<span style='font-size:13px'>{upstream}</span>",
        font    = dict(size=20),
        x=0.0, xanchor="left",
    ),
    font          = dict(size=15, family="sans-serif"),
    width         = W,
    height        = H,
    margin        = dict(l=M_L, r=M_R, t=M_T, b=M_B),
    paper_bgcolor = "white",
)

COLUMN_TITLES   = ["BioSamples", "Evidence", "Tissue", "Setting", "Stress", "Pathogens"]
COLUMN_X_SANKEY = [0.01,          0.18,       0.34,    0.50,      0.67,     0.84]

for title, sx in zip(COLUMN_TITLES, COLUMN_X_SANKEY):
    fig.add_annotation(
        x=s2px(sx), y=0.03, xref="paper", yref="paper",
        text=f"<b>{title}</b>", showarrow=False,
        xanchor="center", yanchor="middle", font=dict(size=16),
    )

# ── SVG: add right-side labels for col5 nodes ─────────────────────────────────

def _svg_text(html, x, y, font_size=15, line_height=18):
    parts  = re.split(r'<br\s*/?>', html, flags=re.I)
    y0     = y - (len(parts) - 1) * line_height / 2
    tspans = []
    for i, part in enumerate(parts):
        bold   = bool(re.search(r'<b>', part, re.I))
        text   = re.sub(r'<[^>]+>', '', part).strip()
        weight = ' font-weight="bold"' if bold else ''
        tspans.append(
            f'<tspan x="{x:.1f}" y="{y0 + i*line_height:.1f}"{weight}>{text}</tspan>'
        )
    return (
        f'<text text-anchor="start" dominant-baseline="central" '
        f'font-family="sans-serif" font-size="{font_size}" fill="#2a2a2a">'
        + ''.join(tspans) + '</text>'
    )

def fix_col5_labels(svg_path):
    plot_w  = W - M_L - M_R
    label_x = M_L + 0.84 * plot_w + 18 + 5
    inserts = [
        _svg_text(html, label_x, M_T + node_y[nidx] * (H - M_T - M_B))
        for nidx, html in col5_annot.items()
    ]
    content = svg_path.read_text()
    svg_path.write_text(content.replace('</svg>', '\n'.join(inserts) + '\n</svg>'))

# ── Export ─────────────────────────────────────────────────────────────────────

OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
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
    print(f"PNG skipped: {e}")
