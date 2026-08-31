#!/usr/bin/env python3
"""
metadata/figures/sample_funnel_v3.py — Sample funnel Sankey, built on meta_classify.py
output (full manuscript text classification), superseding sample_funnel_v2.py.

Reads: metadata/output/meta_classify/data/samples.tsv
Flow:  Full text retrieved -> Tissue -> Setting -> Stress -> Cryptic co-infection

Structural changes from v2 (see memory/research_hypothesis.md for full rationale):
  - No Evidence/DOI stage — meta_classify.py's entry gate already means every row
    here has a full manuscript, so the whole samples.tsv IS the entry population.
  - Setting: field/greenhouse both continue (wild pathogen exposure possible in
    both); growth_chamber/detached_leaf_assay/in_vitro/unclear are terminal
    (each shown separately, not merged — in_vitro specifically is worth seeing,
    since the STAT gate should keep these out but it's worth knowing if any slip
    through).
  - Stress: ALL of biotic/abiotic/none continue now (only unclear drops) — a
    single detected pathogen in an abiotic/none study is still cryptic, since
    nothing was expected at all. See is_cryptic() below for the exact logic.
  - Terminal stage renamed Pathogens -> Cryptic co-infection, and the
    single-vs-multiple STAT split is replaced by a stress-aware cryptic
    determination (not just "more than one pathogen detected"):
      biotic          -> cryptic if STAT found something beyond the named
                         pathogen(s) (pathogen_match_status in
                         {partial_match_plus_undeclared, no_match_stat_found_different,
                         stat_only_no_named})
      abiotic / none  -> cryptic if STAT found ANY pathogen at all (n_pathogens > 0)
                         — nothing was expected, so even one hit is novel

Run from crypt/:  python metadata/figures/sample_funnel_v3.py
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

SAMPLES_TSV = Path("metadata/output/meta_classify/data/samples.tsv")
OUT_HTML    = Path("metadata/output/figures/sankey/sample_funnel_v3.html")
OUT_SVG     = Path("metadata/output/figures/sankey/sample_funnel_v3.svg")
OUT_PNG     = Path("metadata/output/figures/sankey/sample_funnel_v3.png")

# ── Load & classify rows ───────────────────────────────────────────────────────

rows = list(csv.DictReader(open(SAMPLES_TSV), delimiter="\t"))

_CRYPTIC_PMS = {"partial_match_plus_undeclared", "no_match_stat_found_different",
                "stat_only_no_named"}


def ti(r):
    t = r.get("llm_tissue", "")
    return t if t in ("aerial", "non-aerial") else "unclear"


def se(r):
    s = r.get("llm_study_setting", "")
    return s if s in ("field", "greenhouse", "growth_chamber",
                      "detached_leaf_assay", "in_vitro") else "unclear"


def st(r):
    s = r.get("llm_stress", "")
    return s if s in ("biotic", "abiotic", "none") else "unclear"


def is_cryptic(r) -> bool:
    stress = st(r)
    if stress == "biotic":
        return r.get("pathogen_match_status", "") in _CRYPTIC_PMS
    if stress in ("abiotic", "none"):
        try:
            return int(r.get("n_pathogens") or 0) > 0
        except ValueError:
            return False
    return False


# cells[(ti, se, st)] = {"n": int, "cryptic": int}
cells = defaultdict(lambda: {"n": 0, "cryptic": 0})
for r in rows:
    key = (ti(r), se(r), st(r))
    cells[key]["n"] += 1
    if is_cryptic(r):
        cells[key]["cryptic"] += 1

ALL_TI = ["aerial", "non-aerial", "unclear"]
ALL_SE = ["field", "greenhouse", "growth_chamber", "detached_leaf_assay", "in_vitro", "unclear"]
ALL_ST = ["biotic", "abiotic", "none", "unclear"]

TI_ACT = ["aerial"]
ST_ACT = ["biotic", "abiotic", "none"]   # all continue now — only unclear drops

# Setting GROUPS: growth_chamber/detached_leaf_assay/in_vitro merge into one
# "Controlled" node and run the exact same Stress -> Co-infection pipeline as
# field/greenhouse (the side investigation mentioned in memory/research_hypothesis.md
# — is lab co-infection real contamination or STAT resolution error? — now lives
# in the same figure instead of a separate one-off check). "unclear" stays a
# pure dead end — we don't know its setting at all, nothing to split further.
GROUP_MEMBERS = {
    "field":      ["field"],
    "greenhouse": ["greenhouse"],
    "controlled": ["growth_chamber", "detached_leaf_assay", "in_vitro"],
    "unclear":    ["unclear"],
}
SE_GROUPS     = ["field", "greenhouse", "controlled", "unclear"]
SE_GROUPS_ACT = ["field", "greenhouse", "controlled"]


def C(ti_l, se_l, st_l, fld="n"):
    return sum(cells[(t, s, p)][fld] for t in ti_l for s in se_l for p in st_l)


n_total = len(rows)
ti_n    = {t: C([t], ALL_SE, ALL_ST) for t in ALL_TI}
se_n    = {g: C(TI_ACT, GROUP_MEMBERS[g], ALL_ST) for g in SE_GROUPS}
str_n   = {(g, p): C(TI_ACT, GROUP_MEMBERS[g], [p]) for g in SE_GROUPS_ACT for p in ALL_ST}
cr_n    = {(g, p, c): C(TI_ACT, GROUP_MEMBERS[g], [p], fld=c)
           for g in SE_GROUPS_ACT for p in ST_ACT for c in ("n", "cryptic")}

print(f"Total BioSamples (full text already gated): {n_total:,}")
print(f"Tissue   — aerial: {ti_n['aerial']:,}  non-aerial: {ti_n['non-aerial']:,}  unclear: {ti_n['unclear']:,}")
print(f"Setting  — field: {se_n['field']:,}  greenhouse: {se_n['greenhouse']:,}  "
      f"controlled: {se_n['controlled']:,}  unclear: {se_n['unclear']:,}")
for g in SE_GROUPS_ACT:
    print(f"  {g}:")
    for p in ST_ACT:
        n  = str_n[(g, p)]
        cr = cr_n[(g, p, "cryptic")]
        print(f"    {p:<8} n={n:,}  cryptic={cr:,} ({100*cr/max(n,1):.1f}%)")

# ── Colours ─────────────────────────────────────────────────────────────────────

TERMINAL = "#c8d0d4"
TI_COL   = "#1e8449"
SET_COLS = {"field": "#1e8449", "greenhouse": "#2e7d32", "controlled": "#34495e"}
STR_COLS = {
    ("field", "biotic"):  "#1a6b3c", ("field", "abiotic"):  "#52be80", ("field", "none"):  "#a9dfbf",
    ("greenhouse", "biotic"): "#1f3f5b", ("greenhouse", "abiotic"): "#5b8db8", ("greenhouse", "none"): "#a9c4d8",
    ("controlled", "biotic"): "#2c3e50", ("controlled", "abiotic"): "#7f9cb3", ("controlled", "none"): "#c3d2dd",
}
CRYPTIC_COL = "#e67e22"
GREY_NC     = "#c8d0d4"


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
    if val <= 0:
        return
    link_src.append(src); link_tgt.append(tgt)
    link_val.append(val); link_col.append(col); link_lbl.append(lbl)


# Col 0 — entry (no Evidence stage; meta_classify.py's full-text gate already applied)
idx_all = add_node(f"Full text retrieved  {n_total:,}", "#1a5276", 0.01, 0.44)

# Col 1: Tissue
idx_aerial = add_node(f"Aerial  {ti_n['aerial']:,}",         TI_COL,   0.18)
idx_nonaer = add_node(f"Non-aerial  {ti_n['non-aerial']:,}", TERMINAL, 0.18)
idx_ti_unk = add_node(f"Unclear  {ti_n['unclear']:,}",       TERMINAL, 0.18)
add_link(idx_all, idx_aerial, ti_n["aerial"],     rgba(TI_COL, 0.5), f"aerial: {ti_n['aerial']:,}")
add_link(idx_all, idx_nonaer, ti_n["non-aerial"], rgba(TERMINAL, 0.3), f"non-aerial: {ti_n['non-aerial']:,}")
add_link(idx_all, idx_ti_unk, ti_n["unclear"],    rgba(TERMINAL, 0.3), f"unclear: {ti_n['unclear']:,}")

# Col 2: Setting — field/greenhouse/controlled (merged growth_chamber +
# detached_leaf_assay + in_vitro) all continue; unclear is terminal
idx_field = add_node(f"Field  {se_n['field']:,}",           SET_COLS["field"],      0.32)
idx_green = add_node(f"Greenhouse  {se_n['greenhouse']:,}", SET_COLS["greenhouse"], 0.32)
idx_ctrl  = add_node(f"Controlled  {se_n['controlled']:,}", SET_COLS["controlled"], 0.32)
idx_se_unk = add_node(f"Unclear  {se_n['unclear']:,}",      TERMINAL, 0.32)

SET_NODE = {"field": idx_field, "greenhouse": idx_green, "controlled": idx_ctrl}
for g in SE_GROUPS:
    idx = SET_NODE.get(g, idx_se_unk)
    n   = se_n[g]
    col = SET_COLS.get(g, TERMINAL)
    add_link(idx_aerial, idx, n, rgba(col, 0.5 if g in SE_GROUPS_ACT else 0.3), f"{g}: {n:,}")

# Col 3: Stress (within each active setting group)
str_idx: dict = {}
STR_LABELS = {"biotic": "Biotic", "abiotic": "Abiotic", "none": "None"}

for g in SE_GROUPS_ACT:
    for p in ST_ACT:
        n = str_n[(g, p)]
        if n == 0:
            continue
        col = STR_COLS[(g, p)]
        idx = add_node(f"{STR_LABELS[p]}  {n:,}", col, 0.55)
        str_idx[(g, p)] = idx
        add_link(SET_NODE[g], idx, n, rgba(col), f"{g}/{p}")

# Col 4: Co-infection — ONE "Co-infected" and ONE "Not co-infected" node per
# setting group, fed by all 3 stress branches (not split by stress type
# anymore). abiotic/none are structurally ~100% cryptic in this dataset —
# runs.tsv is already gated on STAT detecting a pathogen at all, so
# n_pathogens>0 is true by construction for every row here — so they only ever
# contribute to "Co-infected", never "Not co-infected" (that split only exists
# within the biotic branch, which is the real discriminating statistic).
coinf_idx: dict = {}   # g -> Co-infected node idx
notco_idx: dict = {}   # g -> Not co-infected node idx
coinf_totals: dict = {}
notco_totals: dict = {}

for g in SE_GROUPS_ACT:
    n_cr = cr_n[(g, "biotic", "cryptic")] if (g, "biotic") in str_idx else 0
    n_nc = (str_n[(g, "biotic")] - n_cr) if (g, "biotic") in str_idx else 0
    other_total = sum(str_n[(g, p)] for p in ("abiotic", "none") if (g, p) in str_idx)
    total_coinf = n_cr + other_total
    total_notco = n_nc
    coinf_totals[g] = total_coinf
    notco_totals[g] = total_notco

    if total_coinf > 0:
        coinf_idx[g] = add_node("", CRYPTIC_COL, 0.86)
        if n_cr > 0:
            add_link(str_idx[(g, "biotic")], coinf_idx[g], n_cr, rgba(CRYPTIC_COL, 0.65), f"{g}/biotic/coinfected")
        for p in ("abiotic", "none"):
            if (g, p) in str_idx and str_n[(g, p)] > 0:
                add_link(str_idx[(g, p)], coinf_idx[g], str_n[(g, p)], rgba(CRYPTIC_COL, 0.5), f"{g}/{p}/coinfected")

    if total_notco > 0:
        notco_idx[g] = add_node("", GREY_NC, 0.86)
        add_link(str_idx[(g, "biotic")], notco_idx[g], n_nc, rgba(GREY_NC, 0.25), f"{g}/biotic/not-coinfected")

col4_annot: dict = {}
for g, idx in coinf_idx.items():
    col4_annot[idx] = f"Co-infected  {coinf_totals[g]:,}"
for g, idx in notco_idx.items():
    col4_annot[idx] = f"Not co-infected  {notco_totals[g]:,}"

# ── Y positioning (proportional-power placement, same approach as v2) ─────────

PAD   = 0.01
POWER = 0.55


def place_nodes(pairs: list, pad=PAD):
    adj   = [max(sz, 1) ** POWER for _, sz in pairs]
    total = sum(adj)
    bands, y = [], pad
    for (idx, _), a in zip(pairs, adj):
        h   = (a / total) * (1 - 2 * pad)
        mid = y + h / 2
        node_y[idx] = max(pad, min(1 - pad, mid))
        bands.append((y, y + h))
        y += h
    return bands


ti_bands = place_nodes([
    (idx_aerial, ti_n["aerial"]), (idx_nonaer, ti_n["non-aerial"]), (idx_ti_unk, ti_n["unclear"]),
])

se_bands = place_nodes([
    (idx_field, se_n["field"]), (idx_green, se_n["greenhouse"]),
    (idx_ctrl, se_n["controlled"]), (idx_se_unk, se_n["unclear"]),
])
BAND_OF = {"field": se_bands[0], "greenhouse": se_bands[1], "controlled": se_bands[2]}


def place_stress_in_band(g, top, bot):
    pairs = [(str_idx[(g, p)], str_n[(g, p)]) for p in ST_ACT if (g, p) in str_idx]
    if not pairs:
        return
    adj   = [max(sz, 1) ** POWER for _, sz in pairs]
    total = sum(adj)
    y = top
    for (idx, _), a in zip(pairs, adj):
        h   = (a / total) * (bot - top)
        node_y[idx] = max(PAD, min(1 - PAD, y + h / 2))
        y += h


def place_coinf_in_band(g, top, bot):
    total_coinf, total_notco = coinf_totals.get(g, 0), notco_totals.get(g, 0)
    tot = total_coinf + total_notco
    if tot == 0:
        return
    span = bot - top
    if g in notco_idx:
        node_y[notco_idx[g]] = top + (total_notco / tot) * span / 2
    if g in coinf_idx:
        node_y[coinf_idx[g]] = top + (total_notco / tot) * span + (total_coinf / tot) * span / 2


for g in SE_GROUPS_ACT:
    top, bot = BAND_OF[g]
    place_stress_in_band(g, top, bot)
    place_coinf_in_band(g, top, bot)

node_y[idx_all] = 0.44

print("\nNode index → y (for manual overrides if needed):")
for i, (lbl, y) in enumerate(zip(node_labels, node_y)):
    print(f"  [{i:2d}] y={y:.3f}  {lbl}")

# ── Figure ─────────────────────────────────────────────────────────────────────

W, H     = 1700, 1000
M_L, M_R = 20, 140
M_T, M_B = 110, 60


def s2px(sx):
    return M_L / W + sx * (W - M_L - M_R) / W


pct = lambda n, d: f"{100*n/d:.1f}%" if d else "0%"
upstream = (
    "Upstream: meta_search.py + meta_text.py full-text pipeline "
    f"→ meta_classify.py ({n_total:,} full-text-gated BioSamples).  "
    f"Co-infection: field <b>{coinf_totals['field']:,}/{se_n['field']:,}</b> ({pct(coinf_totals['field'],se_n['field'])})"
    f"  greenhouse <b>{coinf_totals['greenhouse']:,}/{se_n['greenhouse']:,}</b> ({pct(coinf_totals['greenhouse'],se_n['greenhouse'])})"
    f"  controlled <b>{coinf_totals['controlled']:,}/{se_n['controlled']:,}</b> ({pct(coinf_totals['controlled'],se_n['controlled'])})"
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
        text    = "<b>Cryptic co-infection prevalence in field/greenhouse-collected plant RNA-seq</b><br>"
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

COLUMN_TITLES   = ["BioSamples", "Tissue", "Setting", "Stress", "Co-infection"]
COLUMN_X_SANKEY = [0.01,          0.18,     0.32,      0.55,     0.86]

for title, sx in zip(COLUMN_TITLES, COLUMN_X_SANKEY):
    fig.add_annotation(
        x=s2px(sx), y=0.03, xref="paper", yref="paper",
        text=f"<b>{title}</b>", showarrow=False,
        xanchor="center", yanchor="middle", font=dict(size=16),
    )

# ── SVG: add right-side labels for col4 nodes ─────────────────────────────────

def _svg_text(html, x, y, font_size=15, line_height=18):
    parts  = re.split(r'<br\s*/?>', html, flags=re.I)
    y0     = y - (len(parts) - 1) * line_height / 2
    tspans = []
    for i, part in enumerate(parts):
        bold   = bool(re.search(r'<b>', part, re.I))
        text   = re.sub(r'<[^>]+>', '', part).strip()
        weight = ' font-weight="bold"' if bold else ''
        tspans.append(f'<tspan x="{x:.1f}" y="{y0 + i*line_height:.1f}"{weight}>{text}</tspan>')
    return (
        f'<text text-anchor="start" dominant-baseline="central" '
        f'font-family="sans-serif" font-size="{font_size}" fill="#2a2a2a">'
        + ''.join(tspans) + '</text>'
    )


def fix_col4_labels(svg_path):
    plot_w  = W - M_L - M_R
    label_x = M_L + 0.86 * plot_w + 18 + 5
    inserts = [
        _svg_text(html, label_x, M_T + node_y[nidx] * (H - M_T - M_B))
        for nidx, html in col4_annot.items()
    ]
    content = svg_path.read_text()
    svg_path.write_text(content.replace('</svg>', '\n'.join(inserts) + '\n</svg>'))


# ── Export ─────────────────────────────────────────────────────────────────────

OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
fig.write_html(str(OUT_HTML))
print(f"\nWritten: {OUT_HTML}")

fig.write_image(str(OUT_SVG))
fix_col4_labels(OUT_SVG)
print(f"Written: {OUT_SVG}")

try:
    import cairosvg
    cairosvg.svg2png(url=str(OUT_SVG.resolve()), write_to=str(OUT_PNG), scale=2)
    print(f"Written: {OUT_PNG}")
except ImportError:
    print("PNG: pip install cairosvg  (SVG is ready)")
except Exception as e:
    print(f"PNG skipped: {e}")
