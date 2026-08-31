#!/usr/bin/env python3
"""
Literature resolution pipeline Sankey — meta_search.py -> meta_text.py.

True funnel, matching meta_search.py's actual control flow (resolve_bioproject()):
a BioProject only reaches "Web search (Serper)" if the NCBI stage found nothing
at all (Serper only ever fires when pmid_source == "none") — "DOI only (CrossRef)"
is purely an NCBI-XML outcome (a bare <Publication> DOI with no PMID) and never
touches Serper. So the honest sequence is:

  Origin (1,286)
    -> "NCBI (XML/PMC/PubMed)"   peels off (has PMID)
    -> "DOI only (CrossRef)"      peels off (has DOI, no PMID — still NCBI-XML)
    -> "Web search (Serper)"      gets the NCBI-stage remainder; peels off its
                                   own successes; whatever's left is genuinely
                                   exhausted -> "No DOI found"

Every BP that ends up with a DOI (NCBI + DOI-only + Serper) then enters
meta_text.py, which peels off PMC OA / Unpaywall PDF / Manual PDF successes;
the rest is "No OA copy found".

Exactly two true terminals, per design: "Full text retrieved" (fed by PMC OA
+ Unpaywall + Manual) and "No full text available" (fed by "No DOI found" +
"No OA copy found") — provenance stays visible as the flows leading into them,
it just no longer pretends to be a terminal outcome itself.

This figure is upstream of sample_funnel — that one covers what happens
*after* text is in hand (classify_metadata.py's LLM categorisation).

Run from crypt/: python metadata/figures/lit_resolution_sankey.py
"""

import json
import sys
from collections import Counter
from pathlib import Path

try:
    import plotly.graph_objects as go
except ImportError:
    sys.exit("pip install plotly kaleido")

BIOPROJECTS_PATH = Path("metadata/output/meta_search/data/bioprojects.json")
BIOSAMPLES_PATH  = Path("metadata/output/meta_search/data/biosamples.json")
TEXT_CACHE_PATH  = Path("metadata/output/meta_text/data/text_cache.jsonl")
OUT_DIR  = Path("metadata/output/figures/sankey")
OUT_HTML = OUT_DIR / "lit_resolution_sankey.html"
OUT_PNG  = OUT_DIR / "lit_resolution_sankey.png"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Load ──────────────────────────────────────────────────────────────────────

bps     = json.loads(BIOPROJECTS_PATH.read_text())
bs_data = json.loads(BIOSAMPLES_PATH.read_text())

text_cache: dict[str, dict] = {}
with open(TEXT_CACHE_PATH, encoding="utf-8") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        doi, rest = line.split("\t", 1)
        text_cache[doi] = json.loads(rest)   # last line wins (append-only cache)

_ENA = {"SAME", "SAMD"}
total_bps    = len(bps)
ena_bs_count = sum(1 for a in bs_data if a[:4] in _ENA or a[:5] in _ENA)
ena_bs_desc  = sum(1 for a, v in bs_data.items()
                   if (a[:4] in _ENA or a[:5] in _ENA) and v.get("bs_description"))

# ── Classify each BP's provenance, strictly on has_doi + pmid_source ──────────
# (every non-"none" pmid_source has a doi — verified, 0 exceptions — so the
# only ambiguity lives inside "none", which we resolve by doi presence: 14
# migration-era cache entries carry a doi despite pmid_source=="none" and are
# folded into "DOI only" rather than a bucket that's supposed to mean "no doi")

NCBI_SRCS = {"bioproject_xml", "pmc_search", "pubmed_search"}


def _provenance(v: dict) -> str:
    src = v.get("pmid_source") or "none"
    if src in NCBI_SRCS:
        return "NCBI (XML/PMC/PubMed)"
    if src.startswith("serper"):
        return "Web search (Serper)"
    return "DOI only (CrossRef)"   # doi_only→*, bioproject_xml_doi_only, or "none" w/ a doi


def _text_outcome(doi: str, has_full_text: bool) -> str:
    if not has_full_text:
        return "no_oa"
    pdf_url = text_cache.get(doi, {}).get("pdf_url", "")
    if pdf_url.startswith("epmc:"):
        return "pmc_oa"
    if pdf_url.startswith("manual:"):
        return "manual"
    return "unpaywall"


PROV_ORDER = ["NCBI (XML/PMC/PubMed)", "DOI only (CrossRef)", "Web search (Serper)"]
# Fixed internal keys — used for lookups everywhere below. DO NOT rename these
# to change what's shown on the figure; edit TEXT_LABELS instead (it's the
# only place display text lives for this stage).
TEXT_ORDER = ["pmc_oa", "unpaywall", "manual", "no_oa"]
TEXT_LABELS = {
    "pmc_oa":    "PMC",
    "unpaywall": "Unpaywall",
    "manual":    "Manual",
    "no_oa":     "Unavailable",
}

prov_counts: Counter = Counter()
prov_to_text: dict[str, Counter] = {p: Counter() for p in PROV_ORDER}
no_doi_count = 0

for bp, v in bps.items():
    doi = v.get("doi")
    if not doi:
        no_doi_count += 1
        continue
    prov = _provenance(v)
    prov_counts[prov] += 1
    outcome = _text_outcome(doi, bool(v.get("full_text")))
    prov_to_text[prov][outcome] += 1

has_doi_total = sum(prov_counts.values())
assert has_doi_total + no_doi_count == total_bps

text_totals = Counter()
for p in PROV_ORDER:
    for t, n in prov_to_text[p].items():
        text_totals[t] += n

full_text_total = text_totals["pmc_oa"] + text_totals["unpaywall"] + text_totals["manual"]
no_full_text_total = text_totals["no_oa"] + no_doi_count
assert full_text_total + no_full_text_total == total_bps

# NCBI stage remainder = what actually goes on to attempt Web search
ncbi_stage_total   = prov_counts["NCBI (XML/PMC/PubMed)"] + prov_counts["DOI only (CrossRef)"]
web_search_attempt = total_bps - ncbi_stage_total   # = prov_counts["Web search (Serper)"] + no_doi_count

# ── Colours ───────────────────────────────────────────────────────────────────
# Two-tone scheme: script-boundary nodes (meta_search.py / meta_text.py) stay
# neutral dark slate; everything still in the running ("included") is the
# same lighter blue as Database search; everything that has fallen out
# ("excluded") is grey. Link colours just follow their source node.

SCRIPT_COLOR   = "#2c3e50"
INCLUDED_COLOR = "#2980b9"
EXCLUDED_COLOR = "#95a5a6"
SPINE_Y = 0.06


def rgba(hex_color: str, alpha: float = 0.38) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


# ── Build nodes ───────────────────────────────────────────────────────────────
#
# arrangement="fixed" (set below, in the Figure section) means Plotly uses the
# x/y given here EXACTLY — no auto-avoidance of overlaps like "snap" does.
# That's what makes nodes freely repositionable, but it also means the *default*
# y values below have to be real, value-proportional midpoints or nodes will
# clip past each other / off the canvas. place_nodes() computes those defaults;
# after that, every y is just a plain number in `ys` — edit any of them by hand
# for manual placement, nothing downstream depends on the specific value.

labels, colors, xs, ys = [], [], [], []


def add_node(label, color, x, y=0.5):
    i = len(labels)
    labels.append(label); colors.append(color); xs.append(x); ys.append(y)
    return i


def place_nodes(pairs, margin=0.04, gap=0.03):
    """Assign value-proportional y midpoints (with explicit gaps reserved
    between nodes) to a list of (node_idx, value) pairs sharing a column."""
    total_val = sum(v for _, v in pairs) or 1
    total_gap = gap * (len(pairs) - 1)
    avail = max(0.0, (1 - 2 * margin) - total_gap)
    y = margin
    for idx, v in pairs:
        h = (v / total_val) * avail
        ys[idx] = y + h / 2
        y += h + gap


COL_X = {"gateway": 0.03, "dbsearch": 0.22, "websearch": 0.40,
         "spine2": 0.60, "stage2": 0.80, "terminal": 0.96}

# meta_search.py — the literal first node now; no separate raw-count node
gateway_idx = add_node(f"<b>meta_search.py</b><br>{total_bps:,} BPs", SCRIPT_COLOR, COL_X["gateway"])

# Database search (NCBI) — the actual first resolution attempt. Resolved BPs
# (PMID found, or a bare DOI in the XML enriched via CrossRef — both are
# NCBI-XML-derived, so they're not split out separately here) flow straight
# into meta_text.py; everything else falls through to Web search. ENA is NOT
# part of this: meta_search.py's ENA/EBI calls only fetch BioSample attributes
# (tissue, geo_loc, collection_date, host, bs_description) — a separate
# concern from DOI/PMID resolution, shown only in the caption note below.
dbsearch_idx = add_node(f"<b>Database</b><br>{total_bps:,} BPs", INCLUDED_COLOR, COL_X["dbsearch"])

# Web search (Serper) — only ever attempted when NCBI found nothing at all
web_idx = add_node(
    f"<b>Web</b><br>{web_search_attempt:,} BPs attempted"
    f"<br><sup>(NCBI found nothing)</sup>",
    INCLUDED_COLOR, COL_X["websearch"],
)

# meta_text.py spine — every BP with a DOI, regardless of provenance
spine_idx = add_node(f"<b>meta_text.py</b><br>{has_doi_total:,} BPs with a DOI", SCRIPT_COLOR, COL_X["spine2"])

# No DOI found — genuinely exhausted after both search stages
no_doi_idx = add_node(f"<b>No DOI</b><br>{no_doi_count:,} BPs", EXCLUDED_COLOR, COL_X["spine2"])

# Stage 2 — meta_text outcomes (still intermediate, not terminals).
# Keys here are the fixed TEXT_ORDER keys, not display text — see TEXT_LABELS
# above to rename what's shown on the figure.
TEXT_NODE_COLOR = {**{t: INCLUDED_COLOR for t in ("pmc_oa", "unpaywall", "manual")},
                    "no_oa": EXCLUDED_COLOR}
text_idx = {}
for t in TEXT_ORDER:
    n = text_totals[t]
    text_idx[t] = add_node(f"<b>{TEXT_LABELS[t]}</b><br>{n:,} BPs", TEXT_NODE_COLOR[t], COL_X["stage2"])

# The two true terminals
full_text_idx = add_node(
    f"<b>Full text</b><br>{full_text_total:,} BPs ({100*full_text_total/total_bps:.1f}%)",
    INCLUDED_COLOR, COL_X["terminal"],
)
no_full_text_idx = add_node(
    f"<b>No  text</b><br>{no_full_text_total:,} BPs ({100*no_full_text_total/total_bps:.1f}%)",
    EXCLUDED_COLOR, COL_X["terminal"],
)

# ── Position nodes: value-proportional defaults, edit `ys[...]` by hand below
# for manual placement (e.g. ys[web_idx] = 0.25) ───────────────────────────────

place_nodes([(gateway_idx, total_bps)])
place_nodes([(dbsearch_idx, total_bps)])
place_nodes([(web_idx, web_search_attempt)])
place_nodes([(spine_idx, has_doi_total), (no_doi_idx, no_doi_count)])
place_nodes([(text_idx[t], text_totals[t]) for t in TEXT_ORDER])
place_nodes([(full_text_idx, full_text_total), (no_full_text_idx, no_full_text_total)])

# NOTE on stacking order at a shared target (e.g. "No  text" receives both
# "No DOI" and "Unavailable"): Plotly stacks by the SOURCE node's y, smaller
# y first. There's no independent "stacking order" control — moving a
# source's y to change its stack position also moves it within its own
# column, which can look worse than the stacking order it fixes (tried
# exactly that for "Unavailable" here; reverted — it broke the clean
# size-ordering of the meta_text column for a cosmetic gain at the far
# right). Leaving this one as Plotly's default proportional placement gives.

# ── Build links ───────────────────────────────────────────────────────────────

srcs, tgts, vals, lc = [], [], [], []


def add_link(s, t, v, c, alpha=0.38):
    if v <= 0:
        return
    srcs.append(s); tgts.append(t); vals.append(v); lc.append(rgba(c, alpha))


# meta_search.py -> Database search (NCBI)
add_link(gateway_idx, dbsearch_idx, total_bps, SCRIPT_COLOR)

# Database search (NCBI): resolved BPs go straight into meta_text.py, the
# NCBI-exhausted remainder visibly flows into Web search
ncbi_resolved = prov_counts["NCBI (XML/PMC/PubMed)"] + prov_counts["DOI only (CrossRef)"]
add_link(dbsearch_idx, spine_idx, ncbi_resolved, INCLUDED_COLOR, alpha=0.45)
add_link(dbsearch_idx, web_idx, web_search_attempt, INCLUDED_COLOR)

# Web search (Serper): resolved -> meta_text.py, exhausted -> No DOI found
add_link(web_idx, spine_idx, prov_counts["Web search (Serper)"], INCLUDED_COLOR, alpha=0.45)
add_link(web_idx, no_doi_idx, no_doi_count, EXCLUDED_COLOR, alpha=0.42)

# meta_text spine -> stage 2 outcomes
for t in TEXT_ORDER:
    add_link(spine_idx, text_idx[t], text_totals[t], SCRIPT_COLOR, alpha=0.42)

# Stage 2 -> the two true terminals
for t in ("pmc_oa", "unpaywall", "manual"):
    add_link(text_idx[t], full_text_idx, text_totals[t], INCLUDED_COLOR, alpha=0.5)
add_link(text_idx["no_oa"], no_full_text_idx, text_totals["no_oa"], EXCLUDED_COLOR, alpha=0.5)
add_link(no_doi_idx, no_full_text_idx, no_doi_count, EXCLUDED_COLOR, alpha=0.5)

# ── Figure ────────────────────────────────────────────────────────────────────

ena_note = (
    f"ENA BioSample enrichment (EBI API): "
    f"{ena_bs_count:,} SAME*/SAMD* BioSamples · "
    f"{ena_bs_desc:,} with description ({100*ena_bs_desc/ena_bs_count:.0f}%)"
)

fig = go.Figure(go.Sankey(
    # "fixed" = use x/y exactly as given, no auto-repositioning to dodge
    # overlaps. place_nodes() above already reserves proper value-proportional
    # gaps, so this is what makes manual y edits actually stick.
    arrangement="fixed",
    node=dict(
        pad=6,
        thickness=18,
        line=dict(color="white", width=0.5),
        label=labels, color=colors, x=xs, y=ys,
    ),
    link=dict(source=srcs, target=tgts, value=vals, color=lc),
))

fig.update_layout(
    title=dict(
        text=(
            f"Literature resolution — meta_search.py → meta_text.py — {total_bps:,} BioProjects<br>"
            f"<sup>"
            f"Full text retrieved: {full_text_total:,} ({100*full_text_total/total_bps:.1f}%)"
            f" &nbsp;|&nbsp; "
            f"No full text available: {no_full_text_total:,} ({100*no_full_text_total/total_bps:.1f}%)<br>"
            f"{ena_note}"
            f"</sup>"
        ),
        font=dict(size=15), x=0.01,
    ),
    font=dict(size=12, family="Arial"),
    paper_bgcolor="white",
    width=1400, height=700,
    margin=dict(l=10, r=10, t=100, b=30),
)

fig.write_html(OUT_HTML)
print(f"Written {OUT_HTML}")
try:
    fig.write_image(OUT_PNG, scale=2)
    print(f"Written {OUT_PNG}")
except Exception as e:
    print(f"PNG skipped ({e})")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(f"Stage 1 — meta_search  ({total_bps:,} BioProjects)")
print(f"  NCBI (XML/PMC/PubMed)     {prov_counts['NCBI (XML/PMC/PubMed)']:>5,}  ({100*prov_counts['NCBI (XML/PMC/PubMed)']/total_bps:.1f}%)")
print(f"  DOI only (CrossRef)       {prov_counts['DOI only (CrossRef)']:>5,}  ({100*prov_counts['DOI only (CrossRef)']/total_bps:.1f}%)")
print(f"  Web search (Serper) — of {web_search_attempt:,} attempted (NCBI found nothing):")
print(f"    resolved                {prov_counts['Web search (Serper)']:>5,}  ({100*prov_counts['Web search (Serper)']/max(web_search_attempt,1):.1f}% of attempted)")
print(f"    No DOI found            {no_doi_count:>5,}  ({100*no_doi_count/max(web_search_attempt,1):.1f}% of attempted)  <- funnel exit")
print()
print(f"Stage 2 — meta_text  (of {has_doi_total:,} with a DOI)")
for t in TEXT_ORDER:
    n = text_totals[t]
    print(f"  {TEXT_LABELS[t]:<20} {n:>5,}  ({100*n/max(has_doi_total,1):.1f}%)")
print()
print(f"Terminals:")
print(f"  Full text retrieved:      {full_text_total:,}  ({100*full_text_total/total_bps:.1f}%)")
print(f"  No full text available:   {no_full_text_total:,}  ({100*no_full_text_total/total_bps:.1f}%)")
print()
print(f"ENA enrichment: {ena_bs_count:,} SAME*/SAMD* BioSamples · {ena_bs_desc:,} with description")
