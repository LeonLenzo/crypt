#!/usr/bin/env python3
"""Every run is a column, ordered west to east: what was declared against what is there.

    python release/03_classify/figures/samples/plot_samples.py

Reads  ../../data/outputs/supp_data_s5_detections.tsv    the calls
       ../../../01_cohort/data/outputs/supp_data_s2_cohort.tsv   runs, location, declaration
       ../../../01_cohort/figures/country_centroids.tsv   to place a country-only run
Writes samples.png, samples_legend.png, samples_order.tsv  (beside this script)

THE UNIT IS THE RUN, NOT THE BIOPROJECT. A project is not a meaningful division of this
cohort: PRJEB39201 alone spans twelve countries, so grouping by project averages a worldwide
survey into one column and hides the only axis along which the complex actually varies.
Here every run is one pixel-wide column, ordered by COUNTRY (west to east by the mean
longitude of its runs) and by longitude within the country.

Colour is not abundance, it is WHO KNEW. A call is charcoal if the study that generated the
run declared that pathogen, orange if it did not, and grey if nobody has yet curated what the
study declared. The orange is the result and the charcoal is the control, and because the
column is a run rather than a project the two are read against each other sample by sample.

The grey is deliberately not hidden and not dropped. 18 BioProjects have no curated
declaration, and an orange cell there would assert something the curation does not support
while dropping them would quietly shrink the cohort. Grey says "detected, reporting unknown",
which is the honest third state.

Every number is read at draw time; none is transcribed.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _palette import WHEAT_PATH, INK, MUTED             # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from _cohort import is_analysed, describe               # noqa: E402
from detect import declaration                          # noqa: E402

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
CEN = HERE.parent / "01_cohort/figures/country_centroids.tsv"
GRP = HERE / "data/inputs/pathogen_groups.tsv"
FIG = Path(__file__).resolve().parent
csv.field_size_limit(10 ** 7)

NONE, DECL, UNDECL, UNKNOWN = 0, 1, 2, 3
COLOURS = ["#ffffff", "#3d3d3d", WHEAT_PATH, "#c2c2c2"]
CMAP = ListedColormap(COLOURS)

# The archive writes one country several ways. Normalising is not cosmetic: unmerged, 'UK' and
# 'United Kingdom' sort to opposite ends of the panel and split one country into two blocks.
# NAMES NORMALISE TO THE GAZETTEER'S SPELLING, not to the long form: country_centroids.tsv
# keys on 'UK' and 'USA', so mapping the other way produced a name it could not place.
# 'Zimbabwae' is a typo in the curated location field and is corrected here only for drawing;
# the fix belongs in 01_cohort's curation, and `main` prints it so it is not forgotten.
ALIAS = {"United Kingdom": "UK", "United States": "USA", "The Netherlands": "Netherlands",
         "Czechia": "Czech Republic", "Russian Federation": "Russia", "Zimbabwae": "Zimbabwe"}
NULL_LOC = {"", "missing", "unknown", "not applicable", "not collected", "na", "n/a"}
MIN_RUNS = 10           # a taxon needs this many calls across the cohort to get a row
LABEL_MIN = 18          # a country needs this many runs for its name to fit above the block
BAND_GAP = 0.26         # inches of white between bands; the band's name sits in it


def groups():
    """{taxon: (group_order, group)} from data/inputs/pathogen_groups.tsv.

    Rows band by DISEASE TYPE rather than by prevalence. Prevalence order interleaved the
    three rusts with the blotches, so a reader comparing like with like had to hop rows; the
    bands put the organisms that cause one kind of disease together and make the rust block
    legible as a block. The table also carries class and lifestyle, so the banding can be
    switched to either without re-curating.
    """
    out = {}
    with open(GRP) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("taxon\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                out[f[0]] = (int(f[2]), f[1])
    return out


def abbrev(name):
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


def country(loc):
    """Pull the country out of a free-text location field.

    Two shapes occur and both are in the cohort. The INSDC form puts the country first,
    'Serbia: Novi Sad'. A plain address puts it last, 'Rimski Sancevi, Novi Sad, Serbia',
    and reading the first token there yields a locality that places as its own country
    block at the far edge of the panel. Returns "" when the field carries no location.
    """
    loc = loc.strip()
    if loc.lower() in NULL_LOC:
        return ""
    head = loc.split(":")[0].strip() if ":" in loc else loc.split(",")[-1].strip()
    return "" if head.lower() in NULL_LOC else ALIAS.get(head, head)


def load():
    """-> (order, taxa, grid, meta): runs in drawing order and the state of every cell."""
    runs, excluded, decl_of = {}, 0, {}
    with open(S2) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not is_analysed(r):          # lab libraries and the host contrast; see _cohort
                excluded += 1
                continue
            try:
                lat, lon = float(r["lat"]), float(r["lon"])
            except ValueError:
                lat = lon = None
            runs[r["run"]] = {
                "project": r["BioProject"], "country": country(r["location"]),
                "lat": lat, "lon": lon, "status": r["declared_status"]}
            decl_of[r["run"]] = r          # the cohort row IS declaration()'s first argument

    cen = {r["country"]: float(r["lon"]) for r in
           csv.DictReader(open(CEN), delimiter="\t")}

    # Country order is the MEAN LONGITUDE OF ITS OWN RUNS where the runs carry coordinates,
    # and the gazetteer centroid only where none of them do. Using the centroid throughout
    # would place a country by its shape rather than by where the sampling actually happened.
    bycty = defaultdict(list)
    for rid, m in runs.items():
        bycty[m["country"]].append(rid)
    pos, unplaced, unknown_name = {}, [], []
    for c, ids in bycty.items():
        lons = [runs[i]["lon"] for i in ids if runs[i]["lon"] is not None]
        if lons:
            pos[c] = sum(lons) / len(lons)
            if c and c not in cen:
                unknown_name.append((c, len(ids)))   # placed, but the name is still unresolved
        elif c in cen:
            pos[c] = cen[c]
        else:
            unplaced.append(c)          # "" is no location at all and sorts first
            if c:
                unknown_name.append((c, len(ids)))
    ordered = sorted(pos, key=lambda c: pos[c]) + sorted(unplaced)

    # Within a country: by longitude, then latitude, then accession. The accession is the
    # tie-break and is what makes the column order reproducible for runs sharing a site.
    key = lambda i: (runs[i]["lon"] is None, runs[i]["lon"] or 0.0,
                     runs[i]["lat"] or 0.0, i)
    order = [i for c in ordered for i in sorted(bycty[c], key=key)]

    # `declared` is RECOMPUTED here with detect.py's own function rather than read off S5.
    # One definition, imported by the producer and by this figure, so the two cannot drift.
    # It is checked against the stored column and any disagreement is reported loudly: a
    # difference means S5 was written before detect.py's congener fix (2026-11-23) and the
    # published supplementary table is stale until detect.py is re-run on Setonix.
    hit, declared_in, mismatched = defaultdict(set), defaultdict(set), 0
    with open(S5) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["is_wheat_pathogen"] != "yes" or r["run"] not in runs:
                continue
            hit[r["taxon"]].add(r["run"])
            state = declaration(decl_of[r["run"]], r["taxon"], r["genus"], False)[0]
            mismatched += state != r["declared"]
            if state == "declared":
                declared_in[r["taxon"]].add(r["run"])

    grp = groups()
    missing = sorted(t for t in hit if len(hit[t]) >= MIN_RUNS and t not in grp)
    if missing:
        raise SystemExit(f"{GRP} has no group for: {', '.join(missing)}. Add them rather "
                         f"than letting the panel drop a pathogen it detected.")
    # band first, then prevalence within the band
    taxa = sorted((t for t in hit if len(hit[t]) >= MIN_RUNS),
                  key=lambda t: (grp[t][0], -len(hit[t]), t))
    dropped = {t: len(v) for t, v in sorted(hit.items(), key=lambda kv: (-len(kv[1]), kv[0]))
               if len(v) < MIN_RUNS}

    idx = {rid: j for j, rid in enumerate(order)}
    grid = np.zeros((len(taxa), len(order)), dtype=np.uint8)
    for i, t in enumerate(taxa):
        for rid in hit[t]:
            grid[i, idx[rid]] = (DECL if rid in declared_in[t]
                                 else UNKNOWN if runs[rid]["status"] == "not_curated"
                                 else UNDECL)
    return (order, taxa, grid, runs, ordered, bycty, dropped, sorted(unknown_name),
            excluded, [grp[t][1] for t in taxa], mismatched)


def main():
    (order, taxa, grid, runs, ordered, bycty, dropped, unknown_name, excluded,
     bands, mismatched) = load()
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.facecolor": "white",
                         "figure.facecolor": "white"})

    GW, RH = 7.10, 0.235                    # grid width and row height, in inches
    X0, Y0 = 1.42, 0.52                     # Y0 leaves room for the west-east axis beneath
    PADR = 0.95                             # the rightmost country name runs past the grid

    # Consecutive rows sharing a band, top to bottom. The gaps between bands are what makes
    # the banding visible, so they are part of the geometry rather than a drawn line.
    spans = []
    for i, b in enumerate(bands):
        if spans and spans[-1][2] == b:
            spans[-1][1] = i
        else:
            spans.append([i, i, b])
    gh = len(taxa) * RH + (len(spans) - 1) * BAND_GAP
    H = Y0 + gh + 1.30
    fig, ax = plt.subplots(figsize=(X0 + GW + PADR, H))
    top = Y0 + gh

    cw = GW / len(order)
    for k, (i0, i1, name) in enumerate(spans):
        b_top = top - (i0 * RH + k * BAND_GAP)
        b_bot = b_top - (i1 - i0 + 1) * RH
        # imshow rather than 20,000 patches: at ~1,560 columns a rectangle is a fraction of
        # a pixel wide and aliases away. Nearest-neighbour keeps every column present.
        ax.imshow(grid[i0:i1 + 1], cmap=CMAP, vmin=0, vmax=3, interpolation="nearest",
                  aspect="auto", extent=(X0, X0 + GW, b_bot, b_top), zorder=2)
        for i in range(i0, i1 + 1):
            ax.text(X0 - 0.09, b_top - (i - i0 + 0.5) * RH, abbrev(taxa[i]), ha="right",
                    va="center", fontsize=7.4, style="italic", color=INK)
        ax.text(X0 - 0.09, b_top + 0.045, name, ha="right", va="bottom", fontsize=6.8,
                color=MUTED)
        # country hairlines are drawn per band so they do not cross the gaps between bands
        at = 0
        for c in ordered:
            at += len(bycty[c])
            if at < len(order):
                ax.plot([X0 + at * cw] * 2, [b_bot, b_top], color="#ffffff",
                        linewidth=0.8, zorder=3)

    at = 0
    for c in ordered:
        n = len(bycty[c])
        if n >= LABEL_MIN:
            ax.text(X0 + (at + n / 2) * cw, top + 0.07, c or "no location", rotation=55,
                    rotation_mode="anchor", ha="left", va="baseline", fontsize=6.8,
                    color=INK)
        at += n
    ax.annotate("", xy=(X0 + GW, Y0 - 0.17), xytext=(X0, Y0 - 0.17),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, linewidth=0.8,
                                shrinkA=0, shrinkB=0))
    ax.text(X0 + GW / 2, Y0 - 0.24, f"{len(order):,} {describe()} runs, west to east",
            ha="center", va="top", fontsize=7.0, color=MUTED)

    ax.set_xlim(0, X0 + GW + PADR)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "samples.png", dpi=600)

    legend()

    with open(FIG / "samples_order.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["column", "run", "BioProject", "country", "lat", "lon",
                    "declared_status", "n_calls"])
        for j, rid in enumerate(order):
            m = runs[rid]
            w.writerow([j, rid, m["project"], m["country"], m["lat"], m["lon"],
                        m["status"], int((grid[:, j] > 0).sum())])

    state = {k: int((grid == v).sum()) for k, v in
             (("declared", DECL), ("undeclared", UNDECL), ("unknown", UNKNOWN))}
    print(f"wrote {FIG / 'samples.png'}")
    print(f"  {len(order):,} {describe()} runs ({excluded} excluded), "
          f"{len(ordered)} countries, {len(taxa)} taxa drawn")
    print(f"  cells: {state}")
    print(f"  below MIN_RUNS={MIN_RUNS}, not drawn: " +
          ", ".join(f"{abbrev(t)} ({n})" for t, n in dropped.items()))
    noloc = len(bycty.get("", ()))
    print(f"  {noloc} runs carry no location and sit in the unplaced block")
    if mismatched:
        print(f"  !! {mismatched} calls disagree with S5's stored `declared` column. This "
              f"figure recomputes it with detect.py's post-2026-11-23 congener fix; "
              f"supp_data_s5/s6 are STALE until detect.py is re-run on Setonix.")
    if unknown_name:
        print("  COUNTRY NAMES THE GAZETTEER CANNOT RESOLVE -- fix in 01_cohort curation: " +
              ", ".join(f"{c!r} ({n})" for c, n in unknown_name))


def legend():
    """Keys laid out by MEASURING each label, not by estimating from its character count.

    The first version advanced x by a per-character guess and the last key ran off the
    canvas. Matplotlib will report a text's drawn width, so ask it: draw, measure, place the
    next swatch past the measurement.
    """
    SW, GAP, PAD = 0.15, 0.13, 0.30         # swatch, swatch-to-text, key-to-key, in inches
    H = 0.30
    fig, ax = plt.subplots(figsize=(8.0, H))     # provisional; trimmed to content below
    ax.set_xlim(0, 8.0)
    ax.set_ylim(0, H)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    rend = fig.canvas.get_renderer()

    x = 0.08
    for colour, label in ((COLOURS[DECL], "declared by the study"),
                          (COLOURS[UNDECL], "not declared"),
                          (COLOURS[UNKNOWN], "declaration not curated")):
        ax.add_patch(Rectangle((x, 0.09), SW, 0.13, facecolor=colour,
                               edgecolor="#d8d8d8", linewidth=0.5))
        t = ax.text(x + SW + GAP, 0.155, label, ha="left", va="center", fontsize=7.4,
                    color=INK)
        x += SW + GAP + t.get_window_extent(rend).width / fig.dpi + PAD

    W = x - PAD + 0.08
    fig.set_size_inches(W, H)
    ax.set_xlim(0, W)
    fig.savefig(FIG / "samples_legend.png", dpi=300, facecolor="white")


if __name__ == "__main__":
    main()
