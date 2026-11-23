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

import collections
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _palette import MATRIX, INK, MUTED, darken         # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from _cohort import is_analysed, describe               # noqa: E402
from detect import declaration                          # noqa: E402

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
CEN = HERE.parent / "01_cohort/figures/country_centroids.tsv"
GRP = HERE / "data/inputs/pathogen_groups.tsv"
REGIONS = HERE / "data/inputs/country_regions.tsv"
S3REF = HERE.parent / "02_reference/data/outputs/supp_data_s3_reference.tsv"
FIG = Path(__file__).resolve().parent
csv.field_size_limit(10 ** 7)

NONE, DECL, UNDECL, UNKNOWN = 0, 1, 2, 3
DECLARED_COLOUR = "#3d3d3d"     # the study reported it; the one state that is not a category

# A CELL TAKES THE COLOUR THE REFERENCE MATRIX GAVE THAT TAXON. `MATRIX` is keyed on the two
# facts S3 records, where USDA recorded the taxon and what PHI-base scores it as, and is
# imported from 02_reference rather than copied so the two figures cannot drift apart.
#
# The matrix was drawn as inch-wide pills with dark borders; here a cell is one pixel wide,
# and two of its colours (#f0f0f0 for a congener, #c3e8bf for a plant resident) vanish
# against white at that size. LUM_CEILING darkens any colour past it until it reaches the
# ceiling, which preserves the hue and the ordering and only changes what cannot be seen.
# The BAND SWATCH keeps the true matrix colour, so the correspondence is still legible.
# S3 spells the PHI-base axis in full; MATRIX is keyed on plot_matrix.py's short names. The
# two have to be translated, not assumed equal -- the panel raised a KeyError on the first
# wheat pathogen it met, which is the good failure mode for a mismatch like this.
PHI_KEY = {"none": "none", "plant_pathogen": "vagrant", "wheat_pathogen": "wheat_path"}

# TWO BANDING RULES, because the two halves of the panel answer to different things.
# A WHEAT PATHOGEN is banded by the DISEASE it causes, which is the distinction a plant
# pathologist makes and which no column of S3 carries, so it is authored in
# pathogen_groups.tsv. Everything else is banded by its MATRIX SLOT, read straight off S3:
# the band then says exactly what the colour says, nobody has to invent a category, and a
# taxon cannot be placed by opinion. Microdochium majus is the case that forced this -- it
# is a real Microdochium patch agent, but PHI-base does not score it, so under the slot rule
# it sits with the other taxa USDA records on wheat rather than where a reader of the
# literature would file it. That disagreement is the label's, and the figure should show it
# rather than paper over it.
SLOT_BAND = {
    ("wheat", "none"):       (5, "Wheat resident"),
    ("wheat", "vagrant"):    (6, "Pathogen of another plant, recorded on wheat"),
    ("plant", "none"):       (7, "Plant resident"),
    ("plant", "vagrant"):    (8, "Pathogen of another plant"),
}

# Taxa in the ('congener', 'none') slot are NOT DRAWN. They are in the reference only because
# their GENUS qualified, so USDA holds no host record for the taxon itself and there is no
# statement the panel can make about where it belongs. Leon's call, 2026-11-23.
#
# This is a real cost, not only tidying: it removes Dioszegia hungarica, called in 353 runs,
# which 02_reference's scope notes name as part of the wheat resident core. The calls are not
# being judged false, they are being left out of a figure whose every other row carries a
# host record. supp_data_s5 keeps all of them.
DROP_SLOTS = {("congener", "none")}

# The two halves of the panel, bracketed on the far left.
SUPER = {True: "Pathogens", False: "Residents"}
LUM_CEILING = 0.70
SPILL_RATIO = 0.25      # never-alone AND below this share of the congener's coverage
SPILL = "Zymoseptoria brevis"   # what the rule catches here; the caption reads it from this


def _lum(hexc):
    r, g, b = (int(hexc.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def legible(hexc):
    """Darken toward black until the colour is visible as a one-pixel column."""
    out = hexc
    while _lum(out) > LUM_CEILING:
        out = darken(out, 0.94)
    return out


def washed(hexc, f=0.55):
    """The same hue, blended toward white: 'this category, reporting status unknown'."""
    h = hexc.lstrip("#")
    return "#%02x%02x%02x" % tuple(
        round(int(h[i:i + 2], 16) + (255 - int(h[i:i + 2], 16)) * f) for i in (0, 2, 4))

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
# COLUMNS GROUP BY REGION, NOT BY COUNTRY. 41 countries across a 7-inch grid left most
# blocks too narrow to name -- 18 of them, including Ireland and the Netherlands, which the
# caption names -- and the signal is not national anyway: Septoria tritici blotch is a
# maritime northwest European disease and stem rust an East African and Southern Cone one.
# The mapping is authored in data/inputs/country_regions.tsv; a country in the data and
# missing from it is a hard error. The country of every run stays in samples_order.tsv.
#
# A region still needs a minimum width to carry a name, and WHAT FITS IS GEOMETRY. Labels sit
# at LABEL_ROT, so two on one baseline clear each other only if their blocks are at least
# (text height / sin(rot)) apart, about 0.115 in here. Labels are placed on two staggered
# baselines and the separation is MEASURED: one that still would not clear is skipped and
# reported, so a collision is never silent.
LABEL_MIN = 12
LABEL_ROT = 55
LABEL_TIER = 0.30       # vertical offset of the second baseline, in inches
REGION_GAP = 0.085      # gap between two region blocks; the rule sits in the middle of it
BAND_GAP = 0.26         # inches of white between bands; the band's name sits in it

# The host is excluded: every library has one, no study declares it as a pathogen, and a row
# that is charcoal-free and full in every column says nothing. `plant_pathogen` IS drawn even
# though it is nearly empty here -- 5 taxa and 12 calls across the analysed set, none reaching
# MIN_RUNS -- because its emptiness is a result of db2's artefact removal and is reported in
# the caption rather than hidden by never looking.
DRAWN_CLASSES = {"wheat_pathogen", "plant_pathogen", "resident"}



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
            if len(f) >= 7:
                out[f[0]] = (int(f[2]), f[1], f[6] == "yes")
    return out


def abbrev(name):
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


def regions():
    """{country: region} from data/inputs/country_regions.tsv."""
    out = {}
    with open(REGIONS) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("country\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
    return out


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
    s3 = {r["taxon"]: (r["usda_row"], r["phibase_col"]) for r in
          csv.DictReader(open(S3REF), delimiter="\t")}
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
    reg = regions()
    unmapped = sorted({m["country"] for m in runs.values()
                       if m["country"] and m["country"] not in reg})
    if unmapped:
        raise SystemExit(f"{REGIONS} has no region for: {', '.join(unmapped)}. Add them; a "
                         f"silent 'Other' bucket would hide where the runs came from.")
    for m in runs.values():
        m["region"] = reg.get(m["country"], "")

    bycty = defaultdict(list)
    for rid, m in runs.items():
        bycty[m["region"]].append(rid)
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
    grp = groups()
    hit, declared_in, mismatched = defaultdict(set), defaultdict(set), 0
    cluster, cvs, collapsed_of = defaultdict(set), defaultdict(list), defaultdict(set)
    klass, cov, gen_of = {}, defaultdict(dict), {}
    slot_votes = defaultdict(collections.Counter)
    with open(S5) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["run"] not in runs or r["taxon_class"] not in DRAWN_CLASSES:
                continue
            klass[r["taxon"]] = r["taxon_class"]
            cluster[(r["run"], r["genus"])].add(r["taxon"])
            gen_of[r["genus"] if grp.get(r["taxon"], (0, "", False))[2]
                   else r["taxon"]] = r["genus"]
            if r["genus_cov_cv"] not in ("", "nan"):
                cvs[r["genus"]].append(float(r["genus_cov_cv"]))
            # A collapsed genus is ONE row: its species-level calls are one organism spread
            # over near-identical references, so drawing fourteen rows would show fourteen
            # co-infecting congeners that are not there.
            key = r["genus"] if grp.get(r["taxon"], (0, "", False))[2] else r["taxon"]
            collapsed_of[key].add(r["taxon"])
            hit[key].add(r["run"])
            cov[key][r["run"]] = float(r["coverage"])
            slot_votes[key][s3.get(r["taxon"], ("congener", "none"))] += 1
            state = declaration(decl_of[r["run"]], r["taxon"], r["genus"],
                                r["taxon_class"] == "host")[0]
            mismatched += state != r["declared"]
            if state == "declared":
                declared_in[key].add(r["run"])

    # The collapse decision is RE-MEASURED here rather than trusted: a genus whose species
    # are called together in one run at near-identical coverage is one organism, and if an
    # uncollapsed genus starts behaving that way the table is out of date, not the data.
    suspect = []
    for g, cv in cvs.items():
        sizes = [len(v) for (run, gg), v in cluster.items() if gg == g]
        med_cv, med_n = statistics.median(cv), statistics.median(sizes)
        collapsed = any(grp.get(t, (0, "", False))[2] for t in collapsed_of.get(g, ()))
        if med_n >= 3 and med_cv < 0.35 and not collapsed and len(hit.get(g, ())) >= MIN_RUNS:
            suspect.append((g, med_n, med_cv))

    # A collapsed genus has no row in S3, and taking the slot of whichever member happened to
    # be read last made its band depend on file order. Take the MODAL slot of the members,
    # the slot name breaking a tie, so the result is reproducible.
    mrow = {t: (lambda k: (k[0], PHI_KEY[k[1]]))(
                min(v.most_common(), key=lambda kv: (-kv[1], kv[0]))[0])
            for t, v in slot_votes.items()}

    # CONGENER-TRACKING, and it is a FLAG, NEVER A FILTER. The module's standing rule is
    # that a congener cluster is reported and not removed, and an earlier version of this
    # script broke it by deleting the row.
    #
    # The flag marks a taxon never called without a congener and only ever at a small
    # fraction of that congener's coverage. In this cohort it catches one, Zymoseptoria
    # brevis. What it is NOT is k-mer sharing between the two references: Kraken2's LCA
    # already sends every minimizer the two share to the genus node, and in a passing run
    # that node carries about 1.29 M of them while brevis still holds 100 k that are
    # brevis-EXCLUSIVE in db2. Those cannot have come from the Z. tritici reference.
    #
    # The mechanism is reference-strain divergence: field Z. tritici differs from the
    # reference genome, and where a field variant coincides with the brevis reference the
    # k-mer is brevis-exclusive and is called brevis. No coverage floor separates that, the
    # reads being real and the k-mers genuinely exclusive, which is why this is disclosed
    # rather than thresholded. It is the close-relative false positive TIPP-SD (Shen 2026)
    # measures and leaves uncorrected, and release/04_align is what resolves it, from
    # species-unique reads and positional evenness.
    spill = {}
    for t, rr in hit.items():
        if len(rr) < MIN_RUNS:
            continue
        sibs = [x for x in hit if x != t and gen_of.get(x) == gen_of.get(t)]
        if not sibs or (rr - set().union(*(hit[x] for x in sibs))):
            continue                                    # called alone at least once
        best = max(sibs, key=lambda x: len(rr & hit[x]))
        pairs = [cov[t][r] / cov[best][r] for r in rr & hit[best] if cov[best].get(r)]
        if pairs and statistics.median(pairs) < SPILL_RATIO:
            spill[t] = (len(rr), best, statistics.median(pairs))

    is_wp = lambda t: mrow[t][1] == "wheat_path"
    missing = sorted(t for t in hit
                     if len(hit[t]) >= MIN_RUNS and is_wp(t) and t not in grp)
    if missing:
        raise SystemExit(f"{GRP} has no disease band for the wheat pathogen(s): "
                         f"{', '.join(missing)}. Add them rather than letting the panel "
                         f"drop a pathogen it detected.")
    # band first, then prevalence within the band
    # A collapsed genus must have its OWN row in the table. Inheriting the band from
    # whichever species sorted first was order-dependent and silently wrong.
    for t in hit:
        if (len(hit[t]) >= MIN_RUNS and len(collapsed_of[t]) > 1 and is_wp(t)
                and t not in grp):
            raise SystemExit(f"{GRP} collapses {t} to a genus but has no row for the genus "
                             f"itself. Add one; the disease band must be chosen, not "
                             f"inherited.")
    band = lambda t: grp[t][:2] if is_wp(t) else SLOT_BAND[mrow[t]]
    dropped_slot = {t: len(hit[t]) for t in list(hit)
                    if len(hit[t]) >= MIN_RUNS and not is_wp(t) and mrow[t] in DROP_SLOTS}
    for t in dropped_slot:
        del hit[t]
    taxa = sorted((t for t in hit if len(hit[t]) >= MIN_RUNS),
                  key=lambda t: (band(t)[0], -len(hit[t]), t))
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
            excluded, [band(t)[1] for t in taxa],
            [is_wp(t) for t in taxa], mismatched,
            {t: len(collapsed_of[t]) for t in taxa if len(collapsed_of[t]) > 1}, suspect,
            # a band's category is the modal class of the taxa in it, read off the data
            # rather than authored a second time in the grouping table
            {t: collections.Counter(
                klass[x] for x in (collapsed_of[t] or {t})).most_common(1)[0][0]
             for t in taxa}, {t: mrow[t] for t in taxa}, spill, dropped_slot)


def main():
    (order, taxa, grid, runs, ordered, bycty, dropped, unknown_name, excluded,
     bands, supers, mismatched, collapsed, suspect, tclass, mrow, spill,
     dropped_slot) = load()
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.facecolor": "white",
                         "figure.facecolor": "white"})

    GW, RH = 7.10, 0.235                    # grid width and row height, in inches
    X0, Y0 = 1.50, 0.52                     # Y0 leaves room for the west-east axis beneath
    SUPX = 0.10                             # the pathogen/resident bracket, from the left
    BANDX = 0.34                            # band names, indented inside the bracket
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
    H = Y0 + gh + 1.30 + LABEL_TIER
    fig, ax = plt.subplots(figsize=(X0 + GW + PADR, H))
    top = Y0 + gh

    # REGIONS ARE DRAWN AS SEPARATE BLOCKS WITH REAL GAPS BETWEEN THEM. A white hairline
    # over one continuous strip was doing this job before, and at one-pixel columns a line
    # is indistinguishable from an empty run: the separator read as missing data. A gap of
    # actual whitespace cannot be misread, and it also stops the eye carrying a run of
    # detections across a boundary it does not cross.
    #
    # The gaps come OUT of the grid width, so the columns narrow rather than the panel
    # widening, and every region keeps its true share of the remaining space.
    nblk = len(ordered)
    cw = (GW - (nblk - 1) * REGION_GAP) / len(order)
    xs, at = [], 0.0
    for c in ordered:
        xs.append(X0 + at)
        at += len(bycty[c]) * cw + REGION_GAP
    starts = {}
    pos = 0
    for c in ordered:
        starts[c] = pos
        pos += len(bycty[c])

    for k, (i0, i1, name) in enumerate(spans):
        b_top = top - (i0 * RH + k * BAND_GAP)
        b_bot = b_top - (i1 - i0 + 1) * RH
        # imshow rather than 20,000 patches: at ~1,560 columns a rectangle is a fraction of
        # a pixel wide and aliases away. Nearest-neighbour keeps every column present.
        # Each ROW has its own colour now, so the block is built as RGB rather than drawn
        # through one colormap: hue says which matrix category the taxon is, and the three
        # states vary only the treatment of that hue.
        blk = np.ones((i1 - i0 + 1, grid.shape[1], 3))
        for i in range(i0, i1 + 1):
            base = legible(MATRIX[mrow[taxa[i]]])
            for state, col in ((DECL, DECLARED_COLOUR), (UNDECL, base),
                               (UNKNOWN, washed(base))):
                blk[i - i0][grid[i] == state] = to_rgb(col)
        for bi, c in enumerate(ordered):
            c0 = starts[c]
            c1 = c0 + len(bycty[c])
            x1 = xs[bi] + (c1 - c0) * cw
            ax.imshow(blk[:, c0:c1], interpolation="nearest", aspect="auto",
                      extent=(xs[bi], x1, b_bot, b_top), zorder=2)
            # A RULE CENTRED IN THE GAP, not a line butted against the data. Padding on both
            # sides is what keeps it from reading as a column: at one-pixel columns a line
            # drawn flush against the block is indistinguishable from a detection.
            if bi < len(ordered) - 1:
                ax.plot([x1 + REGION_GAP / 2] * 2, [b_bot, b_top], color="#c8c8c8",
                        linewidth=0.7, solid_capstyle="butt", zorder=3)
        for i in range(i0, i1 + 1):
            ty = b_top - (i - i0 + 0.5) * RH
            n_sp = collapsed.get(taxa[i], 0)
            if not n_sp:
                ax.text(X0 - 0.09, ty, abbrev(taxa[i]) + (" †" if taxa[i] in spill else ""),
                        ha="right", va="center", fontsize=7.4, style="italic", color=INK)
                continue
            # A collapsed row is "Alternaria (14 spp.)": the genus is a name and sets italic,
            # the count is not. Drawn as two texts and the suffix measured, so the pair still
            # right-aligns on the same edge as every other row label.
            suffix = f" ({n_sp} spp.)"
            t2 = ax.text(X0 - 0.09, ty, suffix, ha="right", va="center", fontsize=7.0,
                         color=MUTED)
            w = t2.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
            ax.text(X0 - 0.09 - w, ty, taxa[i], ha="right", va="center", fontsize=7.4,
                    style="italic", color=INK)
        # No swatch on the band. Colour belongs to the TAXON now, and a band holds taxa of
        # more than one category -- Stem base and head carries a wheat pathogen and a
        # non-pathogen recorded on wheat -- so a single swatch would contradict its own rows.
        ax.text(BANDX, b_top + 0.045, name, ha="left", va="bottom", fontsize=6.8,
                color=MUTED)

    # PATHOGENS vs RESIDENTS, bracketed down the left. Computed from the rows rather than
    # authored: a band is a pathogen band when its taxa are wheat pathogens.
    runs_of = []
    for i, w in enumerate(supers):
        if runs_of and runs_of[-1][2] == w:
            runs_of[-1][1] = i
        else:
            runs_of.append([i, i, w])
    for i0, i1, w in runs_of:
        k0 = next(k for k, (a, b, _) in enumerate(spans) if a <= i0 <= b)
        k1 = next(k for k, (a, b, _) in enumerate(spans) if a <= i1 <= b)
        s_top = top - (i0 * RH + k0 * BAND_GAP)
        s_bot = top - ((i1 + 1) * RH + k1 * BAND_GAP)
        ax.plot([SUPX, SUPX], [s_bot, s_top], color="#9a9a9a", linewidth=1.1,
                solid_capstyle="butt")
        ax.text(SUPX - 0.08, (s_bot + s_top) / 2, SUPER[w], rotation=90, ha="center",
                va="center", fontsize=7.6, color="#6a6a6a")

    # two staggered baselines, with the separation MEASURED against the font
    import math
    need = (6.8 / 72) / math.sin(math.radians(LABEL_ROT))   # inches between same-tier labels
    last = [-1e9, -1e9]
    skipped = []
    for bi, c in enumerate(ordered):
        n = len(bycty[c])
        if n >= LABEL_MIN:
            x = xs[bi] + n * cw / 2
            tier = next((t for t in (0, 1) if x - last[t] >= need), None)
            if tier is None:
                skipped.append(c or "no location")
            else:
                last[tier] = x
                ax.text(x, top + 0.07 + tier * LABEL_TIER, c or "no location",
                        rotation=LABEL_ROT, rotation_mode="anchor", ha="left",
                        va="baseline", fontsize=6.8, color=INK)
    if skipped:
        print(f"  label skipped, no tier clear: {', '.join(skipped)}")
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
        w.writerow(["column", "run", "BioProject", "country", "region", "lat", "lon",
                    "declared_status", "n_calls"])
        for j, rid in enumerate(order):
            m = runs[rid]
            w.writerow([j, rid, m["project"], m["country"], m["region"], m["lat"], m["lon"],
                        m["status"], int((grid[:, j] > 0).sum())])

    state = {k: int((grid == v).sum()) for k, v in
             (("declared", DECL), ("undeclared", UNDECL), ("unknown", UNKNOWN))}
    print(f"wrote {FIG / 'samples.png'}")
    print(f"  {len(order):,} {describe()} runs ({excluded} excluded), "
          f"{len(ordered)} regions, {len(taxa)} taxa drawn")
    print(f"  cells: {state}")
    top = list(dropped.items())[:8]
    print(f"  below MIN_RUNS={MIN_RUNS}: {len(dropped)} taxa not drawn, largest " +
          ", ".join(f"{abbrev(t)} ({n})" for t, n in top))
    small = sorted(((len(v), c) for c, v in bycty.items() if c and len(v) < LABEL_MIN),
                   reverse=True)
    if small:
        print(f"  drawn but unlabelled, under LABEL_MIN={LABEL_MIN}: {len(small)} regions, "
              f"{sum(n for n, _ in small)} runs -- " +
              ", ".join(f"{c} ({n})" for n, c in small[:8]) + " ...")
    noloc = len(bycty.get("", ()))
    print(f"  {noloc} runs carry no location and sit in the unplaced block")
    if spill:
        print("  FLAGGED † congener-tracking, still drawn (never called without a "
              "congener, and at <"
              f"{SPILL_RATIO:.0%} of its coverage): " +
              ", ".join(f"{t} ({n} runs, all with {p}, ratio {r:.2f})"
                        for t, (n, p, r) in sorted(spill.items())))
    if dropped_slot:
        print(f"  not drawn, {sorted(DROP_SLOTS)} slot (no host record of their own): " +
              ", ".join(f"{t} ({n})" for t, n in sorted(dropped_slot.items(),
                                                        key=lambda kv: -kv[1])))
    if collapsed:
        print("  collapsed to genus: " +
              ", ".join(f"{g} ({n} spp.)" for g, n in sorted(collapsed.items())))
    if suspect:
        print("  !! UNCOLLAPSED GENERA BEHAVING LIKE CLUSTERS -- review the collapse column "
              "in data/inputs/pathogen_groups.tsv: " +
              ", ".join(f"{g} (median {n:.0f} spp/run, CV {c:.2f})" for g, n, c in suspect))
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

    # TWO keys, and the split is the point: inside the grid a colour says who knew, in the
    # left margin it says what kind of organism the band holds. The margin colours are
    # 02_reference's, so a category is the same colour in every figure of the chapter.
    x = 0.08
    wp = legible(MATRIX[("wheat", "wheat_path")])
    for colour, label in ((DECLARED_COLOUR, "declared by the study"),
                          (wp, "not declared"),
                          (washed(wp), "declaration not curated"),
                          (None, "   colour = the taxon's category in the reference matrix")):
        if colour is None:
            t = ax.text(x, 0.155, label, ha="left", va="center", fontsize=7.4, color=MUTED)
            x += t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi + 0.10
            continue
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
