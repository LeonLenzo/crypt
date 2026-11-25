#!/usr/bin/env python3
"""Prepare the per-run panel: every run a column, ordered west to east.

    python release/03_classify/figures/samples/prep_samples.py

THIS SCRIPT DECIDES; plot_samples.R DRAWS. The split is deliberate. Every judgement here --
which runs are analysed, whether a call was declared, which genus collapses, which band and
which colour a taxon takes -- is imported from the pipeline that made it: `declaration` from
detect.py so there is ONE definition of declared, `is_analysed` from _cohort.py, `MATRIX`
from 02_reference. Re-implementing any of that in R would reintroduce exactly the drift
those imports exist to prevent, and a 195-call disagreement with supp_data_s5 is how we
found out it matters. What R gets is a tidy table of cells and the facet structure; what it
does is the layout, which ggplot's faceting does far better than hand-computed inches.

Reads  ../../data/outputs/supp_data_s5_detections.tsv    the calls
       ../../../01_cohort/data/outputs/supp_data_s2_cohort.tsv   runs, location, declaration
       ../../../01_cohort/figures/country_centroids.tsv   to place a country-only run
Writes samples_cells.tsv   one row per DRAWN cell: column, row, state, colour
       samples_rows.tsv    one row per taxon row: label, band, half, colour
       samples_cols.tsv    one row per run column: region, country, project

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

import numpy as np

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _palette import MATRIX, MATRIX_NAME, INK, MUTED, darken  # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from _cohort import is_analysed, describe               # noqa: E402
from detect import declaration, declared_sets, normalise                          # noqa: E402

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
CEN = HERE.parent / "01_cohort/figures/country_centroids.tsv"
GRP = HERE / "data/inputs/pathogen_groups.tsv"
REGIONS = HERE / "data/inputs/country_regions.tsv"
S3REF = HERE.parent / "02_reference/data/outputs/supp_data_s3_reference.tsv"
FIG = Path(__file__).resolve().parent
csv.field_size_limit(10 ** 7)

NONE, DECL, UNDECL, UNKNOWN = 0, 1, 2, 3
# Lightened from #3d3d3d, Leon's call: against a white grid a near-black bar read as the
# heaviest thing on the page, so the declared blocks dominated the figure and the orange
# undeclared calls, which are the result, read as incidental marks on top of them. A mid
# grey still separates unambiguously from every matrix hue without taking the eye first.
DECLARED_COLOUR = "#7a7a7a"     # the study reported it; the one state that is not a category

# A CELL TAKES THE COLOUR THE REFERENCE MATRIX GAVE THAT TAXON. `MATRIX` is keyed on the two
# facts S3 records, where USDA recorded the taxon and what PHI-base scores it as, and is
# imported from 02_reference rather than copied so the two figures cannot drift apart.
#
# The matrix was drawn as inch-wide pills with dark borders; here a cell is one pixel wide,
# and two of its colours (#f0f0f0 for a genus-only taxon, #c3e8bf for a plant resident) vanish
# against white at that size. LUM_CEILING darkens any colour past it until it reaches the
# ceiling, which preserves the hue and the ordering and only changes what cannot be seen.
# The BAND SWATCH keeps the true matrix colour, so the correspondence is still legible.
# S3 spells the PHI-base axis in full; MATRIX is keyed on plot_matrix.py's short names. The
# two have to be translated, not assumed equal -- the panel raised a KeyError on the first
# wheat pathogen it met, which is the good failure mode for a mismatch like this.
PHI_KEY = {"none": "none", "plant_pathogen": "vagrant", "wheat_pathogen": "wheat_path"}

# What each matrix slot is called in the legend. The slot is the figure's hue, so the legend
# has to name it in words rather than leave the reader to infer four greens apart.

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
# A NON-PATHOGEN IS BANDED ON on_wheat, AND NOT ON THE USDA ROW. The USDA row looks like the
# natural axis and is the wrong one. USDA Fungus-Host records fungi OBSERVED ON PLANTS, not
# pathogens, and `route = usda_genus` means only that USDA holds no record under that
# SPECIES name, the taxon having entered because its genus does. It therefore mixes two
# different things and cannot tell them apart: a genuine plant associate USDA has not
# recorded, and a species of a plant-associated genus that lives nowhere near a plant. No
# column in supp_data_s3 separates them.
#
# 02_reference already settled this and the figure should not relitigate it. Its scope is
# read at GENUS precisely because species rank loses 903 of 5,973 fungal detections on field
# wheat, "led by Dioszegia hungarica" -- the taxon an earlier version of this panel dropped
# for carrying the label when it was still called `congener`, which is why it is now called
# `usda_genus`. Reading it as weak evidence would re-impose the species-rank reading the
# reference was built to avoid.
#
# So: two resident bands, split on whether USDA records the SPECIES on wheat. Colour still
# carries the full matrix slot, so the USDA-coverage distinction is visible as a shade
# without being promoted to a claim about the organism.
# Pathogens of other plants are not drawn (see DETAIL_CLASSES), so nothing reaches this
# band now. It is kept rather than deleted because `band` is total over the drawn classes
# and a KeyError is the right failure if that ever changes without this being reconsidered.
OTHER_HOST_BAND = (90, "Other hosts")

# The two halves of the panel, bracketed on the far left. Every drawn taxon row is now a
# pathogen -- residents have no rows -- so the second half holds only the summary ramp.
SUPER = {True: "Pathogens", False: "Pathogens"}
# The resident ramp, pale -> the matrix green -> dark. See ramp3 for why three.
RAMP_STOPS = ("#f4f9f3", "#2e9e52", "#00441b")
RAMP_STOPS_T = lambda t: ramp3(RAMP_STOPS, t)
# The study-target marker takes the DECLARED colour, Leon's call, because it is the same
# kind of fact: both rows say what the study reported. Grey already means "the study knew",
# and a none_declared virus survey is the limiting case of that -- it reported no fungus at
# all -- so giving it a hue of its own would have implied a third category of evidence.
VIRUS_COLOUR = DECLARED_COLOUR
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


def ramp3(stops, t):
    """Interpolate a LIST of hex stops at t in [0,1]. Piecewise-linear, evenly spaced.

    The resident ramp needs three stops, not two. Its counts are strongly bimodal -- 47% of
    runs call no resident at all and 60% call at most one, then a second cluster sits at 13
    to 20 -- and a two-stop pale-to-green ramp ran from near-white to a MID green, so the
    whole scale lived in the top half of the lightness range and neither mode was legible
    against the other. Three stops let the ramp finish genuinely dark while keeping the
    palette's own green as the midpoint, so the hue is still the matrix hue.
    """
    if t <= 0:
        return stops[0]
    if t >= 1:
        return stops[-1]
    n = len(stops) - 1
    i = min(int(t * n), n - 1)
    return ramp(stops[i], stops[i + 1], t * n - i)


def ramp(pale, dark, t):
    """Interpolate pale->dark at t in [0,1], returned as hex."""
    a = [int(pale.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4)]
    b = [int(dark.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4)]
    return "#%02x%02x%02x" % tuple(round(x + (y - x) * t) for x, y in zip(a, b))


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
# EVERY DETECTED PATHOGEN GETS A ROW. Leon's call: the rare ones are the interesting ones,
# and a panel that draws nine of fourteen invites "nine were detected". The threshold was
# never evidential, it was a height budget, and the residents paid that budget instead.
MIN_RUNS = 1            # a taxon needs this many calls across the cohort to get a row
# The two DIAGNOSTICS below keep a real threshold, for a different reason: both are medians
# over the runs a taxon shares with a congener, and a median over two runs is not a measure.
# They decide what the script WARNS about, never what it draws.
MIN_FLAG_RUNS = 10      # below this a congener-spill ratio or a collapse CV is not a number
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

# The host is excluded: every library has one, no study declares it as a pathogen, and a row
# that is charcoal-free and full in every column says nothing. `plant_pathogen` IS drawn even
# though it is nearly empty here -- 5 taxa and 12 calls across the analysed set, none reaching
# MIN_RUNS -- because its emptiness is a result of db2's artefact removal and is reported in
# the caption rather than hidden by never looking.
DRAWN_CLASSES = {"wheat_pathogen", "plant_pathogen", "resident"}

# RESIDENTS GET NO ROW OF THEIR OWN; the per-run ramp at the foot of the panel summarises
# them. Thirteen resident rows were half the panel's height, and the four detection-space
# panels have to share an A4 page with it. Nothing is judged false and nothing leaves the
# data: every resident call is still in supp_data_s5, and the ramp now counts them MORE
# completely than the rows did, because it no longer stops at the taxa clearing MIN_RUNS.
# WHEAT PATHOGENS ONLY. The five pathogens of other plants that MIN_RUNS = 1 brought into
# the panel are out again, Leon's call: each is one to six calls sitting at coverage 0.023
# to 0.045, which is the band the pure-fungal control showed its own artefact reaching
# (<= 0.0336), and four of the five are congeners or close relatives of a wheat pathogen
# that IS present. They are too tenuous to draw as presence and release/04_align is what
# resolves them, from species-unique reads. They stay in supp_data_s5 and the caption names
# them with their coverages.
DETAIL_CLASSES = {"wheat_pathogen"}



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
    s3 = {r["taxon"]: (r["usda_row"], r["phibase_col"], r["on_wheat"]) for r in
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
                "lat": lat, "lon": lon, "status": r["declared_status"],
                "target": r.get("study_target", "")}
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
    # RUNS WITH NO LOCATION ARE DROPPED. The panel's x axis is geography, and a block of
    # columns with no place on it is a column of data the axis cannot describe; it was
    # sitting at the right edge labelled "no location", which reads as a region. They stay
    # in supp_data_s5 and in every count the module makes. Leon's call, 2026-11-23.
    noloc = [k for k, m in runs.items() if not m["region"]]
    for k in noloc:
        del runs[k]

    bycty = defaultdict(list)
    for rid, m in runs.items():
        bycty[m["region"]].append(rid)
    # The gazetteer check is on COUNTRIES, which is where the archive's spelling errors are
    # ('Zimbabwae' was one). Groups are regions now, and a region is placed by its own runs'
    # longitudes, so the centroid fallback only matters if a whole region had no coordinates.
    cc = collections.Counter(m["country"] for m in runs.values() if m["country"])
    unknown_name = sorted((c, n) for c, n in cc.items() if c not in cen)
    pos, unplaced = {}, []
    for c, ids in bycty.items():
        lons = [runs[i]["lon"] for i in ids if runs[i]["lon"] is not None]
        if lons:
            pos[c] = sum(lons) / len(lons)
        else:
            unplaced.append(c)
    ordered = sorted(pos, key=lambda c: pos[c]) + sorted(unplaced)

    # ── within a region: BLOCK BY WHAT THE STUDY DECLARED ─────────────────────
    # Leon's call. Geography alone scattered each survey's runs through its region, so the
    # one comparison the panel exists to support -- what a study looked for against what its
    # reads carry -- had to be made across non-adjacent columns. Blocking by declaration puts
    # every run that was looking for the same organism side by side, and the declared row
    # becomes a solid charcoal bar with the undeclared calls legible against it.
    #
    # THE BLOCKS ARE ORDERED BY THE DISEASE BANDS, not by size, so the blocks run down the
    # rows in the same order the rows run down the panel and the diagonal reads. A
    # declaration naming several organisms takes the band of its FIRST, and the undeclared
    # block sorts last wherever it occurs.
    #
    # Longitude still orders runs INSIDE a block, so the geography is not lost, it is nested
    # one level down. The accession is the final tie-break and is what makes the column order
    # reproducible for runs sharing a site.
    grp = groups()
    band_of = {}
    for rid, m in runs.items():
        spp, gen = declared_sets(decl_of[rid].get("declared_pathogens", ""))
        hits = [grp[t][0] for t in grp
                if normalise(t) in spp or t.split()[0].lower() in gen]
        # no declaration, or one naming nothing we band: last, but kept as its own block
        band_of[rid] = (min(hits) if hits else 10 ** 6,
                        decl_of[rid].get("declared_pathogens", "") or "\uffff")
    key = lambda i: (band_of[i], runs[i]["lon"] is None, runs[i]["lon"] or 0.0,
                     runs[i]["lat"] or 0.0, i)
    order = [i for c in ordered for i in sorted(bycty[c], key=key)]

    # `declared` is RECOMPUTED here with detect.py's own function rather than read off S5.
    # One definition, imported by the producer and by this figure, so the two cannot drift.
    # It is checked against the stored column and any disagreement is reported loudly: a
    # difference means S5 was written before detect.py's congener fix (2026-11-23) and the
    # published supplementary table is stale until detect.py is re-run on Setonix.
    hit, declared_in, mismatched = defaultdict(set), defaultdict(set), 0
    cluster, cvs, collapsed_of = defaultdict(set), defaultdict(list), defaultdict(set)
    klass, cov, gen_of, class_calls = {}, defaultdict(dict), {}, {}
    slot_votes = defaultdict(collections.Counter)
    with open(S5) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["run"] not in runs or r["taxon_class"] not in DRAWN_CLASSES:
                continue
            klass[r["taxon"]] = r["taxon_class"]
            class_calls[(r["run"], r["taxon_class"])] = \
                class_calls.get((r["run"], r["taxon_class"]), 0) + 1
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
            slot_votes[key][s3.get(r["taxon"], ("plant_resident", "none", "no"))] += 1
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
        if med_n >= 3 and med_cv < 0.35 and not collapsed and len(hit.get(g, ())) >= MIN_FLAG_RUNS:
            suspect.append((g, med_n, med_cv))

    # A collapsed genus has no row in S3, and taking the slot of whichever member happened to
    # be read last made its band depend on file order. Take the MODAL slot of the members,
    # the slot name breaking a tie, so the result is reproducible.
    modal = {t: min(v.most_common(), key=lambda kv: (-kv[1], kv[0]))[0]
             for t, v in slot_votes.items()}
    mrow = {t: (k[0], PHI_KEY[k[1]]) for t, k in modal.items()}
    onwheat = {t: k[2] for t, k in modal.items()}

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
        if len(rr) < MIN_FLAG_RUNS:
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
    band = lambda t: grp[t][:2] if is_wp(t) else OTHER_HOST_BAND
    taxa = sorted((t for t in hit if len(hit[t]) >= MIN_RUNS
                   and klass.get(t) in DETAIL_CLASSES),
                  key=lambda t: (band(t)[0], -len(hit[t]), t))
    # WHAT IS STILL NOT DRAWN, and it is no longer the run count. With MIN_RUNS at 1 the
    # only pathogen that can miss a row is one COLLAPSED INTO A GENUS whose other members
    # are residents: Alternaria alternata is a wheat pathogen, the other nineteen Alternaria
    # called here are not, and the genus is drawn as one resident row because the species
    # calls inside it are one unresolved organism. Reporting "0 taxa not drawn" while that
    # happens is the silent kind of wrong, so it is counted and named.
    # `gen_of` is keyed by the COLLAPSED key, so a collapsed species is not in it at all;
    # the collapse flag in pathogen_groups is what says where the species went.
    absorbed = {t: (t.split()[0], len(hit[t.split()[0]])) for t in klass
                if klass[t] in DETAIL_CLASSES and grp.get(t, (0, "", False))[2]
                and t.split()[0] in hit and t.split()[0] not in taxa}
    dropped = {t: len(v) for t, v in sorted(hit.items(), key=lambda kv: (-len(kv[1]), kv[0]))
               if len(v) < MIN_RUNS}

    idx = {rid: j for j, rid in enumerate(order)}
    grid = np.zeros((len(taxa), len(order)), dtype=np.uint8)
    for i, t in enumerate(taxa):
        for rid in hit[t]:
            grid[i, idx[rid]] = (DECL if rid in declared_in[t]
                                 else UNKNOWN if runs[rid]["status"] == "not_curated"
                                 else UNDECL)
    return (order, taxa, grid, runs, ordered, bycty, (dropped, absorbed),
            sorted(unknown_name),
            excluded, [band(t)[1] for t in taxa],
            [is_wp(t) for t in taxa], mismatched, class_calls,
            {t: len(collapsed_of[t]) for t in taxa if len(collapsed_of[t]) > 1}, suspect,
            # a band's category is the modal class of the taxa in it, read off the data
            # rather than authored a second time in the grouping table
            {t: collections.Counter(
                klass[x] for x in (collapsed_of[t] or {t})).most_common(1)[0][0]
             for t in taxa}, {t: mrow[t] for t in taxa}, spill, len(noloc))




def main():
    (order, taxa, grid, runs, ordered, bycty, (dropped, absorbed), unknown_name, excluded,
     bands, supers, mismatched, class_calls, collapsed, suspect, tclass, mrow, spill,
     noloc) = load()

    rows = []
    for i, t in enumerate(taxa):
        base = legible(MATRIX[mrow[t]])
        rows.append({"row": i, "taxon": t,
                     "label": abbrev(t) + (" \u2020" if t in spill else ""),
                     "n_collapsed": collapsed.get(t, 0),
                     "band": bands[i], "half": SUPER[supers[i]],
                     "colour": base, "colour_unknown": washed(base),
                     "italic": 1,
                     "category": MATRIX_NAME[mrow[t]],
                     "n_runs": int((grid[i] > 0).sum())})
    with open(FIG / "samples_rows.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    with open(FIG / "samples_cols.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["col", "run", "BioProject", "country", "region", "lat", "lon",
                    "declared_status", "n_calls"])
        for j, rid in enumerate(order):
            m = runs[rid]
            w.writerow([j, rid, m["project"], m["country"], m["region"], m["lat"], m["lon"],
                        m["status"], int((grid[:, j] > 0).sum())])

    # ONLY NON-ZERO CELLS. The full product is 1,347 x 22 and 93% of it is absence; writing
    # absence as rows would make a 30,000-line table whose only content is that nothing was
    # found there. ggplot draws the panel background instead.
    name = {DECL: "declared", UNDECL: "undeclared", UNKNOWN: "unknown"}
    n_cell = 0
    with open(FIG / "samples_cells.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["col", "row", "state", "fill"])
        for i in range(grid.shape[0]):
            base = legible(MATRIX[mrow[taxa[i]]])
            fill = {DECL: DECLARED_COLOUR, UNDECL: base, UNKNOWN: washed(base)}
            for j in np.nonzero(grid[i])[0]:
                st = int(grid[i, j])
                w.writerow([int(j), i, name[st], fill[st]])
                n_cell += 1

    # ── per-column total, ONE summary row: residents ──────────────────────────
    # A count per RUN cannot be a number on a one-pixel column, so it is a colour.
    #
    # ONLY RESIDENTS GET ONE, Leon's call. The residents have no rows, so without the ramp
    # they are not in the figure at all. The pathogens have nine, so a pathogens-per-run
    # ramp restates a column the reader can already count by eye, and restating it in a
    # third encoding invites the reader to check it against the rows and find a discrepancy:
    # the ramp counted all 14 detected wheat pathogens plus the 5 plant pathogens, the rows
    # draw the 9 above MIN_RUNS, so the two disagree BY DESIGN on 5 taxa and 19 runs.
    #
    # THE COST IS STATED, not hidden: the five wheat pathogens below MIN_RUNS
    # (F. avenaceum 6 runs, O. yallundae 5, F. culmorum 3, B. sorokiniana 3, A. alternata 2)
    # and the five non-wheat plant pathogens now appear nowhere in the panel. They are in
    # supp_data_s5 and the caption names them.
    #
    # THE RAMP COUNTS EVERY RESIDENT CALL IN THE RUN, drawn or not -- which is all of them,
    # residents having no rows.
    halves = [("Residents", MATRIX[("wheat_resident", "none")])]
    want_of = {"Residents": {"resident"}}
    idx_of = {rid: j for j, rid in enumerate(order)}
    per_half = {}
    for hname, hue in halves:
        counts = np.zeros(grid.shape[1], dtype=int)
        for (rid, kl), n in class_calls.items():
            if kl in want_of[hname] and rid in idx_of:
                counts[idx_of[rid]] += n
        per_half[hname] = counts

    extra_rows, extra_cells = [], []
    for hname, hue in halves:
        counts = per_half[hname]
        hi = int(counts.max())
        r = len(rows) + len(extra_rows)
        # The summary row sits at the BOTTOM, in a half of its own, and its label box is NOT
        # tinted: the tint on the taxon boxes says which matrix category the rows are, and a
        # summary of several categories has none to claim.
        extra_rows.append({"row": r, "taxon": f"__total_{hname.lower()}__",
                           # The name goes in the BAND column, not the taxon column, so it
                           # lines up with Rusts / Mildews / Rots rather than with the
                           # binomials. The row itself is unnamed: it is one row, and the
                           # band name beside it already says what it is.
                           #
                           # THE BLANK MUST STILL BE UNIQUE. `label` is the discrete y
                           # scale, so two rows both labelled " " are one factor level and
                           # ggplot stops with "factor level is duplicated". Each unnamed
                           # row gets its own width of blank.
                           "label": " " * (len(extra_rows) + 1),
                           "n_collapsed": 0, "band": hname, "half": "Total",
                           "colour": hue, "colour_unknown": washed(hue), "italic": 0,
                           "category": f"{hname.lower()} per run, 0\u2013{hi}",
                           "n_runs": int((counts > 0).sum())})
        for j, k in enumerate(counts):
            if k:
                extra_cells.append((j, r, "total", RAMP_STOPS_T(int(k) / hi)))
    # ── the study-target marker ───────────────────────────────────────────────
    # Nine of the thirty-one cohort projects are VIRUS SURVEYS, and they are the whole
    # explanation of the panel's largest unexplained feature: a block of columns, mostly in
    # the Americas, with no declared pathogen and often no pathogen call either. Those
    # studies declare no fungus because they were not looking for one, and their libraries
    # are shallower for the same reason. Without the marker a reader sees an absence and
    # has to guess whether it is missing curation, and 'uncurated' is exactly the reading
    # the chapter must not invite.
    #
    # IT IS A CURATED FACT, read from study_target in supp_data_s2 and ultimately from each
    # paper's own statement of aim, NOT inferred from declared_status being none_declared.
    # The two coincide in this cohort and that coincidence is a finding, not a definition.
    # ROW -1, so it sorts ABOVE the pathogens: every ordering in plot_samples.R is taken
    # from this column, so a negative index is all that placing it first requires.
    vrow = -1
    vcols = [j for j, rid in enumerate(order)
             if runs[rid].get("target") == "virus"]
    if vcols:
        extra_rows.append({"row": vrow, "taxon": "__study_target__",
                           "label": " " * (len(extra_rows) + 1),
                           "n_collapsed": 0, "band": "Virus survey", "half": "Study",
                           "colour": VIRUS_COLOUR, "colour_unknown": washed(VIRUS_COLOUR),
                           "italic": 0, "category": "virus survey, no fungus declared",
                           "n_runs": len(vcols)})
        extra_cells.extend((j, vrow, "virus", VIRUS_COLOUR) for j in vcols)

    with open(FIG / "samples_rows.tsv", "a", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writerows(extra_rows)
    with open(FIG / "samples_cells.tsv", "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerows(extra_cells)
    n_cell += len(extra_cells)
    with open(FIG / "samples_ramps.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["half", "k", "max", "fill"])
        for hname, hue in halves:
            hi = int(per_half[hname].max())
            for k in range(hi + 1):
                w.writerow([hname, k, hi,
                            RAMP_STOPS_T(k / hi)])
    print("  per-column totals: " + ", ".join(
        f"{n} 0-{int(per_half[n].max())}" for n, _ in halves))

    print(f"wrote {FIG}/samples_{{cells,rows,cols}}.tsv")
    print(f"  {len(order):,} {describe()} runs ({excluded} excluded, {noloc} with no "
          f"location), {len(ordered)} regions, {len(taxa)} rows, {n_cell:,} cells")
    print(f"  cells: " + ", ".join(
        f"{v} {sum(1 for i in range(grid.shape[0]) for _ in np.nonzero(grid[i] == k)[0])}"
        for k, v in name.items()))
    if spill:
        print("  FLAGGED \u2020 congener-tracking: " +
              ", ".join(f"{t} ({n} runs, all with {p}, ratio {r:.2f})"
                        for t, (n, p, r) in sorted(spill.items())))
    if collapsed:
        print("  collapsed to genus: " +
              ", ".join(f"{g} ({n} spp.)" for g, n in sorted(collapsed.items())))
    if suspect:
        print("  !! UNCOLLAPSED GENERA BEHAVING LIKE CLUSTERS: " +
              ", ".join(f"{g} (median {n:.0f} spp/run, CV {c:.2f})" for g, n, c in suspect))
    if mismatched:
        print(f"  !! {mismatched} calls disagree with S5's stored `declared` column")
    if unknown_name:
        print("  COUNTRY NAMES THE GAZETTEER CANNOT RESOLVE: " +
              ", ".join(f"{c!r} ({n})" for c, n in unknown_name))
    if dropped:
        top = list(dropped.items())[:8]
        print(f"  below MIN_RUNS={MIN_RUNS}: {len(dropped)} taxa not drawn, largest " +
              ", ".join(f"{abbrev(t)} ({n})" for t, n in top))
    if absorbed:
        print(f"  NOT DRAWN, absorbed into a collapsed genus that is not a pathogen row: "
              + ", ".join(f"{abbrev(t)} -> {g} ({n} runs)"
                          for t, (g, n) in sorted(absorbed.items())))


if __name__ == "__main__":
    main()
