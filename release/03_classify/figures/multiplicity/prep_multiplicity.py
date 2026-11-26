#!/usr/bin/env python3
"""Prepare the per-project panel: how many wheat pathogens a run carries.

    python release/03_classify/figures/multiplicity/prep_multiplicity.py

THIS SCRIPT DECIDES; plot_multiplicity.R DRAWS, the same split as figures/samples. Every
bar segment arrives in R as an x range and a hex string, so the stacking order and the
colour scale are decided here and cannot be re-derived differently there.

Reads  ../../data/outputs/supp_data_s6_runs.tsv      n_wheat_pathogen_calls, per run
       ../../data/outputs/supp_data_s5_detections.tsv the calls, for the taxon identities
       ../../data/outputs/arm_resolution.tsv          db2_f_reads, the fungal depth axis
       ../../../01_cohort/data/outputs/supp_data_s2_cohort.tsv  the analysed set and target
       ../samples/samples_cols.tsv                    CHECKED AGAINST, never read as input
Writes multiplicity_bars.tsv      one row per drawn bar segment: x range, colour
       multiplicity_projects.tsv  one row per project: n, rate, depth, label
       multiplicity_ramp.tsv      the ramp and the marker, for the legend

## What this panel is for

Co-infection multiplicity is the chapter's subject and it was in no figure. The per-run
raster shows WHICH pathogens and WHERE; it cannot show how many of them land in one run,
because the pathogens-per-run ramp was removed from it as redundant with the rows it sits
under. This panel is that count, and it is the count per PROJECT.

**A CO-OCCURRENCE NETWORK WAS REJECTED ON MEASUREMENT AND THIS REPLACES IT.** Conditioning
on BioProject reverses the sign of three of the eight commonest pathogen pairs
(P. striiformis with P. triticina: marginal OR 0.28, within-project OR 4.43), so a network
over runs draws where the surveys sampled and reads as a claim that two organisms
associate. See crypt-cooccurrence-confounded and the rejected-approaches note in the
module README. Multiplicity has no such exposure: each project is one observation of its
own rate, and the between-project spread is drawn rather than averaged away.

## The unit is the PROJECT here, and the run in figures/samples

That is deliberate and is the opposite of prep_samples.py's choice, for the opposite
reason. The raster's subject is the pathogen in the sample, and a project spanning twelve
countries averages a worldwide survey into one column. This panel's subject is the SPREAD
between studies: the co-infection rate runs 0% to 85% across the projects with enough runs
to have one, and that heterogeneity is a property of the studies -- what they sampled,
where, at what depth -- not of wheat. Pooling the runs would report a single rate, which
is the number this chapter does not have ([[crypt-prevalence-demoted]]).

## Shares, not counts, and n printed beside every bar

Projects span 4 to 531 analysed runs. Drawn as counts the panel is one project; drawn as
shares every project is comparable and the n carries what each share rests on.

## Two blocks, n >= 10 and n < 10

A 0% rate over 4 runs is not a rate. The three small projects are drawn rather than
dropped -- they are runs the module classified and the cohort holds -- but they sit in a
block of their own so a reader does not read their extremes as the range.

## The depth strip exists to test the obvious objection, and it fails to confirm it

A project with few co-infections might simply have shallow libraries, so the second column
is the median FUNGAL reads per run, which is the depth that matters once the host is
removed (crypt-host-fraction-confound). Over the 22 projects the rank correlation with the
co-infection rate is weak, and four projects with no co-infection at all sit among the
deepest in the cohort, so the spread is not a depth artefact. The strip is drawn because
that is worth showing, not because it explains anything.

## The virus surveys are marked, for the same reason the raster marks them

Five of the six zero-rate projects are virus surveys that declared no fungus because they
were not looking for one. The marker is read from study_target in supp_data_s2, a curated
fact, never inferred from declared_status.
"""

import collections
import csv
import math
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
FIG = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _palette import MATRIX, INK                        # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from _cohort import is_analysed, describe               # noqa: E402
# ramp3 and the location reader come from the raster's own prep, not a second copy: the two
# panels share a page and a ramp drawn by different arithmetic would not look like the same
# figure. prep_samples does nothing at import time.
sys.path.insert(0, str(FIG.parent / "samples"))
from prep_samples import ramp3, country, regions        # noqa: E402

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
S6 = HERE / "data/outputs/supp_data_s6_runs.tsv"
RES = HERE / "data/outputs/arm_resolution.tsv"
S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
COLS = FIG.parent / "samples/samples_cols.tsv"
ARM = "db2"
csv.field_size_limit(10 ** 7)

# 0 is NEUTRAL and 1 to 4+ are the pathogen orange deepening. A single pale-to-dark ramp
# over 0..4 would have made "no wheat pathogen" the palest shade of pathogen; the break
# puts the figure's own statistic -- the share of runs carrying TWO OR MORE -- in the dark
# part of every bar, so the bars are read by how dark their right end is.
NONE_FILL = "#e8e8e8"
RAMP_STOPS = ("#f7c6ae", MATRIX[("wheat_resident", "wheat_path")], "#6b1f08")
TOP = 4                 # the top bin is OPEN: 4 or more, which is 43 runs reaching 6
MIN_N = 10              # below this a project's rate is not a rate; see the module docstring
VIRUS_COLOUR = "#7a7a7a"        # prep_samples.VIRUS_COLOUR, the declared grey


def read(path, comment="#"):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith(comment)),
                                   delimiter="\t"))


def fill(k):
    return NONE_FILL if k == 0 else ramp3(RAMP_STOPS, (k - 1) / (TOP - 1))


def analysed():
    """The runs this panel counts: the same set the raster draws, and asserted equal to it.

    The narrowing is _cohort.is_analysed plus the raster's own location filter. The second
    one has no reason of its own here -- this panel's x axis is a share, not geography --
    and it is applied anyway because the two panels share a page: a reader counting
    pathogen cells down a column of the raster has to be able to find that run in this
    distribution. The cost is 32 runs with no location and it is stated in the caption.

    samples_cols.tsv is the raster's own output and is CHECKED, not read. A check keeps
    each prep script deriving from data/outputs/ while making drift impossible to miss.
    """
    reg = regions()
    runs = {}
    for r in read(S2):
        if not is_analysed(r):
            continue
        if not reg.get(country(r["location"]), ""):
            continue
        runs[r["run"]] = r
    drawn = {r["run"] for r in read(COLS)}
    if drawn and set(runs) != drawn:
        raise SystemExit(
            f"this panel would count {len(runs)} runs and the raster drew {len(drawn)}; "
            f"{len(set(runs) - drawn)} here only, {len(drawn - set(runs))} there only. The "
            f"two panels sit on one page and must hold the same runs. Re-run "
            f"prep_samples.py, or reconcile the narrowing in _cohort.py.")
    return runs


def main():
    runs = analysed()
    s6 = {r["run"]: r for r in read(S6) if r["arm"] == ARM}
    depth = {r["run"]: int(r[f"{ARM}_f_reads"]) for r in read(RES)}
    missing = sorted(set(runs) - set(s6))
    if missing:
        raise SystemExit(f"{len(missing)} analysed runs are not in {S6.name} for arm {ARM}: "
                         f"{', '.join(missing[:5])}")

    n_of = {rid: int(s6[rid]["n_wheat_pathogen_calls"]) for rid in runs}
    by = collections.defaultdict(list)
    for rid, m in runs.items():
        by[m["BioProject"]].append(rid)

    proj = []
    for p, ids in by.items():
        counts = [n_of[i] for i in ids]
        f = [depth[i] for i in ids if i in depth]
        proj.append({
            "project": p, "n": len(ids),
            # THE RATE IS THE SHARE OF RUNS WITH TWO OR MORE, which is what co-infection
            # means here: a run holding its declared pathogen and nothing else is not one.
            "rate": sum(1 for c in counts if c >= 2) / len(ids),
            "mean": statistics.mean(counts), "max": max(counts),
            "f_reads": statistics.median(f) if f else "",
            # log10 of the depth, written here so the panel needs no log scale: the two
            # columns of the figure share one continuous x transform and only their ranges
            # differ, so a log scale on the second would take the first one with it
            "log_f": math.log10(statistics.median(f)) if f else "",
            # a curated fact, read off the study's own statement of aim
            "virus": int(any(runs[i].get("study_target") == "virus" for i in ids)),
            # the block INDEX orders the facets; the block NAME is only ever drawn. A
            # string sort put "n < 10" first, '<' sorting below the '≥' glyph.
            "block_i": 0 if len(ids) >= MIN_N else 1,
            "block": f"n \u2265 {MIN_N}" if len(ids) >= MIN_N else f"n < {MIN_N}"})

    # Order: by rate within block, biggest first, and n breaks a tie so the zero-rate
    # projects are ordered by how much evidence each one's zero rests on.
    proj.sort(key=lambda r: (r["block_i"], -r["rate"], -r["n"]))
    for i, r in enumerate(proj):
        r["row"] = i

    bars = []
    for r in proj:
        counts = collections.Counter(min(n_of[i], TOP) for i in by[r["project"]])
        x = 0.0
        for k in range(TOP + 1):
            share = counts.get(k, 0) / r["n"]
            if share:
                bars.append({"project": r["project"], "row": r["row"], "block": r["block"],
                             "k": k, "klabel": f"{TOP}+" if k == TOP else str(k),
                             "n_runs": counts[k], "share": f"{share:.6f}",
                             # the x range is computed HERE, so R stacks nothing and no
                             # factor-level order can reverse a bar
                             "x0": f"{x:.6f}", "x1": f"{x + share:.6f}", "fill": fill(k)})
            x += share

    def dump(path, rows):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        return path

    dump(FIG / "multiplicity_bars.tsv", bars)
    dump(FIG / "multiplicity_projects.tsv",
         [{"row": r["row"], "project": r["project"], "block": r["block"],
           "block_i": r["block_i"], "n": r["n"],
           "rate": f"{r['rate']:.6f}", "mean": f"{r['mean']:.4f}", "max": r["max"],
           "f_reads": r["f_reads"], "log_f": r["log_f"], "virus": r["virus"],
           # the accession is DRAWN in this colour: a virus survey's is the declared grey,
           # which is the colour the raster marks the same projects in
           "label_colour": VIRUS_COLOUR if r["virus"] else INK} for r in proj])
    # THE FIVE BINS, written whether or not each occurs, the same rule the runs legend
    # follows: the legend states the scheme, and a bin that is empty in this cohort is a
    # fact about the cohort rather than something to leave the reader to infer.
    dump(FIG / "multiplicity_ramp.tsv",
         [{"k": k, "klabel": f"{TOP}+" if k == TOP else str(k), "fill": fill(k),
           "n_runs": sum(1 for rid in runs if min(n_of[rid], TOP) == k)}
          for k in range(TOP + 1)]
         + [{"k": "virus", "klabel": "virus survey, no fungus declared",
             "fill": VIRUS_COLOUR,
             "n_runs": sum(1 for rid in runs
                           if runs[rid].get("study_target") == "virus")}])

    dist = collections.Counter(min(n_of[i], TOP) for i in runs)
    ge2 = sum(v for k, v in dist.items() if k >= 2)
    big = [r for r in proj if r["n"] >= MIN_N]
    taxa = collections.Counter(
        r["taxon"] for r in read(S5)
        if r["arm"] == ARM and r["run"] in runs and r["taxon_class"] == "wheat_pathogen")
    print(f"wrote {FIG}/multiplicity_{{bars,projects,ramp}}.tsv")
    print(f"  {len(runs):,} {describe()} runs over {len(proj)} BioProjects, "
          f"{len(taxa)} wheat pathogens, {sum(taxa.values()):,} calls")
    print("  pathogens per run: " + ", ".join(
        f"{k if k < TOP else str(TOP) + '+'} {dist.get(k, 0)} ({dist.get(k, 0) / len(runs):.1%})"
        for k in range(TOP + 1)))
    print(f"  two or more: {ge2:,} runs ({ge2 / len(runs):.1%})")
    print(f"  {len(big)} projects with n >= {MIN_N}: rate {min(r['rate'] for r in big):.1%} "
          f"to {max(r['rate'] for r in big):.1%}")
    for r in proj:
        print(f"    {r['project']:13s} {r['block']:7s} n={r['n']:4d} rate={r['rate']:6.1%} "
              f"mean={r['mean']:.2f} max={r['max']} "
              f"f_reads={r['f_reads']:>12,.0f}" + ("  virus survey" if r["virus"] else ""))


if __name__ == "__main__":
    main()
