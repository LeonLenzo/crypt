#!/usr/bin/env python3
"""caption.py — the caption for the declared-against-recovered panel, written from the data.

    python release/03_classify/figures/declared/caption.py

Reads  declared_cells.tsv                              (beside this script, written by the plot)
       ../../data/outputs/supp_data_s5_detections.tsv   for the taxa the panel does not draw
       ../../../01_cohort/data/outputs/supp_data_s2_cohort.tsv   the cohort and its curation
Writes declared_caption.txt

The cell counts come from `declared_cells.tsv`, which the plot writes as it draws, so the
caption cannot narrate a different figure. The two other reads are for things that table
cannot carry: how much of the cohort the panel leaves out, and the rows dropped under the
plot's MIN_RUNS.

Nothing here is prose-with-numbers-typed-in: the narrative is a template and every figure in
it is read at write time. Rerun after any pipeline change.
"""

import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
FIG = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _caption import rd, para, write                    # noqa: E402

sys.path.insert(0, str(FIG))
from plot_declared import MIN_RUNS, abbrev              # noqa: E402

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"


def numbers():
    cells = rd(FIG / "declared_cells.tsv")
    cohort = rd(S2)
    calls = rd(S5)

    n = {}
    n["cohort_runs"] = len(cohort)
    n["cohort_projects"] = len({r["BioProject"] for r in cohort})
    drawn_p = {r["BioProject"] for r in cells}
    n["projects"] = len(drawn_p)
    n["runs"] = sum({r["BioProject"]: int(r["n_runs"]) for r in cells}.values())
    # The excluded projects are two different things and the caption must not merge them:
    # `not_curated` means nobody has read the paper, `none_declared` means somebody read it
    # and it named no pathogen. Only the first makes a detection unknown.
    n["nc_projects"] = len({r["BioProject"] for r in cohort
                            if r["declared_status"] == "not_curated"})
    n["nc_runs"] = sum(1 for r in cohort if r["declared_status"] == "not_curated")
    n["nd_projects"] = len({r["BioProject"] for r in cohort
                            if r["declared_status"] == "none_declared"})
    n["nd_runs"] = sum(1 for r in cohort if r["declared_status"] == "none_declared")
    nd_p = {r["BioProject"] for r in cohort if r["declared_status"] == "none_declared"}
    n["nd_calls"] = len({r["run"] for r in calls
                         if r["BioProject"] in nd_p and r["is_wheat_pathogen"] == "yes"})

    # the declared cells: the positive control
    dec = [r for r in cells if r["declared"] == "1"]
    shares = sorted(float(r["share"]) for r in dec)
    n["n_declared"] = len(dec)
    n["dec_at_100"] = sum(1 for s in shares if s >= 0.995)
    n["dec_min"] = shares[0]
    n["dec_min_cell"] = next(f'{abbrev(r["taxon"])} in {r["BioProject"]}'
                             for r in dec if float(r["share"]) == shares[0])
    # NOT summed over pairs: PRJNA718488 declares two targets, so summing its runs once
    # per pair would inflate the denominator by its whole size. Recovery is reported per
    # PAIR, which is the unit the outlines draw.
    n["dec_median"] = shares[len(shares) // 2]

    # the undeclared cells: the result, quoted as a SPREAD over projects, never pooled
    und = [r for r in cells if r["declared"] == "0" and int(r["n_runs_with_call"])]
    n["n_undeclared_cells"] = len(und)
    n["undeclared_taxa"] = len({r["taxon"] for r in und})
    by_tax = defaultdict(list)
    for r in und:
        by_tax[r["taxon"]].append(r)
    lead = max(by_tax, key=lambda t: sum(int(r["n_runs_with_call"]) for r in by_tax[t]))
    lsh = sorted(float(r["share"]) for r in by_tax[lead])
    n["lead"] = abbrev(lead)
    n["lead_runs"] = sum(int(r["n_runs_with_call"]) for r in by_tax[lead])
    n["lead_projects"] = sum(1 for s in lsh if s > 0)
    n["lead_lo"], n["lead_hi"] = lsh[0], lsh[-1]

    # per project, the share of runs carrying any undeclared wheat pathogen
    anyund = defaultdict(set)
    nruns = {r["BioProject"]: int(r["n_runs"]) for r in cells}
    for r in calls:
        if (r["is_wheat_pathogen"] == "yes" and r["declared"] == "undeclared"
                and r["BioProject"] in drawn_p):
            anyund[r["BioProject"]].add(r["run"])
    pr = sorted((len(anyund[p]) / nruns[p], p) for p in drawn_p)
    n["any_lo"], n["any_lo_p"] = pr[0][0], pr[0][1]
    n["any_hi"], n["any_hi_p"] = pr[-1][0], pr[-1][1]
    n["any_runs"] = sum(len(v) for v in anyund.values())

    # rows the panel does not draw, recomputed here rather than passed through the plot
    small = defaultdict(set)
    for r in calls:
        if r["is_wheat_pathogen"] == "yes" and r["BioProject"] in drawn_p:
            small[r["taxon"]].add(r["run"])
    below = sorted(((len(v), t) for t, v in small.items()
                    if len(v) < MIN_RUNS and not any(
                        c["taxon"] == t and c["declared"] == "1" for c in cells)),
                   reverse=True)
    n["below"] = ", ".join(f"{abbrev(t)} ({c})" for c, t in below)
    n["n_below"] = len(below)
    return n


def main():
    n = numbers()
    pct = lambda x: f"{round(x * 100)}%"

    text = para(
        "What Each Study Declared, Against What Its Reads Contain",
        "Caption for declared.png.",
        "",
        "GENERATED by caption.py — do not hand-edit, the numbers will drift.",
        "",
        "",
        "FULL CAPTION",
        "",
        f"Wheat pathogens detected in {n['runs']:,} RNA-seq runs from "
        f"{n['projects']} BioProjects, one column per project and one row per pathogen. "
        f"Cell shading is the share of that project's runs in which the pathogen was "
        f"called. A cell is drawn only where the pathogen was called at least once, so "
        f"blank means no run in that project reached the detection floors. "
        f"The outlined cell is the pathogen the study itself declared.",
        "",
        f"The outlined cells are the positive control and they are read first. Across "
        f"{n['n_declared']} declared study-pathogen pairs, the median pair has the declared "
        f"organism in {pct(n['dec_median'])} of the project's runs, "
        f"{n['dec_at_100']} of {n['n_declared']} pairs reach 100%, and the lowest is "
        f"{pct(n['dec_min'])} ({n['dec_min_cell']}). No outlined cell is empty: there is no "
        f"project whose declared target the screen failed to find. Nothing was tuned per "
        f"project, and the same floors produced every other cell in the panel.",
        "",
        f"Everything off the outlines is undeclared. {n['undeclared_taxa']} pathogens are "
        f"called in projects that did not declare them, over {n['n_undeclared_cells']} "
        f"study-pathogen pairs and counting the rarer taxa listed below that are not given a "
        f"row, and {n['any_runs']:,} runs carry at least one. The leading "
        f"undeclared organism is {n['lead']}, called in {n['lead_runs']} runs across "
        f"{n['lead_projects']} projects and declared by none of them.",
        "",
        f"The result is the SPREAD across projects, not a pooled rate: the share of a "
        f"project's runs carrying any undeclared wheat pathogen runs from "
        f"{pct(n['any_lo'])} ({n['any_lo_p']}) to {pct(n['any_hi'])} ({n['any_hi_p']}), and "
        f"{n['lead']} alone from {pct(n['lead_lo'])} to {pct(n['lead_hi'])} between projects "
        f"that detect it. The project is the unit of analysis, because crop, protocol and "
        f"sampling design are all nested within it, so a percentage pooled over runs would "
        f"be a weighted average over a handful of independent units.",
        "",
        f"Only BioProjects whose declaration has been curated from the study's own text are "
        f"shown. {n['nc_projects']} further projects of the cohort ({n['nc_runs']} runs) "
        f"have not been curated, so their detections are unknown rather than undeclared and "
        f"an outline would assert something the curation does not support; a further "
        f"{n['nd_projects']} ({n['nd_runs']} runs) was curated and declared no pathogen, "
        f"with a wheat pathogen called in {n['nd_calls']} of its runs. All are in the "
        f"supplementary detection table. "
        f"{n['n_below']} further wheat pathogens were called in fewer than {MIN_RUNS} runs "
        f"of the drawn projects and are not given a row: {n['below']}. Every cell, drawn or "
        f"not, is in declared_cells.tsv.",
        "",
        f"Detection is presence, not activity, and a species-level call inside a "
        f"densely represented genus is a call on that genus' k-mer space; both are taken up "
        f"in 04_align.",
    )
    write(FIG / "declared_caption.txt", text)


if __name__ == "__main__":
    main()
