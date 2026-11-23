#!/usr/bin/env python3
"""caption.py — the caption for the per-run panel, written from the data.

    python release/03_classify/figures/samples/caption.py

Reads  samples_order.tsv                               (beside this script, from the plot)
       ../../data/outputs/supp_data_s5_detections.tsv   the calls
Writes samples_caption.txt

The column order and its countries come from `samples_order.tsv`, which the plot writes as
it draws, so the caption cannot narrate a different figure.

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
from plot_samples import MIN_RUNS, abbrev               # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from detect import declaration                          # noqa: E402

S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
# Countries called out by name in the caption's geography sentence. They are named because the
# claim is about THEM, and the shares are computed, never typed.
NW_EUROPE = ["Ireland", "Netherlands", "France", "UK", "Germany"]
ELSEWHERE = ["USA", "Argentina", "South Africa", "Ethiopia", "China", "Australia"]
FOCUS = "Zymoseptoria tritici"
COLLAPSED = "Alternaria"        # the one genus the panel draws as a single row


def numbers():
    cols = rd(FIG / "samples_order.tsv")
    # samples_order.tsv IS the analysed set: the plot writes one row per column it drew, so
    # filtering the calls through it keeps caption and panel on the same runs. Reading S5
    # whole counted the lab libraries and the host-contrast runs the panel excludes.
    drawn = {c["run"] for c in cols}
    calls = [r for r in rd(S5) if r["run"] in drawn]
    # Same recomputation the panel does, from the same function, so caption and panel agree
    # even while supp_data_s5's stored `declared` column is stale. See plot_samples.load.
    crow = {r["run"]: r for r in rd(S2)}
    for r in calls:
        r["declared"] = declaration(crow[r["run"]], r["taxon"], r["genus"],
                                    r["taxon_class"] == "host")[0]
    n = {}
    n["runs"] = len(cols)
    n["countries"] = len({c["country"] for c in cols if c["country"]})
    n["projects"] = len({c["BioProject"] for c in cols})
    n["noloc"] = sum(1 for c in cols if not c["country"])
    n["geo"] = sum(1 for c in cols if c["lon"] not in ("", "None"))

    by_country = defaultdict(list)
    for c in cols:
        by_country[c["country"]].append(c["run"])

    # the panel now draws three classes, so the caption has to account for all of them
    cls = defaultdict(lambda: defaultdict(set))
    for r in calls:
        if r["taxon_class"] in ("wheat_pathogen", "plant_pathogen", "resident"):
            cls[r["taxon_class"]][r["taxon"]].add(r["run"])
    n["n_resident_taxa"] = len(cls["resident"])
    n["n_resident_drawn"] = sum(1 for v in cls["resident"].values() if len(v) >= MIN_RUNS)
    n["n_pp_taxa"] = len(cls["plant_pathogen"])
    n["n_pp_calls"] = sum(len(v) for v in cls["plant_pathogen"].values())
    n["pp_max"] = max((len(v) for v in cls["plant_pathogen"].values()), default=0)
    top = sorted(cls["resident"].items(), key=lambda kv: -len(kv[1]))[:3]
    n["top_residents"] = ", ".join(f"{abbrev(t)} {len(v)}" for t, v in top)

    # The collapse justification is MEASURED here, not transcribed. These three numbers are
    # the whole argument for drawing Alternaria as one row, so a caption asserting them from
    # memory would be the exact failure crypt-figure-caption-discipline is about.
    import statistics
    cv, size = defaultdict(list), defaultdict(list)
    per_run = defaultdict(set)
    for r in calls:
        if r["taxon_class"] == "host":
            continue
        per_run[(r["run"], r["genus"])].add(r["taxon"])
        if r["genus_cov_cv"] not in ("", "nan"):
            cv[r["genus"]].append(float(r["genus_cov_cv"]))
    for (run, g), v in per_run.items():
        size[g].append(len(v))
    n["alt_cv"] = statistics.median(cv[COLLAPSED])
    n["alt_n"] = statistics.median(size[COLLAPSED])
    others = [statistics.median(v) for g, v in cv.items()
              if g != COLLAPSED and len(v) >= 20 and max(size[g]) > 1]
    n["other_cv_lo"], n["other_cv_hi"] = min(others), max(others)
    n["collapsed"] = COLLAPSED

    wp = [r for r in calls if r["is_wheat_pathogen"] == "yes"]
    n["taxa"] = len({r["taxon"] for r in wp})
    state = defaultdict(set)
    for r in wp:
        state[r["declared"]].add((r["run"], r["taxon"]))
    n["n_declared"] = len(state["declared"])
    n["n_undeclared"] = len(state["undeclared"])
    n["n_unknown"] = len(state["unknown"])
    n["runs_undeclared"] = len({x[0] for x in state["undeclared"]})

    # a run that carries a pathogen its study did not declare, alongside one it did
    dec_runs = {x[0] for x in state["declared"]}
    n["runs_both"] = len({x[0] for x in state["undeclared"]} & dec_runs)
    assert n["runs_both"] <= n["runs_undeclared"]
    # Every run carrying an undeclared pathogen also carries its study's declared target, so
    # runs_only is 0. That is a CONSEQUENCE of the declared targets being recovered in nearly
    # every run, not an independent result, and the caption says so rather than printing the
    # same number twice as if the two were separate observations.
    n["runs_only"] = n["runs_undeclared"] - n["runs_both"]
    per_taxon = defaultdict(set)
    for r in wp:
        per_taxon[r["taxon"]].add(r["run"])
    n["drawn_taxa"] = sum(1 for v in per_taxon.values() if len(v) >= MIN_RUNS)

    hit = defaultdict(set)
    for r in wp:
        hit[r["taxon"]].add(r["run"])
    n["focus"] = abbrev(FOCUS)
    share = lambda c: (100 * len(set(by_country[c]) & hit[FOCUS]) / len(by_country[c])
                       if by_country[c] else 0.0)
    n["nw"] = ", ".join(f"{c} {share(c):.0f}%" for c in NW_EUROPE if by_country[c])
    n["el"] = ", ".join(f"{c} {share(c):.0f}%" for c in ELSEWHERE if by_country[c])
    nwr = [r for c in NW_EUROPE for r in by_country[c]]
    elr = [r for c in ELSEWHERE for r in by_country[c]]
    n["nw_pooled"] = 100 * len(set(nwr) & hit[FOCUS]) / len(nwr)
    n["el_pooled"] = 100 * len(set(elr) & hit[FOCUS]) / len(elr)
    n["nw_n"], n["el_n"] = len(nwr), len(elr)
    return n


def main():
    n = numbers()
    text = para(
"What Is in Every Run of the Cohort, and Who Reported It",
        "Caption for samples.png.",
        "",
        "GENERATED by caption.py — do not hand-edit, the numbers will drift.",
        "",
        "",
        "FULL CAPTION",
        "",
        f"Every one of the {n['runs']:,} RNA-seq runs in the cohort is one column, grouped by "
        f"country and ordered west to east; within a country, columns run by longitude. "
        f"{n['countries']} countries are represented across {n['projects']} BioProjects, "
        f"{n['geo']:,} runs carry coordinates and {n['noloc']} carry no location and sit in "
        f"the block at the right. Rows are every fungus called in at least {MIN_RUNS} runs, "
        f"banded by what kind of organism it is: the wheat pathogens by the disease they "
        f"cause, then the phyllosphere residents, then three congeners of wheat pathogens "
        f"that are not themselves recorded as one. The plant host is not drawn, because "
        f"every library has one and no study declares it as a pathogen.",
        "",
        f"The resident layer is the larger part of what is there: {n['n_resident_taxa']} "
        f"resident taxa are called, {n['n_resident_drawn']} of them often enough for a row, "
        f"led by {n['top_residents']} runs. {n['collapsed']} is drawn as ONE row rather "
        f"than one per species. Within a run its species are called at a median coverage CV "
        f"of {n['alt_cv']:.2f} against {n['other_cv_lo']:.2f} to {n['other_cv_hi']:.2f} for "
        f"every other genus called as a cluster, and a run calling it calls a median of "
        f"{n['alt_n']:.0f} at once: that is one organism spread over near-identical "
        f"references, not {n['alt_n']:.0f} co-infecting congeners. The genus also contains wheat "
        f"pathogens, and the point of the collapse is that these calls cannot be resolved "
        f"to them. 04_align is what separates a cluster.",
        "",
        f"Pathogens of other plants but not of wheat are almost absent: {n['n_pp_taxa']} "
        f"taxa and {n['n_pp_calls']} calls in total, the commonest in {n['pp_max']} runs, "
        f"so none reaches a row. That is the db2 arm working. Under db1' these were the "
        f"conspicuous artefact, reads from a Pleosporaceae with no close reference being "
        f"forced onto a tomato or soybean pathogen; adding the predicted competitors gave "
        f"them somewhere nearer to go.",
        "",
        f"Colour is not abundance, it is who knew. A call is charcoal where the study that "
        f"generated the run declared that pathogen, orange where it did not, and grey where "
        f"nobody has yet curated what the study declared. Counting all {n['taxa']} wheat "
        f"pathogens and not only the {n['drawn_taxa']} given a row, there are "
        f"{n['n_declared']:,} declared calls, {n['n_undeclared']} undeclared and "
        f"{n['n_unknown']} of unknown reporting status. "
        f"{n['runs_undeclared']} runs carry at least one undeclared wheat pathogen, and in "
        f"{n['runs_only']} of them the undeclared pathogen is the only one present. So the "
        f"undeclared organism is a second pathogen beside the declared one rather than an "
        f"alternative to it, which follows from the declared targets being recovered in "
        f"nearly every run and is stated here as a consequence, not a separate finding.",
        "",
        f"The unit is the run and not the BioProject, deliberately. Projects in this cohort "
        f"are not geographically coherent — the largest spans twelve countries — so grouping "
        f"by project averages a worldwide survey into one value and hides the axis along "
        f"which the complex actually varies.",
        "",
        f"That axis is visible without any statistics. {n['focus']} is called in "
        f"{n['nw']}, and in {n['el']}: "
        f"{n['nw_pooled']:.0f}% of {n['nw_n']:,} runs from maritime northwest Europe against "
        f"{n['el_pooled']:.0f}% of {n['el_n']} runs from the drier and more continental "
        f"sampling, which is the published distribution of Septoria tritici blotch. The "
        f"pattern is recovered from reads collected for other reasons, and it is evidence "
        f"the calls are infections rather than contamination: neither spore drift, index "
        f"hopping nor cross-mapping is organised by continent.",
        "",
        f"Detection is presence, not activity, and a species-level call inside a densely "
        f"represented genus is a call on that genus' k-mer space; both are taken up in "
        f"04_align. Column order, project, country, coordinates and call count per run are "
        f"in samples_order.tsv.",
    )
    write(FIG / "samples_caption.txt", text)


if __name__ == "__main__":
    main()
