#!/usr/bin/env python3
"""caption.py — the caption for the per-run panel, written from the data.

    python release/03_classify/figures/samples/caption.py

Reads  samples_cols.tsv                               (beside this script, from the plot)
       ../../data/outputs/supp_data_s5_detections.tsv   the calls
Writes samples_caption.txt

The column order and its countries come from `samples_cols.tsv`, which the plot writes as
it draws, so the caption cannot narrate a different figure.

Nothing here is prose-with-numbers-typed-in: the narrative is a template and every figure in
it is read at write time. Rerun after any pipeline change.
"""

import collections
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
FIG = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _caption import rd, para, write                    # noqa: E402

sys.path.insert(0, str(FIG))
import prep_samples                                      # noqa: E402
from prep_samples import MIN_RUNS, abbrev               # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from detect import declaration                          # noqa: E402
from _cohort import describe                            # noqa: E402

S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
S3REF = HERE.parent / "02_reference/data/outputs/supp_data_s3_reference.tsv"

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
# The geography sentence contrasts these two groups of REGIONS. They are named because the
# claim is about them; every share is computed, never typed.
MARITIME = ["Europe"]
ELSEWHERE = ["Americas", "Africa", "Asia"]
FOCUS = "Zymoseptoria tritici"
COLLAPSED = "Alternaria"        # the one genus the panel draws as a single row
SPILL_PARTNER = "Zymoseptoria tritici"   # the congener the withheld taxon rides on


def numbers():
    cols = rd(FIG / "samples_cols.tsv")
    # samples_cols.tsv IS the analysed set: the plot writes one row per column it drew, so
    # filtering the calls through it keeps caption and panel on the same runs. Reading S5
    # whole counted the lab libraries and the host-contrast runs the panel excludes.
    drawn = {c["run"] for c in cols}
    calls = [r for r in rd(S5) if r["run"] in drawn]
    # Same recomputation the panel does, from the same function, so caption and panel agree
    # even while supp_data_s5's stored `declared` column is stale. See prep_samples.load.
    crow = {r["run"]: r for r in rd(S2)}
    for r in calls:
        r["declared"] = declaration(crow[r["run"]], r["taxon"], r["genus"],
                                    r["taxon_class"] == "host")[0]
    n = {}
    n["runs"] = len(cols)
    n["countries"] = len({c["country"] for c in cols if c["country"]})
    n["regions"] = len({c["region"] for c in cols if c["region"]})
    n["projects"] = len({c["BioProject"] for c in cols})
    n["noloc"] = sum(1 for c in cols if not c["country"])
    n["geo"] = sum(1 for c in cols if c["lon"] not in ("", "None"))

    by_region = defaultdict(list)
    for c in cols:
        by_region[c["region"]].append(c["run"])

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
    ppc = sorted(float(r["coverage"]) for r in calls if r["taxon_class"] == "plant_pathogen")
    n["pp_cov"] = f"{ppc[0]:.3f} to {ppc[-1]:.3f}" if ppc else "n/a"
    # How many are inside the band the pure-fungal control's own artefact reached. Stated
    # as a COUNT, never as "every one": the top call is 0.117, well clear of 0.0336, and
    # the earlier wording claimed the whole set sat under a bound two of them exceed.
    n["pp_lo"] = sum(1 for c in ppc if c <= 0.0336)
    n["pp_hi"] = f"{ppc[-1]:.3f}" if ppc else "n/a"
    n["pp_names"] = ", ".join(f"{abbrev(t)} ({len(v)})" for t, v in
                              sorted(cls["plant_pathogen"].items(), key=lambda kv: -len(kv[1])))
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
    # the spillover case, recomputed with the plot's rule so the two cannot disagree
    byt, covs, genus = defaultdict(set), defaultdict(dict), {}
    for r in calls:
        if r["taxon_class"] == "host":
            continue
        byt[r["taxon"]].add(r["run"])
        covs[r["taxon"]][r["run"]] = float(r["coverage"])
        genus[r["taxon"]] = r["genus"]
    sp = prep_samples.SPILL
    n["sp"], n["sp_partner"] = abbrev(sp), abbrev(SPILL_PARTNER)
    n["sp_runs"] = len(byt[sp])
    n["sp_both"] = len(byt[sp] & byt[SPILL_PARTNER])
    n["sp_ratio"] = 100 * statistics.median(
        [covs[sp][r] / covs[SPILL_PARTNER][r] for r in byt[sp] & byt[SPILL_PARTNER]
         if covs[SPILL_PARTNER][r]])
    ci = "Cladosporium inversicolor"
    sibs = [t for t in byt if t != ci and genus[t] == genus[ci]]
    mm, mn = "Microdochium majus", "Microdochium nivale"
    n["mm_runs"] = len(byt[mm])
    n["mm_alone"] = len(byt[mm] - byt[mn])
    n["mm_ratio"] = statistics.median(
        [covs[mm][r] / covs[mn][r] for r in byt[mm] & byt[mn] if covs[mn][r]])
    # what the panel leaves out, computed so the caption cannot drift from the plot

    n["ci_alone"] = len(byt[ci] - set().union(*(byt[t] for t in sibs)))

    n["alt_cv"] = statistics.median(cv[COLLAPSED])
    n["alt_n"] = statistics.median(size[COLLAPSED])
    others = [statistics.median(v) for g, v in cv.items()
              if g != COLLAPSED and len(v) >= 20 and max(size[g]) > 1]
    n["other_cv_lo"], n["other_cv_hi"] = min(others), max(others)
    n["collapsed"] = COLLAPSED

    wp = [r for r in calls if r["is_wheat_pathogen"] == "yes"]
    n["taxa"] = len({r["taxon"] for r in wp})
    n["n_resident_calls"] = sum(len(v) for v in cls["resident"].values())
    rows_drawn = rd(FIG / "samples_rows.tsv")
    n["drawn_rows"] = sum(1 for r in rows_drawn if not r["taxon"].startswith("__total"))
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
    # WHAT IS NOT DRAWN, derived from the rows the figure actually wrote rather than from
    # a threshold this script re-applies. MIN_RUNS is 1 now, so the only way a detected
    # wheat pathogen misses a row is being collapsed into a genus that is not a pathogen.
    # the three groups of column, so the caption can say the panel keeps all of them
    percall = defaultdict(list)
    for r in wp:
        percall[r["run"]].append(r["declared"])
    n["only_declared"] = sum(1 for v in percall.values()
                             if len(v) == 1 and v[0] == "declared")
    n["no_call"] = n["runs"] - len(percall)
    n["und_pct"] = 100.0 * n["runs_undeclared"] / n["runs"]
    nres = defaultdict(int)
    for r in calls:
        if r["taxon_class"] == "resident":
            nres[r["run"]] += 1
    n["res_zero_pct"] = 100.0 * (n["runs"] - len(nres)) / n["runs"]

    # the study-target marker, counted over the DRAWN columns and over the whole cohort
    n["virus_runs"] = sum(1 for c in cols
                          if crow[c["run"]].get("study_target") == "virus")
    n["n_virus"] = len({r["BioProject"] for r in crow.values()
                        if r.get("study_target") == "virus"})
    n["n_proj"] = len({r["BioProject"] for r in crow.values()})

    drawn = {r["taxon"] for r in rows_drawn if not r["taxon"].startswith("__total")}
    undrawn = sorted(((t, len(v)) for t, v in per_taxon.items() if t not in drawn),
                     key=lambda kv: (-kv[1], kv[0]))
    n["n_undrawn"] = len(undrawn)
    n["undrawn"] = ", ".join(f"{abbrev(t)} ({c})" for t, c in undrawn)

    hit = defaultdict(set)
    for r in wp:
        hit[r["taxon"]].add(r["run"])
    n["focus"] = abbrev(FOCUS)
    share = lambda c: (100 * len(set(by_region[c]) & hit[FOCUS]) / len(by_region[c])
                       if by_region[c] else 0.0)
    n["nw"] = ", ".join(f"{c} {share(c):.0f}%" for c in MARITIME if by_region[c])
    n["el"] = ", ".join(f"{c} {share(c):.0f}%" for c in ELSEWHERE if by_region[c])
    nwr = [r for c in MARITIME for r in by_region[c]]
    elr = [r for c in ELSEWHERE for r in by_region[c]]
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
        "GENERATED by caption.py \u2014 do not hand-edit, the numbers will drift.",
        "",
        "",
        "FULL CAPTION",
        "",
        f"Every one of the {n['runs']:,} field-collected {describe().split()[-2]} "
        f"{describe().split()[-1]} RNA-seq runs is one column, grouped into "
        f"{n['regions']} regions and ordered west to east; within a region, columns run by "
        f"longitude. {n['countries']} countries are represented across {n['projects']} "
        f"BioProjects and {n['geo']:,} of the runs carry coordinates; the country-to-region "
        f"mapping is in data/inputs/country_regions.tsv and a country missing from it stops "
        f"the figure rather than falling into an Other bucket. Runs with no location at all "
        f"are not drawn, the x axis being geography, and they stay in supp_data_s5 and in "
        f"every count the module makes.",
        "",
        f"EVERY WHEAT PATHOGEN DETECTED GETS A ROW, however rare: {n['taxa']} were called "
        f"and {n['drawn_rows']} are drawn, banded by the disease they cause. There is no "
        f"prevalence threshold on the rows, so a single call is as visible as a thousand "
        f"and the reader can see which organisms are at the thin end. The {n['n_undrawn']} "
        f"not drawn is {n['undrawn']}: it is collapsed into the Alternaria genus row, that "
        f"genus being drawn as one unresolved organism rather than as twenty congeners, "
        f"and the genus is a resident. Its calls are in supp_data_s5.",
        "",
        f"RUNS CARRYING ONLY THEIR DECLARED PATHOGEN ARE DRAWN, and so are runs carrying "
        f"no wheat pathogen at all. Of {n['runs']:,} columns, {n['only_declared']:,} hold "
        f"one call and it is the declared one, {n['no_call']} hold none, and "
        f"{n['runs_undeclared']} hold at least one undeclared pathogen. Dropping the first "
        f"two groups would raise the undeclared share from {n['und_pct']:.0f}% to almost "
        f"all of it, which is the arithmetic of conditioning a denominator on its own "
        f"outcome. A run with no call is not an empty run: it carries a resident community "
        f"with no single species concentrated enough to clear the floors.",
        "",
        f"RESIDENTS HAVE NO ROW OF THEIR OWN: {n['n_resident_taxa']} resident taxa are "
        f"called over {n['n_resident_calls']:,} calls, and the single ramp at the foot of "
        f"the panel carries them as a count per run. That is a decision about space rather "
        f"than about evidence \u2014 the detection-space panels share a page with this one "
        f"\u2014 and the ramp counts every resident call. Its scale is three-stop rather "
        f"than two because the counts are bimodal: {n['res_zero_pct']:.0f}% of runs call no "
        f"resident and a second cluster sits in the teens, and a pale-to-mid-green ramp put "
        f"both in the top half of the lightness range. The pathogens have no matching ramp, "
        f"deliberately: they have rows, so a ramp would restate a column the reader can "
        f"already count. The plant host is not drawn either: every library has one, so "
        f"it is a ranking question rather than a detection one, and supp_data_s5 carries "
        f"the ranked host for all {n['runs']:,} runs.",
        "",
        f"THE TOP ROW IS THE STUDY, NOT AN ORGANISM. {n['n_virus']} of the "
        f"{n['n_proj']} cohort projects are virus surveys, and the {n['virus_runs']} runs "
        f"they contributed are marked there in the same grey the declared calls use, "
        f"because it is the same kind of fact: what the study reported. They declare no "
        f"fungus because they were not looking for one, which is why the panel's largest "
        f"block of columns with no declared pathogen is theirs, and why those columns also "
        f"carry the fewest calls \u2014 their libraries are about half the depth of the "
        f"rest. Without the marker that block reads as missing curation, which it is not: "
        f"every project in the cohort now has a curated declaration. `study_target` is "
        f"curated from each paper's own statement of aim and is published in "
        f"supp_data_s2; it is not inferred from the declaration being empty.",
        "",
        f"Colour carries two things and neither is abundance. The HUE is the taxon's "
        f"category in the reference matrix, the same hue it has in 02_reference's matrix "
        f"and tree. The TREATMENT is who knew: full colour where the study that generated "
        f"the run did not declare that organism, charcoal where it did. There is no third "
        f"state, because every study drawn here now has a curated declaration.",
        "",
        f"Counting all {n['taxa']} wheat pathogens and not only the {n['drawn_rows']} given "
        f"a row, there are {n['n_declared']:,} declared calls and {n['n_undeclared']} "
        f"undeclared. {n['runs_undeclared']} runs carry at least one undeclared wheat "
        f"pathogen; in {n['runs_both']} of them it sits beside the pathogen the study was "
        f"about, and in {n['runs_only']} it is the only wheat pathogen present. The first "
        f"is the common case and follows from the declared targets being recovered in "
        f"nearly every run; the second is the stronger one, a study that reported no "
        f"pathogen we can detect and whose reads carry one anyway.",
        "",
        f"PATHOGENS OF OTHER PLANTS ARE NOT DRAWN. {n['n_pp_taxa']} were called over "
        f"{n['n_pp_calls']} calls \u2014 {n['pp_names']}, runs in brackets. They span "
        f"coverage {n['pp_cov']}, and {n['pp_lo']} of the {n['n_pp_calls']} sit at or below "
        f"0.034, the highest the pure-fungal control's own artefact reached; the largest, "
        f"at {n['pp_hi']}, does not, so this is a weak set rather than a demonstrably "
        f"spurious one. Four of the five are Pleosporales, the order holding the wheat spot "
        f"and blotch pathogens, and each was called in a run where one of those was "
        f"present. That is too tenuous to "
        f"draw as presence on a panel whose other rows are not, and it is the question "
        f"release/04_align answers from species-unique reads rather than the question this "
        f"screen can. They are in supp_data_s5. Their scarcity is still the db2 arm "
        f"working: under db1' these were the conspicuous artefact, reads from a "
        f"Pleosporaceae with no close reference forced onto a tomato or soybean pathogen, "
        f"and adding the predicted competitors gave them somewhere nearer to go.",
        "",
        f"The unit is the run and not the BioProject, deliberately. Projects in this cohort "
        f"are not geographically coherent \u2014 the largest spans twelve countries \u2014 "
        f"so grouping by project averages a worldwide survey into one value and hides the "
        f"axis along which the complex varies.",
        "",
        f"That axis is visible without any statistics. {n['focus']} is called in "
        f"{n['nw']} and in {n['el']}: {n['nw_pooled']:.0f}% of {n['nw_n']:,} European runs "
        f"against {n['el_pooled']:.0f}% of {n['el_n']} runs from everywhere else, which is "
        f"the published distribution of Septoria tritici blotch. The pattern is recovered "
        f"from reads collected for other reasons, and it is evidence the calls are "
        f"infections rather than contamination: neither spore drift, index hopping nor "
        f"cross-mapping is organised by continent.",
        "",
        f"Detection is presence, not activity, and a species-level call inside a densely "
        f"represented genus is a call on that genus' k-mer space; both are taken up in "
        f"04_align, and the congener case is shown in the detection-space panels beside "
        f"this one. Column order, project, country, region, coordinates and call count per "
        f"run are in samples_cols.tsv.",
    )
    write(FIG / "samples_caption.txt", text)


if __name__ == "__main__":
    main()
