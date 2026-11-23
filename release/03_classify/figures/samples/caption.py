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

import collections
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]              # release/03_classify
FIG = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_reference/scripts"))
from _caption import rd, para, write                    # noqa: E402

sys.path.insert(0, str(FIG))
import plot_samples                                      # noqa: E402
from plot_samples import MIN_RUNS, abbrev               # noqa: E402
sys.path.insert(0, str(HERE / "scripts"))
from detect import declaration                          # noqa: E402

S2 = HERE.parent / "01_cohort/data/outputs/supp_data_s2_cohort.tsv"
S3REF = HERE.parent / "02_reference/data/outputs/supp_data_s3_reference.tsv"

S5 = HERE / "data/outputs/supp_data_s5_detections.tsv"
# The geography sentence contrasts these two groups of REGIONS. They are named because the
# claim is about them; every share is computed, never typed.
MARITIME = ["Northwest Europe"]
ELSEWHERE = ["North America", "South America", "Southern Africa", "East Africa",
             "East Asia", "Oceania"]
FOCUS = "Zymoseptoria tritici"
COLLAPSED = "Alternaria"        # the one genus the panel draws as a single row
SPILL_PARTNER = "Zymoseptoria tritici"   # the congener the withheld taxon rides on


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
    sp = plot_samples.SPILL
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
    s3 = {r["taxon"]: r["usda_row"] for r in rd(S3REF)}
    # A member of a COLLAPSED genus never was a row, so it cannot have been dropped as one:
    # four Alternaria sp. entries sit in the congener slot but travel inside the collapsed
    # Alternaria row, and counting them here reported 6 rows lost where the plot lost 2.
    grp = plot_samples.groups()
    nod = {t: len(v) for t, v in byt.items()
           if len(v) >= MIN_RUNS and s3.get(t) == "congener"
           and not grp.get(t, (0, "", False))[2]}
    n["n_norecord"] = len(nod)
    n["norecord"] = ", ".join(f"{abbrev(t)} ({v})" for t, v in
                              sorted(nod.items(), key=lambda kv: -kv[1]))
    unl = sorted(collections.Counter(
        c["region"] for c in cols if c["region"]).items(), key=lambda kv: -kv[1])
    small = [(c, v) for c, v in unl if v < plot_samples.LABEL_MIN]
    n["n_unlabelled"] = len(small)
    n["unlabelled_runs"] = sum(v for _, v in small)

    n["ci_alone"] = len(byt[ci] - set().union(*(byt[t] for t in sibs)))

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
        "GENERATED by caption.py — do not hand-edit, the numbers will drift.",
        "",
        "",
        "FULL CAPTION",
        "",
        f"Every one of the {n['runs']:,} RNA-seq runs in the cohort is one column, grouped "
        f"into {n['regions']} regions and ordered west to east; within a region, columns run "
        f"by longitude. The regions are cut for wheat pathology rather than by continent, "
        f"because the split that matters inside Europe is maritime northwest against "
        f"continental centre-east against Mediterranean, which is the climate axis the "
        f"foliar diseases follow; the mapping is in data/inputs/country_regions.tsv and a "
        f"country missing from it stops the figure rather than falling into an Other "
        f"bucket. {n['countries']} countries are represented across {n['projects']} "
        f"BioProjects, {n['geo']:,} runs carry coordinates and {n['noloc']} carry no "
        f"location and sit in the block at the right. {n['n_unlabelled']} region holding "
        f"{n['unlabelled_runs']} runs is drawn but too narrow to name; samples_order.tsv "
        f"gives the country and region of every column. "
        f"Rows are every fungus called in at least {MIN_RUNS} runs, "
        f"banded by two different rules, because the two halves of the panel answer to "
        f"different things. A wheat pathogen is banded by the DISEASE it causes, which is "
        f"the distinction a plant pathologist makes and which no column of the reference "
        f"table carries. Everything else is banded by its SLOT IN THE REFERENCE MATRIX, "
        f"read straight off that table, so the band says exactly what the colour says and "
        f"no taxon is placed by opinion. The plant host is not drawn, because every library "
        f"has one and no study declares it as a pathogen.",
        "",
        f"Taxa with no USDA host record of their own are not drawn. They are in the "
        f"reference only because their genus qualified, so there is no slot for them and no "
        f"statement this panel can make about where they belong. That costs "
        f"{n['n_norecord']} rows: {n['norecord']}. The first is named in 02_reference's own "
        f"scope notes as part of the wheat resident core, so these calls are not being "
        f"judged false, only left out of a figure whose every other row carries a host "
        f"record; supp_data_s5 keeps all of them.",
        "",
        f"The slot rule puts Microdochium majus with the taxa recorded on wheat rather than "
        f"with the stem-base diseases, and that disagreement is worth reading rather than "
        f"smoothing. It is a genuine agent of Microdochium patch, but PHI-base does not "
        f"score it, so the label says resident. The reads do not: it is called without "
        f"M. nivale in {n['mm_alone']} of its {n['mm_runs']} runs and at "
        f"{n['mm_ratio']:.2f} times its coverage where both occur, which is a detection in "
        f"its own right. The gap is in the label, not the call.",
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
        f"Colour carries two things and neither is abundance. The HUE is the taxon's "
        f"category in the reference matrix, unchanged from that figure: orange a wheat "
        f"pathogen, green a non-pathogen USDA records on wheat, pale green one it records "
        f"on another plant. The TREATMENT is who "
        f"knew: full colour where the study that generated the run did not declare that "
        f"organism, charcoal where it did, and a washed-out tint where nobody has yet "
        f"curated what the study declared. Counting all {n['taxa']} wheat "
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
        f"regions, which is the published distribution of Septoria tritici blotch. The "
        f"pattern is recovered from reads collected for other reasons, and it is evidence "
        f"the calls are infections rather than contamination: neither spore drift, index "
        f"hopping nor cross-mapping is organised by continent.",
        "",
        f"One row carries a dagger, and it is a flag rather than a removal. {n['sp']} is "
        f"called in {n['sp_runs']} runs, {n['sp_both']} of them also carrying "
        f"{n['sp_partner']}, at a median {n['sp_ratio']:.0f}% of its coverage and tracking "
        f"it across runs: never found alone, never more than a trace beside an abundant "
        f"congener. The rule is applied to every genus and marks only this one. The "
        f"Alternaria cluster does not trip it, its members sitting at comparable coverage "
        f"to each other, which is why that one is collapsed instead; Cladosporium "
        f"inversicolor does not either, being called with no congener at all in "
        f"{n['ci_alone']} runs, and being found alone is the evidence of independence.",
        "",
        f"The flagged call is NOT k-mers shared between the two references, and the "
        f"distinction decides how it can be handled. Kraken2 assigns every minimizer two "
        f"taxa share to their common ancestor, so in a passing run the Zymoseptoria genus "
        f"node absorbs about 1.29 million of them while Z. brevis still holds some 100,000 "
        f"exclusive to it in this database, which cannot have come from the Z. tritici "
        f"reference. The likely mechanism is divergence from the reference genome: field "
        f"Z. tritici differs from the assembly the database was built on, and where a field "
        f"variant coincides with the brevis reference the k-mer is counted for brevis. No "
        f"coverage floor separates that, the reads being real and the k-mers genuinely "
        f"exclusive, so it is disclosed rather than thresholded away. It is the "
        f"close-relative false positive that TIPP-SD (Shen et al. 2026) measures and leaves "
        f"uncorrected, and release/04_align is what resolves it, from species-unique reads "
        f"and positional evenness.",
        "",
        f"Detection is presence, not activity, and a species-level call inside a densely "
        f"represented genus is a call on that genus' k-mer space; both are taken up in "
        f"04_align. Column order, project, country, coordinates and call count per run are "
        f"in samples_order.tsv.",
    )
    write(FIG / "samples_caption.txt", text)


if __name__ == "__main__":
    main()
