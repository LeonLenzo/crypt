#!/usr/bin/env python3
"""resolve.py — turn the cohort into a fetch PLAN. Decides nothing about what is on disk.

Writes two tables:

    run_list.tsv    one row per cohort run: where it belongs, whether it is fetchable, and
                    whether it is in the pilot sample for its project
    fetch_plan.tsv  one row per BioProject: the numbers you read to decide what to launch

THE TWO INPUTS, and they must agree. Membership comes from `cohort_samples.tsv`: a run is in
the cohort if and only if it has a row there. But that table is keyed on RUN and carries no
BioSample column, so BioSample, LibraryStrategy, LibrarySelection and spots come from
`runs.tsv`. Any run present in one and not the other is a fault in 01_literature's release and
is raised, never defaulted: a missing BioSample would silently become its own sampling unit and
quietly change the pilot design.

## Two things this script decides, and they are different kinds of decision

**1. fetch_eligible — what the ASSAY can support.** An upstream drop, and legitimately so,
because it is about what a library can answer and not about what we expect to find in it:

    AMPLICON    single-locus PCR. 74 runs in PRJEB31334 at 0.48 M reads against 21 M for the
                shotgun runs in the same project. A CDS database scored against an rDNA
                amplicon produces confident nonsense, and all 55 cohort amplicon runs happen
                to be ones we have not downloaded, so this costs nothing and saves garbage.
    ncRNA-Seq   66 runs in PRJNA609211, size-selected small RNA, mean read 49.3 nt against
                234.1 nt for the RNA-Seq runs in the SAME project. Note the reason is NOT that
                49 nt is too short for k=35 -- it is not, that gives 15 k-mers. It is that the
                prep EXCLUDES the fungal mRNA this method detects, so a zero there is
                uninformative; and that k-mer coverage of a reference scales with read length,
                so a 49 nt and a 234 nt library have non-comparable detection floors for the
                same fungal abundance. Mixing them breaks the one-pipeline-for-every-sample
                rule that the whole cohort rests on.

**2. pilot — how much of a big project to FETCH FIRST.** NOT an exclusion, and the distinction
is the one `_exclusions.py` is built around ("exclude here, never upstream; nothing is skipped
because of what we expect to find"). Most of the large projects are physiology studies where
we have no reason to expect any fungal signal, and committing 3.2 TiB and a full classification
pass to one of them before knowing is the wrong order. So a seeded random sample of each large
project is fetched and classified first, and the decision to continue is made on measurements.

Three rules keep that honest:

  * **The pilot limits what we FETCH, never what we classify.** Everything on disk is
    classified. PRJEB39201 is 527/538 already downloaded, so its pilot flag is inert.
  * **The sample is RANDOM within the project, seeded, over BioSamples.** Random because the
    question is a proportion ("does this project contain fungal signal, and in what fraction of
    its samples"), and a random sample of runs is an unbiased estimate of that project's rate
    where a convenience sample is not. Over BioSamples, kept whole, because co-infection is
    called per BioSample and a half-sampled BioSample has a different detection sensitivity
    from a fully-sampled one.
  * **The denominator is the runs CLASSIFIED, and 05_detect must count REPORTS to get it.**
    `fetch_mode` here says what we intended to pull and is NOT that denominator: PRJEB39201
    trips the pilot threshold on run count but already has 527 of its 538 runs on disk, so it
    will be classified essentially whole. Coverage is a fact about which reports exist, not
    about what we planned to fetch, and a column in the fetch plan cannot know it. 100 runs
    with no detection bounds that project's rate below ~3% (95%, Clopper-Pearson one-sided on
    0/100); it does not license dropping the project from the pooled denominator. Doing that
    would make the pooled rate rise by construction, which is the failure mode the exclusion
    rule exists to prevent.

THREE PROJECTS CANNOT BE SAMPLED OVER BioSamples: PRJNA383416 registered 740 runs under ONE
BioSample, PRJNA721299 36 under one, PRJNA1181034 12 under one. Sampling BioSamples there draws
exactly one unit. Those fall back to sampling RUNS, recorded in `sample_unit`, and the same
registration defect means they cannot carry a per-BioSample co-infection claim at all -- which
is a 05_detect problem, flagged here because this is where it is discovered.

    python 02_classify/03_cohort/resolve.py
    python 02_classify/03_cohort/resolve.py --pilot 150
"""

import argparse
import csv
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import COHORT_SAMPLES, FETCH_PLAN, LIT_RUNS, RUN_LIST, SILVER  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _util import atomic_write  # noqa: E402

# Fixed so the sample is reproducible from the table alone. Changing it re-draws every pilot
# and orphans any reads already fetched under the old draw, so it is a constant, not an option.
SEED = 20261019

# A project gets a pilot if it is large in EITHER sense. Run count alone is not enough:
# PRJNA395908 is 112 runs but 634 GiB at 5.66 GiB/run, where PRJDB7234's whole 100-run pilot
# costs 25 GiB at 0.25 GiB/run. A 20x spread in per-run size means "how many runs" and "how
# much commitment" are different questions, and the pilot exists to bound the second.
PILOT_THRESHOLD = 150     # eligible runs
PILOT_GIB = 300           # or estimated volume, whichever trips first
PILOT_RUNS = 100          # 0/100 bounds a project's rate below ~3% at 95%

# What a CDS Kraken2 pass over shotgun RNA-seq can and cannot be scored on. See the module
# docstring: this is an assay judgment, not an expectation about content.
ELIGIBLE_STRATEGIES = {"RNA-Seq"}
INELIGIBLE_REASON = {
    "AMPLICON": "single-locus PCR, not shotgun: a CDS database cannot be scored on it",
    "ncRNA-Seq": "size-selected small RNA: excludes fungal mRNA, and read length makes the "
                 "detection floor non-comparable with the RNA-Seq runs in the same project",
}

# The strata that can be analysed as soon as reads exist, because the study itself went looking
# for disease. Under the two-evidence-standards rule a finding here is a SECOND pathogen in a
# known-infected sample; in `unselected` it is a FIRST one. Both count, and neither is dropped.
PATHOLOGY_SELECTION = {"disease-selected", "inoculated", "symptom-avoided"}

COLS = ["Run", "BioProject", "BioSample", "crop", "host", "setting", "sampling_selection",
        "tissue", "spots", "library_strategy", "library_selection",
        "fetch_eligible", "ineligible_reason", "pathology", "pilot", "sample_unit"]

PLAN_COLS = ["BioProject", "crops", "runs", "biosamples", "eligible", "ineligible",
             "pathology", "pilot_runs", "pilot_biosamples", "fetch_mode", "sample_unit",
             "est_gib", "pilot_gib", "note"]


def load() -> tuple[list[dict], dict]:
    """The cohort, joined to its library metadata. Raises on any disagreement."""
    cohort = list(csv.DictReader(open(COHORT_SAMPLES), delimiter="\t"))
    lit = {r["Run"]: r for r in csv.DictReader(open(LIT_RUNS), delimiter="\t")}
    missing = [r["Run"] for r in cohort if r["Run"] not in lit]
    if missing:
        sys.exit(f"{len(missing)} cohort runs absent from {LIT_RUNS.name}: {missing[:5]}\n"
                 "The two 01_literature tables disagree about which runs exist. Fix the "
                 "release, do not default the missing fields.")
    blank = [r["Run"] for r in cohort if not lit[r["Run"]].get("BioSample")]
    if blank:
        sys.exit(f"{len(blank)} cohort runs have no BioSample: {blank[:5]}\n"
                 "A blank BioSample becomes its own sampling unit and silently changes the "
                 "pilot design. Fix it upstream.")
    return cohort, lit


def pilot_set(runs_by_sample: dict, quota: int, rng: random.Random) -> set:
    """A seeded random sample of whole BioSamples, accumulated to a run quota.

    The quota is in RUNS and not in bytes, even though bytes are what the pilot is rationing.
    The decision the pilot informs is a PROPORTION ("in what fraction of this project's samples
    is there fungal signal"), and its confidence interval is a function of the number of samples
    scored, not of how many gigabytes they took. A byte quota would give a tight bound on
    shallow projects and a useless one on deep projects, which is backwards: the deep ones are
    exactly where the commitment is large enough to need the bound. Per-project pilot cost is
    reported in the plan instead, so the launch order can be chosen on it.

    BioSamples are kept WHOLE: once one is drawn, all of its runs come. So the quota is a
    target rather than an exact count, and the loop stops before overshooting rather than
    after, except that it always takes at least one unit.
    """
    units = sorted(runs_by_sample)        # sorted first, so the draw is seed-determined only
    rng.shuffle(units)
    picked, n = set(), 0
    for u in units:
        rs = runs_by_sample[u]
        if n and n + len(rs) > quota:
            continue                     # skip, do not stop: a later smaller unit may fit
        picked.update(rs)
        n += len(rs)
        if n >= quota:
            break
    return picked


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", type=int, default=PILOT_RUNS,
                    help=f"runs to fetch first from a large project (default {PILOT_RUNS})")
    ap.add_argument("--pilot-threshold", type=int, default=PILOT_THRESHOLD,
                    help=f"projects with more eligible runs than this get a pilot "
                         f"(default {PILOT_THRESHOLD})")
    ap.add_argument("--pilot-gib", type=float, default=PILOT_GIB,
                    help=f"or estimated GiB above which a project gets a pilot "
                         f"(default {PILOT_GIB})")
    ap.add_argument("--on-disk", type=Path,
                    help="file of run accessions already on scratch, one per line. Only "
                         "affects the printed cost summary; the tables stay "
                         "machine-independent by design")
    args = ap.parse_args()

    cohort, lit = load()

    by_proj = defaultdict(list)
    for r in cohort:
        by_proj[r["BioProject"]].append(r)

    # Mean bytes per spot, measured over the 3,225 runs already on scratch (median 171.8).
    # Used only to SIZE the plan; nothing decides anything on it.
    BYTES_PER_SPOT = 171.8

    rows, plan = [], []
    for proj in sorted(by_proj):
        prs = by_proj[proj]
        elig = [r for r in prs
                if lit[r["Run"]]["LibraryStrategy"] in ELIGIBLE_STRATEGIES]

        # Sampling unit: BioSample, unless the project's registration makes that meaningless.
        runs_by_sample = defaultdict(list)
        for r in elig:
            runs_by_sample[lit[r["Run"]]["BioSample"]].append(r["Run"])
        unit = "biosample"
        note = ""
        if len(elig) > 1 and len(runs_by_sample) == 1:
            unit = "run"
            note = (f"{len(elig)} runs registered under ONE BioSample: cannot be sampled or "
                    f"scored per BioSample")
            runs_by_sample = {r["Run"]: [r["Run"]] for r in elig}

        est = sum(int(lit[x["Run"]]["spots"] or 0) for x in elig) * BYTES_PER_SPOT
        needs_pilot = (len(elig) > args.pilot_threshold
                       or est / 2**30 > args.pilot_gib)
        if needs_pilot:
            rng = random.Random(f"{SEED}:{proj}")   # per project, so one project's size
            pilot = pilot_set(runs_by_sample, args.pilot, rng)   # cannot shift another's draw
        else:
            pilot = {r["Run"] for r in elig}        # fetch it whole; the flag is then inert

        for r in prs:
            l = lit[r["Run"]]
            strat = l["LibraryStrategy"]
            ok = strat in ELIGIBLE_STRATEGIES
            rows.append({
                "Run": r["Run"], "BioProject": proj, "BioSample": l["BioSample"],
                "crop": r["crop"], "host": r["host"], "setting": r["setting"],
                "sampling_selection": r["sampling_selection"], "tissue": r["tissue"],
                "spots": l["spots"], "library_strategy": strat,
                "library_selection": l["LibrarySelection"],
                "fetch_eligible": "yes" if ok else "no",
                "ineligible_reason": "" if ok else INELIGIBLE_REASON.get(
                    strat, f"library strategy {strat} is not shotgun RNA-seq"),
                "pathology": "yes" if r["sampling_selection"] in PATHOLOGY_SELECTION else "no",
                "pilot": "yes" if r["Run"] in pilot else "no",
                "sample_unit": unit,
            })

        pilot_bs = len({lit[x]["BioSample"] for x in pilot})
        pilot_est = sum(int(lit[x]["spots"] or 0) for x in pilot) * BYTES_PER_SPOT
        plan.append({
            "BioProject": proj,
            "crops": ",".join(sorted({r["crop"] for r in prs})),
            "runs": len(prs), "biosamples": len({lit[r["Run"]]["BioSample"] for r in prs}),
            "eligible": len(elig), "ineligible": len(prs) - len(elig),
            "pathology": sum(1 for r in prs
                             if r["sampling_selection"] in PATHOLOGY_SELECTION),
            "pilot_runs": len(pilot), "pilot_biosamples": pilot_bs,
            # What we INTENDED to pull. Not a denominator -- see the docstring: 05_detect
            # counts reports, because this cannot know what is already on disk.
            "fetch_mode": "pilot" if needs_pilot else "whole",
            "sample_unit": unit,
            "est_gib": f"{est / 2**30:.1f}",
            "pilot_gib": f"{pilot_est / 2**30:.1f}",
            "note": note,
        })

    SILVER.mkdir(parents=True, exist_ok=True)
    for path, cols, data in ((RUN_LIST, COLS, rows), (FETCH_PLAN, PLAN_COLS, plan)):
        with atomic_write(path) as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
        print(f"-> {path.name}  ({len(data)} rows)")

    # ── summary ─────────────────────────────────────────────────────────────────────────────
    elig = [r for r in rows if r["fetch_eligible"] == "yes"]
    pil = [r for r in elig if r["pilot"] == "yes"]
    print(f"\n{len(rows)} cohort runs over {len(plan)} BioProjects")
    print(f"  fetch-eligible          {len(elig):6d}")
    for strat, n in Counter(r["library_strategy"] for r in rows
                            if r["fetch_eligible"] == "no").most_common():
        print(f"    dropped, {strat:12s} {n:6d}")
    print(f"  first fetch (pilot)     {len(pil):6d}"
          f"   of which pathology-stratum {sum(1 for r in pil if r['pathology'] == 'yes')}")
    print(f"  pathology stratum       {sum(1 for r in rows if r['pathology'] == 'yes'):6d}")

    on_disk = set()
    if args.on_disk:
        on_disk = {ln.strip() for ln in open(args.on_disk) if ln.strip()}
    want = {r["Run"] for r in elig if r["pilot"] == "yes"}
    spots = {r["Run"]: int(r["spots"] or 0) for r in rows}

    # Ordered by what the FIRST fetch costs, because that is the launch order: the cheap
    # unrepresented crops come first on both run count and bytes.
    print(f"\n{'BioProject':14s} {'runs':>5s} {'elig':>5s} {'pilot':>6s} {'mode':>6s} "
          f"{'pilotGiB':>9s} {'toGetGiB':>9s}  crops")
    for p in sorted(plan, key=lambda x: float(x["pilot_gib"])):
        todo = sum(spots[r["Run"]] for r in rows
                   if r["BioProject"] == p["BioProject"] and r["Run"] in want
                   and r["Run"] not in on_disk) * BYTES_PER_SPOT / 2**30
        print(f"{p['BioProject']:14s} {p['runs']:5d} {p['eligible']:5d} {p['pilot_runs']:6d} "
              f"{p['fetch_mode']:>6s} {float(p['pilot_gib']):9.1f} {todo:9.1f}  "
              f"{p['crops'][:30]}")
    tot_pilot = sum(p["pilot_runs"] for p in plan)
    print(f"\n  every eligible run:              {sum(float(p['est_gib']) for p in plan) / 1024:6.2f} TiB"
          f"   ({len(elig)} runs)")
    print(f"  the first fetch, from zero:      {sum(float(p['pilot_gib']) for p in plan) / 1024:6.2f} TiB"
          f"   ({tot_pilot} runs)")
    if on_disk:
        todo_runs = [r for r in elig if r["pilot"] == "yes" and r["Run"] not in on_disk]
        print(f"  the first fetch, NET of scratch:  "
              f"{sum(spots[r['Run']] for r in todo_runs) * BYTES_PER_SPOT / 2**40:5.2f} TiB"
              f"   ({len(todo_runs)} runs to pull, "
              f"{tot_pilot - len(todo_runs)} already held)")
    flagged = [p for p in plan if p["note"]]
    if flagged:
        print(f"\n{len(flagged)} projects need a 05_detect caveat:")
        for p in flagged:
            print(f"  {p['BioProject']:14s} {p['note']}")


if __name__ == "__main__":
    main()
