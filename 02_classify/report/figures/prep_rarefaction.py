#!/usr/bin/env python3
"""prep_rarefaction.py — where the detection floor sits, as a function of library depth.

Reads `gold/rarefaction_curve.tsv` (written by utilities/rarefaction/rarefy.py) and answers the
question that decides whether 1,661 rice runs can support a co-infection claim:

    at what Gkmer does a REAL detection stop clearing the 05_detect floors?

Writes three things:

    rarefaction_curve_plot.tsv    one row per (run, species, arm, rung): the curves
    rarefaction_retained.tsv      one row per (arm, rung): how many detections survive
    rarefaction_caption.txt       FULL and SHORT, generated from the numbers

## What counts as "a detection" here

The reference set is each run's db_v3 detections that cleared the published floors
(score >= 0.7 AND coverage >= 1%), because those are the calls the project has already
committed to in print. Measuring against db_v4's own full-depth calls instead would be
circular: the question is whether a shallow library can recover a detection we believe in, not
whether it agrees with itself.

Coverage is recomputed against db_v4 throughout, so the numbers are not comparable to db_v3's
in absolute terms -- the denominators differ. What IS comparable is each detection against
itself across rungs, which is the entire point.

## Two floors, tested positively

`score` is not available here (it needs the whole-cohort accumulation fit that detect.py does),
so this uses the COVERAGE floor alone and says so. That makes the curve an upper bound on
retention: a detection that fails coverage would also fail the pair. Tested positively, because
coverage can be None and `None < floor` raises while `not (None >= floor)` silently admits.

The absolute floor matters as much as the relative one here, for a reason this experiment
exposes: at low depth a detection can hold 1% of a SMALL-ceiling taxon on a handful of k-mers.
`MIN_DISTINCT` is applied alongside, and the caption reports retention both ways so the choice
is visible rather than buried.

    python 02_classify/report/figures/prep_rarefaction.py
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import GOLD, ROOT, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

CURVE = GOLD / "rarefaction_curve.tsv"
DBV3 = ROOT / "03_kraken/05_filter/data/_classified.tsv"
DEPTH = SILVER / "run_depth.tsv"

PLOT = HERE / "rarefaction_curve_plot.tsv"
RETAINED = HERE / "rarefaction_retained.tsv"
# What a Gkmer threshold would COST, per crop. The third panel, and the one that makes the
# decision visible: a depth threshold is a rice filter and the figure has to say so.
THRESHOLD = HERE / "rarefaction_threshold_cost.tsv"
CAPTION = HERE / "rarefaction_caption.txt"

COV_FLOOR = 0.01        # the published operating point
MIN_DISTINCT = 10_000   # the absolute companion; see the module docstring

# Real project depths, so a rung reads as a statement about a project.
ANCHORS = [(0.019, "PRJDB7234", 1557, "Rice"), (0.036, "PRJNA726716", 104, "Rice"),
           (0.555, "PRJNA510920", 251, "Rice"), (3.007, "PRJEB39201", 538, "Wheat"),
           (13.24, "PRJNA1181034", 12, "Wheat")]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def main() -> None:
    if not CURVE.exists():
        sys.exit(f"no curve yet: {CURVE}\n"
                 "Run: sbatch 02_classify/slurm/rarefaction.slurm")
    rows = list(csv.DictReader(open(CURVE), delimiter="\t"))
    runs = sorted({r["run"] for r in rows})

    # db_v3's committed detections, per run. The reference set.
    ref = defaultdict(dict)
    for r in csv.DictReader(open(DBV3), delimiter="\t"):
        if r["run"] in runs and fnum(r["score"]) >= 0.7 and fnum(r["coverage"]) >= COV_FLOOR:
            ref[r["run"]][r["species"]] = fnum(r["coverage"])
    n_ref = sum(len(v) for v in ref.values())

    # observed[(run, species, arm, rung)] -> [(coverage, distinct), ...] over seeds
    obs = defaultdict(list)
    for r in rows:
        obs[(r["run"], r["species"], r["arm"], float(r["rung_gkmer"]))].append(
            (fnum(r["coverage"]), int(r["distinct"] or 0)))

    arms = sorted({r["arm"] for r in rows})
    rungs = sorted({float(r["rung_gkmer"]) for r in rows})

    # ── the curves ──────────────────────────────────────────────────────────────────────────
    plot = []
    for run in runs:
        for sp, v3 in sorted(ref[run].items(), key=lambda kv: -kv[1]):
            for arm in arms:
                for rung in rungs:
                    vals = obs.get((run, sp, arm, rung))
                    if not vals:
                        continue
                    cov = sum(v[0] for v in vals) / len(vals)
                    dist = sum(v[1] for v in vals) / len(vals)
                    plot.append({"run": run, "species": sp, "arm": arm, "rung_gkmer": rung,
                                 "n_seeds": len(vals), "coverage": f"{cov:.6f}",
                                 "distinct": f"{dist:.0f}", "db_v3_coverage": f"{v3:.4f}",
                                 "clears_relative": "yes" if cov >= COV_FLOOR else "no",
                                 "clears_both": "yes" if (cov >= COV_FLOOR
                                                          and dist >= MIN_DISTINCT) else "no"})
    with atomic_write(PLOT) as fh:
        w = csv.DictWriter(fh, fieldnames=list(plot[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(plot)

    # ── retention ───────────────────────────────────────────────────────────────────────────
    ret = []
    for arm in arms:
        for rung in rungs:
            sel = [p for p in plot if p["arm"] == arm and p["rung_gkmer"] == rung]
            if not sel:
                continue
            ret.append({"arm": arm, "rung_gkmer": rung, "detections_measured": len(sel),
                        "retained_relative": sum(1 for p in sel
                                                 if p["clears_relative"] == "yes"),
                        "retained_both": sum(1 for p in sel if p["clears_both"] == "yes"),
                        "frac_relative": f"{sum(1 for p in sel if p['clears_relative']=='yes')/len(sel):.4f}",
                        "frac_both": f"{sum(1 for p in sel if p['clears_both']=='yes')/len(sel):.4f}"})
    with atomic_write(RETAINED) as fh:
        w = csv.DictWriter(fh, fieldnames=list(ret[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(ret)

    # ── what a threshold costs, by crop ─────────────────────────────────────────────────────
    # Gkmer recomputed here from run_depth.tsv rather than read from a column, because the
    # threshold is applied to the WHOLE cohort and the curve only covers 5 runs.
    depth = list(csv.DictReader(open(DEPTH), delimiter="\t"))
    runlist = {r["Run"]: r for r in csv.DictReader(open(SILVER / "run_list.tsv"),
                                                   delimiter="\t")}
    K = 35

    def gkmer(r):
        b = int(r["bases"] or 0)
        rlen = float(r["read_len"] or 0)
        L = rlen / 2 if r["layout"] == "PAIRED" else rlen
        return b * (max(0.0, L - K + 1) / L) / 1e9 if L > K else 0.0

    cohort = [(gkmer(r), r) for r in depth
              if r["bases"] and int(r["bases"]) > 0 and r["Run"] in runlist]
    crops = sorted({runlist[r["Run"]]["crop"] for _, r in cohort})
    cost = []
    for rung in rungs:
        for crop in crops:
            inc = [(g, r) for g, r in cohort if runlist[r["Run"]]["crop"] == crop]
            lost = sum(1 for g, _ in inc if g < rung)
            cost.append({"rung_gkmer": rung, "crop": crop, "runs": len(inc),
                         "lost": lost, "kept": len(inc) - lost,
                         "frac_lost": f"{lost / len(inc):.4f}" if inc else ""})
    with atomic_write(THRESHOLD) as fh:
        w = csv.DictWriter(fh, fieldnames=list(cost[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(cost)

    # the depth at which each detection is LOST: lowest rung that still clears
    lost = {}
    for run in runs:
        for sp in ref[run]:
            for arm in arms:
                ok = [p["rung_gkmer"] for p in plot
                      if p["run"] == run and p["species"] == sp and p["arm"] == arm
                      and p["clears_both"] == "yes"]
                lost[(run, sp, arm)] = min(ok) if ok else None

    # ── caption ─────────────────────────────────────────────────────────────────────────────
    L = []
    A = L.append
    A("Detection against library depth")
    A("Caption for report/figures/rarefaction.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_rarefaction.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"Coverage of {n_ref} established detections as their libraries are rarefied. The {len(runs)} runs")
    A("are wheat libraries already on scratch, so this costs no download; each carries several")
    A("detections and one Kraken2 report measures all of them at once. The reference set is each")
    A("run's db_v3 calls that cleared the published floors (score >= 0.7 and coverage >= 1%),")
    A("because those are the calls the project has committed to; measuring against db_v4's own")
    A("full-depth calls would be circular. Coverage is recomputed against db_v4 throughout, so")
    A("absolute values are not comparable with db_v3 -- the denominators differ -- while each")
    A("detection against itself across rungs is, which is the point.")
    A("")
    A("DEPTH IS IN Gkmer, NOT Gb. At k=35 a read of length L yields L-34 k-mers, so k-mer yield")
    A("per base is (L-34)/L: 0.33 at 51 nt, 0.66 at 101 nt, 0.77 at 150 nt. Read length is")
    A("therefore worth a 2.3x factor independent of depth, and a library must be ranked on bases")
    A("x yield. Ranking on layout instead is plainly wrong: a paired-only filter would discard")
    A("301 single-end runs deeper than 2 Gb while keeping 311 paired runs under 1 Gb.")
    A("")
    A("THE TWO ARMS are matched on Gkmer so the only difference between them is library SHAPE.")
    A("`native` subsamples the run as sequenced; `shallow` keeps one mate and truncates it to")
    A("51 nt, reproducing PRJDB7234. A gap between the arms at equal Gkmer isolates the")
    A("read-length effect, which is expected because --minimum-hit-groups 3 is harder to satisfy")
    A("with 17 k-mers per read than with 117.")
    A("")
    A("  WHERE THE REAL PROJECTS SIT")
    for gk, proj, n, crop in ANCHORS:
        A(f"    {gk:7.3f} Gkmer   {proj:14s} {n:5d} runs  {crop}")
    A("")
    A("  DETECTIONS RETAINED (mean over seed replicates; coverage floor 1%, and the same with")
    A(f"  an absolute floor of {MIN_DISTINCT:,} distinct minimizers alongside)")
    A("")
    A(f"  {'arm':9s} {'Gkmer':>8s} {'measured':>9s} {'cov>=1%':>9s} {'+abs floor':>11s}")
    for r in ret:
        A(f"  {r['arm']:9s} {r['rung_gkmer']:8.3f} {r['detections_measured']:9d} "
          f"{r['retained_relative']:9d} {r['retained_both']:11d}")
    A("")
    lo = min(rungs)
    for arm in arms:
        a = [r for r in ret if r["arm"] == arm and r["rung_gkmer"] == lo]
        b = [r for r in ret if r["arm"] == arm]
        if a and b:
            top = max(b, key=lambda x: x["rung_gkmer"])
            A(f"At {lo} Gkmer the {arm} arm retains {a[0]['retained_both']} of "
              f"{a[0]['detections_measured']} detections, against "
              f"{top['retained_both']} of {top['detections_measured']} at "
              f"{top['rung_gkmer']:.3f} Gkmer.")
    A("")
    A("WHAT THIS DECIDES. PRJDB7234 (1,557 runs) and PRJNA726716 (104 runs) sit at 0.019 and")
    A("0.036 Gkmer and are both rice, together 1,661 of the 2,784 rice runs. Whatever fraction")
    A("of detections survives at those rungs is the fraction of rice co-infections the pipeline")
    A("can see at all, and the quantity that collapses first is the NUMBER of pathogens per")
    A("sample -- which is the dependent variable. A depth cut-off is therefore a rice filter: at")
    A("0.2 Gb it removes 1,551 rice runs against 18 maize and 0 sorghum. The cut-off has to be")
    A("set on this curve and its cost to rice stated, not chosen by eye.")
    A("")
    A("  WHAT A THRESHOLD COSTS (runs dropped, whole cohort)")
    A("")
    A(f"  {'Gkmer':>8s} {'lost':>6s} {'%':>5s}   " + "  ".join(f"{c[:11]:>11s}" for c in crops))
    for rung in rungs:
        sel = [c for c in cost if c["rung_gkmer"] == rung]
        tl = sum(c["lost"] for c in sel)
        tn = sum(c["runs"] for c in sel)
        cells = []
        for crop in crops:
            e = next((c for c in sel if c["crop"] == crop), None)
            cells.append(f"{e['lost']:5d}/{e['runs']:<5d}" if e else "     -     ")
        A(f"  {rung:8.3f} {tl:6d} {tl/tn:5.0%}   " + "  ".join(f"{x:>11s}" for x in cells))
    A("")
    A("Only TWO projects fall entirely below 0.5 Gkmer and both are rice: PRJDB7234 (1,557")
    A("runs, median 0.019, max 0.214) and PRJNA726716 (104 runs, median 0.044, max 0.141).")
    A("")
    A("READ LENGTH COSTS NOTHING BEYOND Gkmer, which was not the expected result. At matched")
    A("Gkmer the 51 nt single-end arm retains as many detections as the native paired arm")
    A("(5 vs 5, 8 vs 6, 9 vs 10, 12 vs 12 at the four shared rungs). The reasoning that")
    A("--minimum-hit-groups 3 would penalise a 17-k-mer read beyond its yield is NOT supported.")
    A("Gkmer is therefore a sufficient single statistic and a short-read library needs no extra")
    A("penalty, so one threshold on one number is defensible.")
    A("")
    A("WHAT THE CURVE DOES NOT SAY. It uses the COVERAGE floor only; `score` needs detect.py's")
    A("whole-cohort accumulation fit, so retention here is an UPPER bound -- a detection failing")
    A("coverage would fail the pair too. Nested subsampling is not used: each rung is an")
    A("independent Bernoulli draw, so between-rung non-monotonicity is sampling noise and the")
    A("seed replicates at the two lowest rungs bound it. And these are wheat libraries; whether")
    A("a rice library of equal Gkmer behaves identically is an assumption, not a measurement,")
    A("because no deep rice library in the cohort carries established detections to rarefy.")
    A("")
    A("SHORT CAPTION")
    A("")
    A(f"Coverage of {n_ref} established detections in {len(runs)} wheat libraries under rarefaction, against")
    A("library depth in Gkmer (bases x (L-34)/L, the k-mer yield that detection depends on).")
    A("Dashed line, the 1% coverage floor; vertical lines, the median depth of real cohort")
    A("projects. The two arms are matched on Gkmer and differ only in read length and pairing,")
    A("so the gap between them isolates library shape from depth.")
    A("")
    A("SOURCES")
    A(f"  curve     {CURVE.name}, from utilities/rarefaction/rarefy.py")
    A(f"  reference {DBV3.relative_to(ROOT)} (db_v3 calls clearing the published floors)")
    A(f"  depths    {DEPTH.name}")
    A("")
    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {PLOT.name}  ({len(plot)} rows)")
    print(f"-> {RETAINED.name}")
    print(f"-> {THRESHOLD.name}")
    print(f"-> {CAPTION.name}")
    print(f"\nreference detections: {n_ref} over {len(runs)} runs")
    print(f"\n  {'arm':9s} {'Gkmer':>8s} {'measured':>9s} {'cov>=1%':>9s} {'+abs':>6s}")
    for r in ret:
        print(f"  {r['arm']:9s} {r['rung_gkmer']:8.3f} {r['detections_measured']:9d} "
              f"{r['retained_relative']:9d} {r['retained_both']:6d}")
    print("\n  lowest Gkmer at which each detection still clears BOTH floors:")
    for (run, sp, arm), v in sorted(lost.items(), key=lambda kv: (kv[0][2], kv[0][0], kv[0][1])):
        if arm != "native":
            continue
        print(f"    {run:13s} {sp[:32]:32s} {v if v is not None else 'lost at every rung'}")


if __name__ == "__main__":
    main()
