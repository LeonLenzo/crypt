#!/usr/bin/env python3
"""prep_depth_funnel.py -- why each crop's runs leave the analysis, and the depth behind it.

The inclusion funnel drops rice entirely, and "rice is 0%" is not a result: it is a statement
about what the archive holds. This writes the two things needed to show that honestly, the
per-run sequencing depth and the reason each run is excluded.

## THE FUNNEL, REPRODUCED EXACTLY

Four nested filters, applied in this order, checked against the counts recorded on 2026-10-27:

    cohort as curated                 7,688 runs  41 projects
    disease_exposure == open          6,382       31
    + not inoculated                  6,339       29
    + setting == field                6,235       28
    + Gkmer >= 0.463                  3,680       27

**The counts moved +100 on 2026-10-27 and the check is what surfaced it**, which is the only
reason the targets are hardcoded. PRJEB84390's `disease_exposure` was corrected from `unknown`
to `open` (Leon's call): the 100 runs are field-collected, naturally Pt-infected wheat from 11
countries, so exposure occurred by construction. The earlier `unknown` applied rule 2 to a
project whose evidence was already in `provenance.tsv` under setting/tissue/selection and simply
never wired to this field. 100 wheat runs joined the retained set, all of them already fetched
and classified. Previous counts, for the record: 6,282 / 6,239 / 6,135 / 3,580.

`setting == "field"` is EXACT, not a prefix match. 104 runs are
`field (chamber-raised, transferred outdoors)`, which is not grown in the field throughout, and
a `startswith` admits them and misses the step by exactly that many.

## DEPTH IS Gkmer, NOT READS OR BASES

ENA's `read_len` is bases per SPOT, so a paired library's per-read length is half of it. At
k=35 a 51 nt read yields 17 k-mers and a 150 nt read 116, a factor of 7 on the same base count.
Rice is 1,636 single-end 51 nt runs out of 2,784, which is most of why it sits where it does.

Writes depth_funnel.tsv (one row per cohort run) and depth_funnel_summary.tsv.

    python 02_classify/report/figures/prep_depth_funnel.py
"""

import csv
import statistics as st
import textwrap
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import COHORT_SAMPLES, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

FLOOR = 0.463
CURATION = MODULE.parent / "01_literature/seeds/curation.tsv"
TSV = HERE / "depth_funnel.tsv"
SUMMARY = HERE / "depth_funnel_summary.tsv"
# The funnel as a TABLE, per crop, one row per gate. The printed funnel below checks the
# totals but is not an output, so the per-crop losses existed only in whoever ran it. This is
# the table the Methods needs: a reader cannot reconstruct the analysed set from prose, and
# every gate falls on the crops very unevenly (every rain-exclusion run is rice, every
# inoculated run is maize or wheat).
STEPS = HERE / "depth_funnel_steps.tsv"
# The caption for Supplementary Table S4, GENERATED from the counts so it cannot go stale.
# The funnel is a RESULT, not bookkeeping (Leon, 2026-10-27): it is the only place a reader can
# see that rice leaves the analysis entirely and why, and that each gate falls on one crop.
# Captions are generated here for the reason in `crypt-figure-caption-discipline`: a hand-edited
# caption beside a regenerated table is how a stale number survives a rebuild.
CAPTION = HERE / "depth_funnel_caption.txt"

# Grouped for the figure. The raw reason is kept in its own column so nothing is lost.
GROUP = {
    "retained": "retained",
    "below depth floor": "below the depth floor",
    "inoculated": "disease suppressed or manipulated by design",
    "not field throughout": "disease suppressed or manipulated by design",
    "exposure: protected": "disease suppressed or manipulated by design",
    "exposure: sheltered": "disease suppressed or manipulated by design",
    "exposure: unknown": "exposure undetermined",
}


def gkmer(d):
    """Gkmer for one run. read_len is bases per SPOT, so halve it for a paired library."""
    try:
        bases, L = float(d["bases"]), float(d["read_len"])
    except (TypeError, ValueError, KeyError):
        return None
    if d.get("layout") == "PAIRED":
        L = L / 2
    return bases * (L - 34) / L / 1e9 if L > 35 else None


def main() -> None:
    coh = list(csv.DictReader(open(COHORT_SAMPLES), delimiter="\t"))
    dep = {r["Run"]: r for r in csv.DictReader(open(SILVER / "run_depth.tsv"), delimiter="\t")}
    expo = {r["BioProject"]: r["value"]
            for r in csv.DictReader(open(CURATION), delimiter="\t")
            if r["field"] == "disease_exposure"}

    rows = []
    for r in coh:
        d = dep.get(r["Run"], {})
        g = gkmer(d)
        e = expo.get(r["BioProject"], "unknown")
        if e != "open":
            reason = f"exposure: {e}"
        elif r["sampling_selection"] == "inoculated":
            reason = "inoculated"
        elif r["setting"] != "field":
            reason = "not field throughout"
        elif g is None or g < FLOOR:
            reason = "below depth floor"
        else:
            reason = "retained"
        rows.append({
            "Run": r["Run"], "BioProject": r["BioProject"], "crop": r["crop"],
            "setting": r["setting"], "sampling_selection": r["sampling_selection"],
            "disease_exposure": e,
            "layout": d.get("layout", ""), "read_len": d.get("read_len", ""),
            "gkmer": f"{g:.5f}" if g is not None else "",
            "reason": reason, "fate": GROUP[reason],
        })

    with atomic_write(TSV) as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    by = defaultdict(list)
    for r in rows:
        by[r["crop"]].append(r)
    with atomic_write(SUMMARY) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["crop", "n_runs", "n_projects", "retained", "pct_below_floor",
                    "median_gkmer", "q25_gkmer", "q75_gkmer", "dominant_layout"])
        for crop, v in sorted(by.items(), key=lambda kv: -len(kv[1])):
            g = sorted(float(x["gkmer"]) for x in v if x["gkmer"])
            lay = Counter(f"{x['layout']} {x['read_len']}" for x in v).most_common(1)[0]
            w.writerow([crop, len(v), len({x["BioProject"] for x in v}),
                        sum(1 for x in v if x["reason"] == "retained"),
                        f"{100 * sum(1 for x in g if x < FLOOR) / len(g):.1f}",
                        f"{st.median(g):.3f}", f"{g[int(.25 * (len(g) - 1))]:.3f}",
                        f"{g[int(.75 * (len(g) - 1))]:.3f}", f"{lay[0]} (n={lay[1]})"])

    # ONE definition of the gates, used by both the table and the check below. They were
    # written twice before and the printed check is the only thing that would have caught a
    # drift between them, which is backwards.
    GATES = (
        ("cohort as curated", lambda r: True, 7688),
        ("exposure is open", lambda r: r["disease_exposure"] == "open", 6382),
        ("not inoculated",
         lambda r: r["disease_exposure"] == "open"
         and r["sampling_selection"] != "inoculated", 6339),
        ("field throughout",
         lambda r: r["disease_exposure"] == "open"
         and r["sampling_selection"] != "inoculated" and r["setting"] == "field", 6235),
        (f"depth >= {FLOOR} Gkmer", lambda r: r["reason"] == "retained", 3680),
    )
    crops = sorted({r["crop"] for r in rows},
                   key=lambda c: -sum(1 for r in rows if r["crop"] == c))
    with atomic_write(STEPS) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["step", "gate", "measure"] + crops + ["TOTAL", "n_projects"])
        prev = None
        for i, (lbl, keep, _t) in enumerate(GATES):
            surv = [r for r in rows if keep(r)]
            w.writerow([i, lbl, "surviving"]
                       + [sum(1 for r in surv if r["crop"] == c) for c in crops]
                       + [len(surv), len({r["BioProject"] for r in surv})])
            if prev is not None:
                lost = [r for r in prev if not keep(r)]
                w.writerow([i, lbl, "lost at this gate"]
                           + [sum(1 for r in lost if r["crop"] == c) for c in crops]
                           + [len(lost), ""])
            prev = surv
        # The raw reason, which is finer than `fate`: the grouped "disease suppressed or
        # manipulated by design" is FOUR different exclusions, and they do not fall alike.
        for reason in ("exposure: unknown", "exposure: protected", "exposure: sheltered",
                       "inoculated", "not field throughout", "below depth floor", "retained"):
            sub = [r for r in rows if r["reason"] == reason]
            if sub:
                w.writerow(["", reason, "final reason"]
                           + [sum(1 for r in sub if r["crop"] == c) for c in crops]
                           + [len(sub), len({r["BioProject"] for r in sub})])

    # ── the caption, every number computed ──────────────────────────────────────────────────
    def n(**kw):
        return sum(1 for r in rows if all(r[k] == v for k, v in kw.items()))

    surv = {lbl: sum(1 for r in rows if keep(r)) for lbl, keep, _ in GATES}
    per_crop = {c: n(crop=c, reason="retained") for c in crops}
    n_proj = len({r["BioProject"] for r in rows})
    n_proj_ret = len({r["BioProject"] for r in rows if r["reason"] == "retained"})
    rice_dep = sorted(float(r["gkmer"]) for r in rows
                      if r["crop"] == "Rice" and r["reason"] == "below depth floor" and r["gkmer"])
    rice_design = (n(crop="Rice", reason="exposure: protected")
                   + n(crop="Rice", reason="exposure: sheltered")
                   + n(crop="Rice", reason="not field throughout"))
    rice_unk = n(crop="Rice", reason="exposure: unknown")

    # Composed as sentences and wrapped at the end. Writing the template with its own line
    # breaks around the interpolations wrapped mid-sentence and read as a broken caption.
    body = [
        "Runs surviving each of the four inclusion gates, applied in the order shown, by crop. "
        "Gates are NESTED, not independent: every run is attributed to the FIRST gate it fails, "
        "so a run excluded for disease exposure is not counted again at the depth floor. The "
        f"curated cohort of {surv['cohort as curated']:,} runs across {n_proj} BioProjects "
        f"reduces to {surv[GATES[-1][0]]:,} runs across {n_proj_ret} BioProjects.",

        "Each gate falls on the crops unevenly, and that is the result the table carries. All "
        f"{n(reason='exposure: protected') + n(reason='exposure: sheltered')} "
        "fungicide-protected and rain-excluded runs are rice, as are all "
        f"{n(reason='not field throughout')} runs not grown in the field throughout; all "
        f"{n(reason='inoculated')} inoculated runs are maize or wheat. The depth floor is the "
        "only gate that touches four crops.",

        f"RICE LEAVES THE ANALYSIS ENTIRELY: {per_crop.get('Rice', 0)} of {n(crop='Rice'):,} "
        f"rice runs are retained. {rice_design} are excluded by a documented feature of the "
        "study design (fungicide application, rain exclusion, or not field-grown throughout), "
        f"{rice_unk} are recorded as disease exposure undetermined, and {len(rice_dep):,} fall "
        "below the depth floor.",

        "The depth exclusion does not depend on where the floor is set. Those runs have a "
        f"median {rice_dep[len(rice_dep)//2]:.3f} Gkmer and a maximum of {rice_dep[-1]:.3f} "
        f"against a floor of {FLOOR}, and none is within 2x of it, so any floor at or above "
        f"{rice_dep[-1] * 1.05:.2f} Gkmer excludes exactly the same runs. Rice is therefore "
        "reported as outside the analysable set rather than as a crop in which no co-infection "
        "was found.",

        "Depth is Gkmer, the number of k-mer start positions a library provides at k = 35, not "
        "reads or bases. ENA read_len is bases per spot, so a paired library's per-read length "
        "is half of it, and at k = 35 a 51 nt read yields 17 k-mers against 116 for a 150 nt "
        "read.",

        "Retained by crop: "
        + ", ".join(f"{c} {per_crop[c]:,}" for c in crops if per_crop[c]) + ".",
    ]
    with atomic_write(CAPTION) as fh:
        fh.write("Supplementary Table S4. The inclusion funnel, per crop.\n")
        fh.write("Caption for report/figures/depth_funnel_steps.tsv "
                 "and depth_funnel.{png,pdf}\n\n")
        fh.write("GENERATED by report/figures/prep_depth_funnel.py -- do not hand-edit.\n\n\n")
        fh.write("FULL CAPTION\n\n")
        fh.write("\n\n".join(textwrap.fill(b, 88) for b in body) + "\n")

    # The funnel, printed so a mismatch against the recorded counts is visible immediately.
    print("  funnel:")
    for lbl, keep, target in GATES:
        n = sum(1 for r in rows if keep(r))
        print(f"    {lbl:<22}{n:>6}   {'ok' if n == target else f'MISMATCH, expected {target}'}")
    print("\n  retained by crop:",
          dict(Counter(r["crop"] for r in rows if r["reason"] == "retained")))


main()
