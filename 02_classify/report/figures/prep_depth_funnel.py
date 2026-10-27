#!/usr/bin/env python3
"""prep_depth_funnel.py -- why each crop's runs leave the analysis, and the depth behind it.

The inclusion funnel drops rice entirely, and "rice is 0%" is not a result: it is a statement
about what the archive holds. This writes the two things needed to show that honestly, the
per-run sequencing depth and the reason each run is excluded.

## THE FUNNEL, REPRODUCED EXACTLY

Four nested filters, applied in this order, checked against the counts recorded on 2026-10-23:

    cohort as curated                 7,688 runs  41 projects
    disease_exposure == open          6,282       30
    + not inoculated                  6,239       28
    + setting == field                6,135       27
    + Gkmer >= 0.463                  3,580       26

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
        ("exposure is open", lambda r: r["disease_exposure"] == "open", 6282),
        ("not inoculated",
         lambda r: r["disease_exposure"] == "open"
         and r["sampling_selection"] != "inoculated", 6239),
        ("field throughout",
         lambda r: r["disease_exposure"] == "open"
         and r["sampling_selection"] != "inoculated" and r["setting"] == "field", 6135),
        (f"depth >= {FLOOR} Gkmer", lambda r: r["reason"] == "retained", 3580),
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

    # The funnel, printed so a mismatch against the recorded counts is visible immediately.
    print("  funnel:")
    for lbl, keep, target in GATES:
        n = sum(1 for r in rows if keep(r))
        print(f"    {lbl:<22}{n:>6}   {'ok' if n == target else f'MISMATCH, expected {target}'}")
    print("\n  retained by crop:",
          dict(Counter(r["crop"] for r in rows if r["reason"] == "retained")))


main()
