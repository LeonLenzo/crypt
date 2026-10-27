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
# THE FOUR CROPS THE PAPER COVERS. Leon, 2026-10-27: the supplementary table has to be true to
# the paper, which claims wheat, maize, rice and sorghum, so "other cereal" and "barley" are not
# columns in it. They are not a rounding error either: dropping them removes 58 runs, 56 of them
# previously retained, and two whole projects (PRJNA721299, 36 Setaria italica KERNEL runs, and
# PRJNA1232637, 1 Hordeum vulgare run). The retained "barley" was 2 runs of Hordeum jubatum, a
# wild grass rather than the crop. Cohort 7,688 -> 7,630 and retained 3,680 -> 3,624.
PAPER_CROPS = ("Wheat", "Maize", "Rice", "Sorghum")
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

    # ── THE PUBLICATION TABLE: three stages, two filters of different KIND ────────────────
    # Leon's design, 2026-10-27. The two filters are not the same sort of thing and the table
    # says so: "Naturally Infected" is METHODOLOGICAL (what the originating study did to the
    # plants, or failed to report) and "Library Depth" is TECHNICAL (what the sequencing can
    # support). Collapsing them into one list of five reasons hid that distinction, which is
    # the one a reader needs in order to judge whether an exclusion is about the world or
    # about us.
    #   Filter rows carry the LOSS, not the survivors, with a Total Remaining row at the foot.
    # Restricted to PAPER_CROPS.
    pc = [r for r in rows if r["crop"] in PAPER_CROPS]

    def methodological(r):
        return (r["disease_exposure"] == "open"
                and r["sampling_selection"] != "inoculated"
                and r["setting"] == "field")

    def technical(r):
        return bool(r["gkmer"]) and float(r["gkmer"]) >= FLOOR

    s1 = pc
    s2 = [r for r in s1 if methodological(r)]
    s3 = [r for r in s2 if technical(r)]

    def cells(sub):
        return [sum(1 for r in sub if r["crop"] == c) for c in PAPER_CROPS]

    def lost(a, b):
        keep = {id(r) for r in b}
        return [r for r in a if id(r) not in keep]

    l2, l3 = lost(s1, s2), lost(s2, s3)
    with atomic_write(STEPS) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([""] + list(PAPER_CROPS) + ["TOTAL"])
        w.writerow(["Literature Search"] + cells(s1) + [len(s1)])
        w.writerow(["Naturally Infected"] + [-x if x else 0 for x in cells(l2)] + [-len(l2)])
        w.writerow(["Library Depth"] + [-x if x else 0 for x in cells(l3)] + [-len(l3)])
        w.writerow(["Total Remaining"] + cells(s3) + [len(s3)])
        # Not part of the manuscript table, but the caption quotes these and they must be
        # checkable against the same run-level file rather than retyped.
        w.writerow([])
        w.writerow(["# sub-reasons within the methodological filter"]
                   + list(PAPER_CROPS) + ["TOTAL", "n_projects"])
        for reason in ("exposure: unknown", "exposure: protected", "exposure: sheltered",
                       "inoculated", "not field throughout"):
            sub = [r for r in pc if r["reason"] == reason]
            if sub:
                w.writerow([reason] + cells(sub)
                           + [len(sub), len({r["BioProject"] for r in sub})])

    # ── the caption: what the table IS, not what it shows ───────────────────────────────────
    # Leon, 2026-10-27: "the caption doesn't need to describe the data, just what it is". So it
    # defines the three stages, the two kinds of filter and the depth unit, and stops. The
    # findings the earlier draft narrated (rice leaving entirely, each filter falling on one
    # crop) are RESULTS and belong in the Results text, not here. The sub-reasons are listed
    # because they define what the methodological filter is; their counts are data and are in
    # the table's second block.
    body = [
        "Runs removed at each stage of the inclusion funnel, for the four cereal crops this "
        "study covers. Literature Search is the curated cohort. Naturally Infected is a "
        # lowercase in the caption: the distinction is carried by the sentence, and caps read as
        # shouting on a printed page. The emphasis belongs in the code comment, not the paper.
        "methodological filter, removing runs whose originating study did not provide an "
        "unmanipulated, naturally infected field sample: disease exposure not determinable from "
        "the paper, fungicide protection, rain exclusion, inoculation rather than natural "
        "infection, or plants not field-grown throughout. Library Depth is a technical filter, "
        f"removing libraries below {FLOOR} Gkmer, the count of k-mer start positions a library "
        "provides at k = 35. Filter rows give runs removed; Total Remaining is the analysed set.",
    ]
    with atomic_write(CAPTION) as fh:
        fh.write("Supplementary Table T4. The inclusion funnel.\n")
        fh.write("Caption for report/figures/depth_funnel_steps.tsv\n\n")
        fh.write("GENERATED by report/figures/prep_depth_funnel.py -- do not hand-edit.\n\n\n")
        fh.write("FULL CAPTION\n\n")
        fh.write("\n\n".join(textwrap.fill(b, 88) for b in body) + "\n")

    # The self-check still runs over the WHOLE cohort, not PAPER_CROPS, because its targets are
    # the recorded cohort funnel and the point is to catch curation drift.
    GATES = (
        ("cohort as curated", lambda r: True, 7688),
        ("exposure is open", lambda r: r["disease_exposure"] == "open", 6382),
        ("+ not inoculated",
         lambda r: r["disease_exposure"] == "open"
         and r["sampling_selection"] != "inoculated", 6339),
        ("+ field throughout",
         lambda r: r["disease_exposure"] == "open"
         and r["sampling_selection"] != "inoculated" and r["setting"] == "field", 6235),
        ("+ depth floor", lambda r: r["reason"] == "retained", 3680),
    )

    # The funnel, printed so a mismatch against the recorded counts is visible immediately.
    print("  funnel:")
    for lbl, keep, target in GATES:
        n = sum(1 for r in rows if keep(r))
        print(f"    {lbl:<22}{n:>6}   {'ok' if n == target else f'MISMATCH, expected {target}'}")
    print("\n  retained by crop:",
          dict(Counter(r["crop"] for r in rows if r["reason"] == "retained")))


main()
