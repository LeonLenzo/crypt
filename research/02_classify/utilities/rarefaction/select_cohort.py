#!/usr/bin/env python3
"""select_cohort.py -- choose the runs for the cohort-wide rarefaction, and say why.

    python 02_classify/utilities/rarefaction/select_cohort.py            # the cohort arm
    python 02_classify/utilities/rarefaction/select_cohort.py --controls # the control arm
    python 02_classify/utilities/rarefaction/select_cohort.py --ceiling 3.0

Writes a comma-separated run list to stdout, for `rarefy.py --select --runs "$(...)"`. It does
not touch rarefy.py's measurement, which is unchanged and already array-shaped.

WHY THIS REPLACES FIVE HAND-PICKED LIBRARIES. The depth floor removes 2,321 runs, more than
every study-design criterion put together, and until now it rested on a curve measured in five
WHEAT libraries carrying 22 detections. Whether the same curve holds for maize or sorghum was
an assumption the caption had to concede. There is no longer any reason for that: 2,908
analysed runs clear 2 Gkmer, near-evenly split across the three crops, and classifying them
costs roughly 1,300 core-hours against 256,000 SU remaining. Leon, 2026-10-30.

THE CONSTRAINT IS STATISTICAL, NOT COMPUTATIONAL, and it decides the design. A library can
only be subsampled DOWNWARD. Rarefy everything and the 3 Gkmer rung is measured only on runs
already deeper than 3 Gkmer, so the denominator changes at every rung -- and changes
non-randomly, because deep runs skew to disease-selected wheat. The curve would then confound
depth with which libraries happen to be deep.
  So the set is fixed at the CEILING: every run that can be carried through every rung, and no
others. At 2.0 Gkmer that is 2,908 runs over 21 projects, wheat 1,012 / maize 974 /
sorghum 922. Raising the ceiling to 3.0 adds a rung but collapses wheat to 616 and the balance
with it; lowering it to 1.0 keeps more runs but cannot probe above the floor it is meant to
justify.

TWO ARMS, AND THE CONTROL ARM IS THE POINT. The cohort arm asks: at what depth do we stop
recovering the calls we would have made at full depth? That is self-referential -- it cannot
speak to what full depth itself misses -- and the original caption flagged exactly that.
The control arm asks the same question of libraries with KNOWN ground truth: a single
inoculated pathogen plus its host, from PRJEB29046, PRJEB58154 and PRJNA613154. Losing a
detection there is losing one known to be real. Run both on the same rung ladder and the
control arm anchors what the cohort arm can only assume.
"""

import argparse
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import SILVER  # noqa: E402

TABLES = MODULE.parent / "01_literature" / "report" / "tables"
REPORTS = MODULE / "data" / "bronze" / "reports"
COINFECTION = MODULE / "data" / "gold" / "coinfection.tsv"
S2 = TABLES / "Supplementary Data S2.tsv"
S3 = TABLES / "Supplementary Data S3.tsv"
CROPS = ("Wheat", "Maize", "Rice", "Sorghum")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ceiling", type=float, default=2.0,
                    help="Gkmer; every selected run must reach this, so the denominator is "
                         "fixed across rungs (default 2.0)")
    ap.add_argument("--controls", action="store_true",
                    help="the control arm instead: known single-pathogen libraries")
    ap.add_argument("--null-sample", type=int, default=300,
                    help="how many ZERO-detection runs to carry for the spurious-gain "
                         "measurement, sampled deterministically per crop (default 300; "
                         "0 to measure recovery only)")
    ap.add_argument("--all-zero", action="store_true",
                    help="keep every zero-detection run instead of sampling. Doubles the "
                         "compute to measure a rate that a few hundred runs already pin down")
    ap.add_argument("--no-require-report", action="store_true",
                    help="skip the check that every selected run is already classified. "
                         "Only for inspecting a selection; the measurement needs the reports")
    ap.add_argument("--max-per-project", type=int, default=0,
                    help="cap per BioProject, 0 for no cap. A cap trades the fixed "
                         "denominator for balance and should be a deliberate choice")
    args = ap.parse_args()

    depth = {r["Run"]: r for r in
             csv.DictReader(open(SILVER / "run_depth.tsv"), delimiter="\t")}

    if args.controls:
        cd = SILVER / "control_depth.tsv"
        if not cd.exists():
            sys.exit("no control_depth.tsv; run utilities/controls/control_depth.py")
        rows = [(r["Run"], r["BioProject"], "control") for r in
                csv.DictReader(open(cd), delimiter="\t")]
        # Controls are not filtered on the ceiling: there are only 26 and every one is worth
        # measuring. rarefy.py skips any rung above a run's own depth, so a thin control
        # simply contributes to fewer rungs rather than being excluded.
        print(f"  control arm: {len(rows)} runs over "
              f"{len({b for _, b, _ in rows})} accessions", file=sys.stderr)
    else:
        if not (S2.exists() and S3.exists()):
            sys.exit("missing S2/S3; run report/tables/build_s2.py and build_s3.py")
        s2 = {r["Run"]: r for r in csv.DictReader(open(S2), delimiter="\t")}
        s3 = {r["Run"]: r for r in csv.DictReader(open(S3), delimiter="\t")}
        rows = []
        for run, r in s2.items():
            t = s3.get(run)
            if not t or t.get("analysed") != "1":
                continue
            if r["crop"] not in CROPS:
                continue
            g = depth.get(run, {}).get("bases")
            try:
                gk = float(t["gkmer"]) if t["gkmer"] else 0.0
            except ValueError:
                gk = 0.0
            if gk >= args.ceiling and g:
                rows.append((run, r["BioProject"], r["crop"]))

    # EVERY SELECTED RUN MUST ALREADY BE CLASSIFIED, and this is checked rather than assumed.
    # The cohort arm's reference is each run's own full-depth calls, so a run with no report
    # has nothing to be compared against and would silently contribute zero recovered
    # detections at every rung -- indistinguishable from a run that genuinely lost them, and
    # pulling the curve down in proportion to how many such runs there are. The analysed flag
    # this selector reads is a property of S2 and S3 alone; it says nothing about whether a
    # classify job ever ran. Those two sets had drifted by 757 runs undetected before
    # 2026-11-02, which is why the check exists and why it is fatal rather than a warning.
    if not args.controls and not args.no_require_report:
        have = {f.name.split(".report")[0] for f in REPORTS.rglob("*.report*")}
        unclassified = [r for r, _, _ in rows if r not in have]
        if unclassified:
            import collections as _c
            bps = _c.Counter(s2[r]["BioProject"] for r in unclassified)
            sys.exit(f"{len(unclassified)} of {len(rows)} selected runs have no report, so "
                     f"they have no full-depth reference.\n  by BioProject: {dict(bps)}\n"
                     f'  classify them first:  PROJECT="{" ".join(sorted(bps))}"\n'
                     "  (run utilities/rarefaction/audit_analysed.py for the full picture)")

    # SPLIT ON WHETHER THERE IS ANYTHING TO RECOVER, measured 2026-11-02 and not anticipated
    # when this selector was written. The cohort arm's reference is each run's own full-depth
    # calls, so a run with none contributes 0/0 at every rung. Of 2,908 runs selected on depth
    # alone, 1,593 (54.8%) were in that position, and they were not spread evenly: wheat
    # 1,010/1,012 carry a detection against maize 215/974 and sorghum 90/922. Selecting on
    # depth therefore balanced the crops in the denominator and left the curve 77%
    # wheat-determined, which is the exact assumption the cohort-wide run exists to retire.
    #   The zero-detection runs are kept, in a bounded sample, because they measure the other
    # side of the floor. A run with no full-depth call that GAINS one at a shallow rung is a
    # spurious gain, so these rows give a false-discovery rate per rung where the recovery rows
    # give sensitivity. That is a better use of them than discarding them, and a few hundred
    # pin the rate down as well as 1,593 would at half the compute.
    if not args.controls:
        if not COINFECTION.exists():
            sys.exit("no coinfection.tsv; run 05_detect/detect.py first")
        nfungi = {r["run"]: int(r["n_fungi"] or 0) for r in
                  csv.DictReader(open(COINFECTION), delimiter="\t")}
        recover = [x for x in rows if nfungi.get(x[0], 0) > 0]
        zero = [x for x in rows if nfungi.get(x[0], 0) == 0]
        if args.all_zero or not args.null_sample:
            keep_zero = zero if args.all_zero else []
        else:
            # Deterministic and per-crop, so the null sample cannot itself be crop-skewed and
            # so re-running the selector returns the same set. Sorted, then strided: no RNG,
            # no seed to record, and reproducible from the inputs alone.
            keep_zero = []
            by_crop = collections.defaultdict(list)
            for x in zero:
                by_crop[x[2]].append(x)
            per = max(1, args.null_sample // max(1, len(by_crop)))
            for crop in sorted(by_crop):
                z = sorted(by_crop[crop])
                step = max(1, len(z) // per)
                keep_zero.extend(z[::step][:per])
        print(f"  recovery rows {len(recover):,}  "
              f"(by crop: {dict(collections.Counter(c for _, _, c in recover))})",
              file=sys.stderr)
        print(f"  spurious-gain rows {len(keep_zero):,} sampled from {len(zero):,} "
              f"zero-detection runs", file=sys.stderr)
        rows = recover + keep_zero

    if args.max_per_project:
        seen, keep = {}, []
        for run, bp, c in sorted(rows):
            seen[bp] = seen.get(bp, 0) + 1
            if seen[bp] <= args.max_per_project:
                keep.append((run, bp, c))
        rows = keep

    by_crop = collections.Counter(c for _, _, c in rows)
    print(f"  {len(rows):,} runs over {len({b for _, b, _ in rows})} projects: "
          f"{dict(by_crop)}", file=sys.stderr)
    print(",".join(r for r, _, _ in sorted(rows)))


if __name__ == "__main__":
    main()
