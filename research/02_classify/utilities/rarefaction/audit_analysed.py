#!/usr/bin/env python3
"""audit_analysed.py -- does every analysed run actually have a report, and vice versa?

    python 02_classify/utilities/rarefaction/audit_analysed.py
    python 02_classify/utilities/rarefaction/audit_analysed.py --reads-list runs_on_scratch.txt

This is precondition 4 of RUN_COHORT.md, which carried it as a placeholder comment. It exists
because the two sets are built from completely different things and had drifted by 757 runs
without anything noticing. The analysed set is a claim made by the SUPPLEMENT: S2 says the
study design did not exclude the run and S3 says it cleared the depth floor. The report set is
a fact about DISK: a file written by a classify job that may have timed out, or never been
submitted with the environment it needed. Nothing reconciles the two, so a figure can quote an
analysed cohort of 3,541 while the calls behind it rest on 2,784.

The asymmetry matters and the script reports both directions separately:

  ANALYSED WITHOUT A REPORT is a hole in the result. Every number downstream -- the detection
  rate, the rarefaction denominator, the funnel's final bar -- is computed over whatever
  happens to be classified, so this is the gap that silently shrinks the denominator.

  A REPORT WITHOUT AN ANALYSED FLAG is not an error and is expected to be large. Runs excluded
  on study design or depth were often classified before the exclusion existed, and the reports
  are kept: deleting them would mean reclassifying from scratch if a criterion is revisited,
  which has already happened twice. They must simply never enter a denominator.

--reads-list takes a file of run accessions, one per line, as produced on Setonix by

    find .../bronze/reads -name '*.fastq.gz' -printf '%f\n' \
      | sed -E 's/_[12]\.fastq\.gz$//; s/\.fastq\.gz$//' | sort -u

With it the script also separates "no reads" from "reads but not classified", which is the
difference between resubmitting a fetch and resubmitting a classify.
"""

import argparse
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
ROOT = MODULE.parent

TABLES = ROOT / "01_literature" / "report" / "tables"
S2 = TABLES / "Supplementary Data S2.tsv"
S3 = TABLES / "Supplementary Data S3.tsv"
REPORTS = MODULE / "data" / "bronze" / "reports"
CROPS = ("Wheat", "Maize", "Rice", "Sorghum")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads-list", type=Path,
                    help="file of run accessions that have reads on disk, one per line")
    ap.add_argument("--reports", type=Path, default=REPORTS, help="reports root")
    args = ap.parse_args()

    for p in (S2, S3):
        if not p.exists():
            sys.exit(f"missing {p.name}; run report/tables/build_s2.py and build_s3.py")

    s2 = {r["Run"]: r for r in csv.DictReader(open(S2), delimiter="\t")}
    s3 = {r["Run"]: r for r in csv.DictReader(open(S3), delimiter="\t")}

    # The analysed predicate is S3's `analysed`, never `passes_depth_floor`: the latter admits
    # the one MinION run that clears the floor only because Gkmer overstates usable sequence
    # for long reads. Crop and methodological exclusion come from S2.
    analysed = {
        run for run, r in s2.items()
        if r["selected"] == "1" and r["crop"] in CROPS
        and not r["methodological_exclusion"]
        and s3.get(run, {}).get("analysed") == "1"
    }

    # A report may be <Run>.report or <Run>.report.tsv, flat or under <BioProject>/.
    have = {p.name.split(".report")[0] for p in args.reports.rglob("*.report*")}

    missing = analysed - have
    extra = have - analysed

    print(f"analysed cohort (S2 + S3)        {len(analysed):,}")
    print(f"reports on disk                  {len(have):,}")
    print(f"  analysed WITH a report         {len(analysed & have):,}")
    print(f"  analysed WITHOUT a report      {len(missing):,}"
          f"{'   <- the hole in the result' if missing else '   OK'}")
    print(f"  reports not in the cohort      {len(extra):,}   (expected; excluded runs)")

    if missing:
        by_bp = collections.Counter(s2[r]["BioProject"] for r in missing)
        print("\n  analysed but unclassified, by BioProject:")
        for bp, n in by_bp.most_common():
            print(f"    {bp:<16} {n:>5}")

        if args.reads_list:
            on = {l.strip() for l in open(args.reads_list) if l.strip()}
            no_reads = missing - on
            no_class = missing & on
            print(f"\n  of those, reads ON disk but no report  {len(no_class):,}"
                  "   -> resubmit classify")
            print(f"            reads NOT on disk            {len(no_reads):,}"
                  "   -> resubmit fetch")
            for label, s in (("classify", no_class), ("fetch", no_reads)):
                if s:
                    bps = collections.Counter(s2[r]["BioProject"] for r in s)
                    print(f"    {label:<9} PROJECT=\"{' '.join(sorted(bps))}\"")

    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
