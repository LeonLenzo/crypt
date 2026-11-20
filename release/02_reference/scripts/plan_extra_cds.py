#!/usr/bin/env python3
"""plan_extra_cds.py - the db2 addition: predicted CDS that is not in the build plan.

    python release/02_reference/scripts/plan_extra_cds.py
    python release/02_reference/scripts/plan_extra_cds.py --check

Reads  research/predict/data/worklist.tsv              the 135 genomes, accession + ncbi_taxid
       research/predict/data/predicted/<acc>/coding.fa  vipsania's CDS, deduplicated
Writes data/curation/extra_cds_db2.tsv                  accession, taxid, repo-relative path

WHY THIS IS A SEPARATE FILE AND NOT A TIER IN build_plan.tsv. `plan_build.py` emits only CDS
that exists to download from NCBI, which predicted CDS by definition does not, so the plan
cannot carry it without the plan becoming a thing that describes files it cannot find.
`build_db.py --extra-cds` exists for exactly this and takes three columns; this script is what
fills them. db2 is then the same two tiers db1' used plus this file, which is what makes
db2 - db1' a clean subtraction: one plan, one difference.

WHY curation/ AND NOT outputs/. It is a decision about what goes into a database, same as
`build_plan.tsv`, and it must be identical across any rebuild for the subtraction to hold. It
is small and tracked for that reason.

THE WORKLIST IS AUTHORITATIVE, NOT THE DIRECTORY LISTING. `predicted/` holds 137 genomes and
the worklist 135: the two extras are both *Pichia fermentans* assemblies from the
`pfermentans_alt` side-check, which was about whether an alternate assembly changed the BUSCO
score, not about entering a database. Taking the directory listing as the input would silently
add them. So the join runs the other way, and an accession in the worklist with no prediction
is a hard failure rather than a row quietly dropped.

PREDICTED CDS IS A COMPETITOR, NEVER EVIDENCE. It may absorb reads; it must never support a
positive call about the species it came from. That is a reporting rule and nothing here can
enforce it.
"""

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
REPO = MODULE.parent.parent
PREDICT = REPO / "research" / "predict" / "data"
WORKLIST = PREDICT / "worklist.tsv"
PREDICTED = PREDICT / "predicted"
OUT = MODULE / "data" / "curation" / "extra_cds_db2.tsv"
COLS = ["accession", "taxid", "path"]


def rd(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="resolve and report, write nothing")
    args = ap.parse_args()

    if not WORKLIST.exists():
        sys.exit(f"no worklist at {WORKLIST} - this reads research/predict's output, which "
                 f"only exists on Setonix")

    rows, bad = [], []
    for r in rd(WORKLIST):
        acc = (r.get("accession") or "").strip()
        taxid = (r.get("ncbi_taxid") or "").strip()
        if not acc:
            bad.append(("<no accession>", "worklist row has no accession"))
            continue
        if not taxid.isdigit():
            bad.append((acc, f"no usable ncbi_taxid ({taxid!r})"))
            continue
        cds = PREDICTED / acc / "coding.fa"
        if not cds.is_file():
            bad.append((acc, f"no prediction at {cds}"))
            continue
        if not cds.stat().st_size:
            bad.append((acc, "coding.fa is empty"))
            continue
        # REPO-RELATIVE, so the table describes the same files on both machines. The first
        # generated copy held absolute scratch paths, which made a tracked decision table
        # machine-specific; build_db.py resolves a relative path against the repo root.
        rows.append({"accession": acc, "taxid": taxid,
                     "path": str(cds.relative_to(REPO))})

    nbytes = sum((PREDICTED / r["accession"] / "coding.fa").stat().st_size for r in rows)
    print(f"worklist {len(rd(WORKLIST))} rows -> {len(rows)} resolved "
          f"({nbytes / 1e9:.2f} GB of predicted CDS)")

    extras = sorted({d.name for d in PREDICTED.iterdir() if d.is_dir()}
                    - {r["accession"] for r in rows})
    if extras:
        print(f"  {len(extras)} predicted genome(s) deliberately NOT included (absent from "
              f"the worklist): {', '.join(extras)}")

    if bad:
        print(f"\n{len(bad)} worklist row(s) could not be resolved:", file=sys.stderr)
        for acc, why in bad:
            print(f"  {acc}  {why}", file=sys.stderr)
        sys.exit(1)

    if args.check:
        print("\n--check: nothing written")
        return

    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".tsv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(OUT)
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
