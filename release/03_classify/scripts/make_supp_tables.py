#!/usr/bin/env python3
"""Build the manuscript's detection supplementary tables from the calls.

    python release/03_classify/scripts/make_supp_tables.py
    python release/03_classify/scripts/make_supp_tables.py --reported-arm db1p

Reads  data/outputs/detections_<arm>.tsv    the calls, both arms where both exist
       data/raw/candidates_<arm>.tsv        for the per-run candidate count
       ../01_cohort/data/outputs/supp_data_s2_cohort.tsv   every cohort run, called or not
Writes data/outputs/supp_data_s5_detections.tsv  one row per call in the reported arm
       data/outputs/supp_data_s6_runs.tsv        one row per cohort run

The reported arm is **db2**, the fullest database, and S5 carries `also_in_db1p` so a reviewer
can see which calls depend on the predicted references being present as competitors. The
predicted taxa themselves are never calls; that rule lives in detect.py.

S6 is one row per COHORT run and not per classified run. A run with no call is not missing
data: it is a run whose fungal signal is diffuse, and the median such run in the db_v4 pass
still carried 139,238 fungal clade reads. Dropping them would condition the denominator on the
outcome, which is the one exclusion this project has explicitly rejected: it takes the headline
contrast from 43.2% against 8.4% to 43.3% against 32.6%, deleting the finding by filter.
Report them as "no dominant fungus", never as "no fungi".

`coinfected` is species-level and `coinfected_genus` is the genus-level read. Quote both or
neither: a genus-level count collapses congeners that may be one organism's reads spread over
several references, and a species-level count keeps them apart. Which is right depends on the
run, and 04_align is what settles it.
"""

import argparse
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
RELEASE = MODULE.parent
DATA = MODULE / "data"
COHORT = RELEASE / "01_cohort" / "data" / "outputs" / "supp_data_s2_cohort.tsv"

ARMS = ("db1p", "db2")
S5_COLS = ["run", "BioProject", "host", "setting", "arm", "taxid", "taxon", "genus",
           "taxon_class", "reads", "distinct", "ceiling", "coverage", "saturation",
           "redundancy", "judgeable", "on_wheat", "is_plant_pathogen", "is_wheat_pathogen",
           "declared", "declared_rank", "declared_status", "genus_calls", "genus_cov_cv",
           "also_in_db1p"]
S6_COLS = ["run", "BioProject", "host", "setting", "arm", "classified", "n_candidates",
           "n_calls", "n_genera", "n_pathogen_calls", "n_wheat_pathogen_calls",
           "n_undeclared_calls", "declared_pathogens", "declared_status", "coinfected",
           "coinfected_genus", "host_called", "top_taxon", "top_coverage",
           "congener_clusters"]


def read(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def yes(v: str) -> bool:
    return v in ("yes", "1", "true")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reported-arm", default="db2", choices=ARMS)
    args = ap.parse_args()
    arm = args.reported_arm
    other = "db1p" if arm == "db2" else "db2"

    det_path = DATA / "outputs" / f"detections_{arm}.tsv"
    if not det_path.is_file():
        sys.exit(f"no {det_path.name}; run detect.py --arm {arm} first")
    dets = read(det_path)

    other_path = DATA / "outputs" / f"detections_{other}.tsv"
    other_calls = set()
    if other_path.is_file():
        other_calls = {(d["run"], d["taxid"]) for d in read(other_path)}
    else:
        print(f"note: {other_path.name} not present, so also_in_{other} is left blank")

    # ---- S5: one row per call ----
    s5 = []
    for d in dets:
        s5.append({**d, f"also_in_{other}": ("" if not other_calls else
                                             int((d["run"], d["taxid"]) in other_calls))})
    cols = [c if c != "also_in_db1p" else f"also_in_{other}" for c in S5_COLS]
    out5 = DATA / "outputs" / "supp_data_s5_detections.tsv"
    with open(out5, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(s5)

    # ---- S6: one row per cohort run ----
    by_run = collections.defaultdict(list)
    for d in dets:
        by_run[d["run"]].append(d)
    cand_path = DATA / "raw" / f"candidates_{arm}.tsv"
    n_cand = collections.Counter()
    classified = set()
    if cand_path.is_file():
        for c in read(cand_path):
            n_cand[c["run"]] += 1
            classified.add(c["run"])
    else:
        print(f"note: {cand_path.name} not present (it lives on scratch), so n_candidates "
              "is blank and `classified` is read off the detections only")
        classified = set(by_run)

    s6 = []
    for r in read(COHORT):
        run = r["run"]
        calls = [d for d in by_run.get(run, []) if d["taxon_class"] != "host"]
        genera = {d["genus"] for d in calls}
        path_calls = [d for d in calls if yes(d["is_plant_pathogen"])]
        wheat_calls = [d for d in calls if yes(d["is_wheat_pathogen"])]
        undeclared = [d for d in calls if d["declared"] == "undeclared"]
        clusters = sorted({d["genus"] for d in calls if int(d["genus_calls"] or 1) >= 3})
        top = max(calls, key=lambda d: float(d["coverage"]), default=None)
        s6.append({
            "run": run, "BioProject": r["BioProject"], "host": r.get("host", ""),
            "setting": r.get("setting", ""), "arm": arm,
            "classified": int(run in classified),
            "n_candidates": n_cand.get(run, "") if cand_path.is_file() else "",
            "n_calls": len(calls), "n_genera": len(genera),
            "n_pathogen_calls": len(path_calls),
            "n_wheat_pathogen_calls": len(wheat_calls),
            "n_undeclared_calls": len(undeclared),
            "declared_pathogens": r.get("declared_pathogens", ""),
            "declared_status": r.get("declared_status", ""),
            "coinfected": int(len(path_calls) >= 2),
            "coinfected_genus": int(len({d["genus"] for d in path_calls}) >= 2),
            "host_called": int(any(d["taxon_class"] == "host" for d in by_run.get(run, []))),
            "top_taxon": top["taxon"] if top else "",
            "top_coverage": top["coverage"] if top else "",
            "congener_clusters": ";".join(clusters),
        })
    out6 = DATA / "outputs" / "supp_data_s6_runs.tsv"
    with open(out6, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=S6_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(s6)

    done = [r for r in s6 if r["classified"]]
    withcall = [r for r in done if r["n_calls"]]
    print(f"-> {out5.relative_to(MODULE)}  {len(s5):,} calls in {arm}")
    print(f"-> {out6.relative_to(MODULE)}  {len(s6):,} cohort runs, {len(done):,} classified, "
          f"{len(withcall):,} with a fungal call")
    print(f"   runs with a call and a curated declaration: "
          f"{sum(1 for r in withcall if r['declared_status'] != 'not_curated'):,}; "
          f"with an uncurated declaration: "
          f"{sum(1 for r in withcall if r['declared_status'] == 'not_curated'):,}")
    print("No rate out of these tables is a prevalence. See CLAUDE.md.")


if __name__ == "__main__":
    main()
