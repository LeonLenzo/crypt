#!/usr/bin/env python3
"""build_s3.py -- Supplementary Data S3: library statistics and the depth filter.

    python 01_literature/report/tables/build_s3.py

WHY THIS IS A SEPARATE TABLE. Leon, 2026-10-29. S2 holds the cohort as the literature defines
it, and every exclusion there is a statement about what the authors did, traceable to a quoted
sentence. The depth filter is not that. It is a statement about what we could measure: the
same experiment sequenced deeper would pass it, and it was applied long after the cohort was
assembled. Mixing the two made a sequencing decision look like a property of the experiment.

So S3 carries the technical side: what each library actually is, and whether it cleared the
floor. One row per run, joining to S2 on Run.

THE FLOOR IS IN Gkmer, NOT READS. A k-mer cannot start within the last k-1 bases of a read,
so a library of short reads converts a smaller share of its bases into usable k-mers than a
library of long ones: at k=35, a 51 nt read yields 0.33 of its bases and a 150 nt read 0.77.
Counting reads or bases would rank a deep short-read library above a shallower long-read one
that carries more usable sequence. `gkmer` is the number the floor is applied to and
`usable_fraction` is printed beside it so the penalty is visible rather than implied.
"""

import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _util import atomic_write  # noqa: E402

SILVER = MODULE / "data" / "silver"
RELEASE = MODULE / "release"
CLASSIFY = MODULE.parent / "02_classify"
FUNNEL = CLASSIFY / "report" / "figures" / "depth_funnel" / "depth_funnel.tsv"
OUT = HERE / "Supplementary Data S3.tsv"

K = 35
# 0.5 Gkmer, set at the knee of the rarefaction curve: detection retention climbs
# steeply to 0.5 (23% of established detections at 0.019 Gkmer, 41% at 0.1, 55% at
# 0.21, 68% at 0.5) and then flattens, gaining only nine points by 2.0 Gkmer for a
# fourfold increase in depth. It was 0.463 until 2026-10-30, a value that first
# appears in TODO-RESULTS.md on 2026-10-26 with no derivation recorded anywhere and
# which matches nothing measurable: percentile 36.3 of the depth distribution, no gap
# near it, no project median. Moving to 0.5 costs 63 runs, all maize.
#   CAVEAT, and it belongs in the Methods: the curve is five WHEAT libraries, so
# applying it to maize and sorghum is an assumption. And it is a floor on TOTAL
# depth, not fungal depth, which is bimodal and cannot be used because it is an
# outcome of the classification rather than a property of the library.
FLOOR = 0.5   # Gkmer; the depth floor applied in the inclusion funnel

COLS = ["Run", "BioProject", "platform", "instrument_model", "library_layout",
        "read_length", "n_fastq", "reads", "bases", "usable_fraction", "gkmer",
        "depth_floor_gkmer", "passes_depth_floor", "analysed", "fastq_bytes"]


def main() -> None:
    cohort_bp = {r["BioProject"] for r in
                 csv.DictReader(open(RELEASE / "cohort_bioprojects.tsv"), delimiter="\t")}
    depth = {r["Run"]: r for r in
             csv.DictReader(open(CLASSIFY / "data" / "silver" / "run_depth.tsv"),
                            delimiter="\t")}
    fun = {r["Run"]: r for r in csv.DictReader(open(FUNNEL), delimiter="\t")} \
        if FUNNEL.exists() else {}
    cur = {r["Run"] for r in
           csv.DictReader(open(RELEASE / "cohort_samples.tsv"), delimiter="\t")}

    rows = []
    for r in csv.DictReader(open(SILVER / "runs.tsv"), delimiter="\t"):
        if r["BioProject"] not in cohort_bp or r["Run"] not in cur:
            continue
        d = depth.get(r["Run"], {})
        f = fun.get(r["Run"], {})
        # ENA reports read_len as bases PER SPOT, so a paired library's per-mate length is
        # half of it. Getting this backwards doubles every depth estimate.
        try:
            per_spot = float(d.get("read_len") or 0)
        except ValueError:
            per_spot = 0.0
        layout = d.get("layout", "")
        L = per_spot / 2 if layout == "PAIRED" else per_spot
        usable = (L - K + 1) / L if L > K else 0.0
        try:
            gk = float(f["gkmer"]) if f.get("gkmer") else None
        except ValueError:
            gk = None
        rows.append({
            "Run": r["Run"], "BioProject": r["BioProject"],
            "platform": d.get("platform", ""), "instrument_model": d.get("model", ""),
            "library_layout": layout,
            "read_length": f"{L:.0f}" if L else "",
            "n_fastq": d.get("n_fastq", ""),
            "reads": d.get("reads", ""), "bases": d.get("bases", ""),
            "usable_fraction": f"{usable:.3f}" if usable else "",
            "gkmer": f"{gk:.5f}" if gk is not None else "",
            "depth_floor_gkmer": FLOOR,
            "passes_depth_floor": "" if gk is None else int(gk >= FLOOR),
            # A run is analysed only if it cleared BOTH filters; the methodological one is
            # S2's and is repeated here as a single flag rather than re-litigated.
            "analysed": int(f.get("fate") == "retained"),
            "fastq_bytes": d.get("fastq_bytes", ""),
        })
    rows.sort(key=lambda x: (x["BioProject"], x["Run"]))

    with atomic_write(OUT) as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    withgk = [x for x in rows if x["gkmer"]]
    passed = [x for x in withgk if x["passes_depth_floor"] == 1]
    print(f"-> {OUT.name}  ({len(COLS)} columns)")
    print(f"   runs with library statistics      {len(rows):,}")
    print(f"   ... clearing the {FLOOR} Gkmer floor  {len(passed):,}")
    print(f"   ... analysed (both filters)       {sum(x['analysed'] for x in rows):,}")
    pl = collections.Counter(x["platform"] for x in rows if x["platform"])
    print(f"   platforms: {dict(pl.most_common(4))}")
    rl = collections.Counter(x["read_length"] for x in rows if x["read_length"])
    print(f"   read lengths (per mate): {dict(rl.most_common(6))}")
    tb = sum(int(x["fastq_bytes"]) for x in rows if (x["fastq_bytes"] or "").isdigit())
    print(f"   total library size: {tb / 1e12:.2f} TB")


if __name__ == "__main__":
    main()
