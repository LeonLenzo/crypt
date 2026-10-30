#!/usr/bin/env python3
"""prep_browser.py -- pack the per-run detections into one JSON the HTML viewer reads.

    python 02_classify/report/figures/detect_runs/prep_browser.py

WHY. Rendering the cohort to PNG took 45 minutes and 472 MB, and a 23-species change to
taxon_status.tsv invalidated every image the same afternoon. The panel is cheap to draw and
expensive to store, so ship the DATA and draw in the browser. Leon, 2026-10-29.

The output is ~7 MB for 2,862 runs, which opens from disk with no server, no R and no
install -- the viewer is one HTML file beside it.

WHAT IS KEPT, and the one thing that is not:

  * Every DECISION-RELEVANT row in full: called, or cleared the absolute floor. ~12 per run.
    These are the rows where the call was actually a decision, so the viewer can re-run the
    floors on them client-side and show what a different threshold would do.
  * The deep rejected cloud, DOWNSAMPLED to CLOUD_N per run and carrying only (x, y). It is
    there to show that the floors reject far more than they admit, which is a shape, not a
    set of identities. Downsampling is stated in the viewer rather than hidden.

Numbers are rounded, but NOT on the rows a floor is applied to. The cloud goes to 3
significant figures because nothing is decided from it; the decision rows go to 6, because at
4 a coverage of 0.0199953 rounds to 0.02 and passes a 0.02 floor it should fail. A rounding
that can move a value across a threshold is a wrong answer, not a smaller file.
"""

import csv
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _util import atomic_write  # noqa: E402

SRC = HERE / "detect_runs_all_proposed.tsv"
if not SRC.exists():
    SRC = HERE / "detect_runs_proposed.tsv"
OUT = HERE / "browser" / "data.json"
# The TEMPLATE is _template.html and the thing you open is index.html, which is generated
# with the data inlined. That way round on purpose: index.html is what anyone double-clicks,
# and the template fetching data.json over file:// is blocked by Chrome as a cross-origin
# request from a null origin. Naming the working file index.html removes the choice.
TPL = HERE / "browser" / "_template.html"
SELF = HERE / "browser" / "index.html"

CLOUD_N = 140          # rejected points kept per run, for the shape of the cloud
MIN_DISTINCT = 10_000  # 05_detect's absolute floor
COV_FLOOR = 0.02
SAT_FLOOR = 0.05
SEED = 20261029        # the downsample is seeded, so the picture is reproducible

TIERS = ["host", "pathogen_crop", "pathogen_other",
         "resident_crop", "resident_other", "no_record"]
TROPHIC = ["", "biotroph", "necrotroph", "weak", "opportunist"]
EVID = ["", "A", "B", "C"]


def r6(x):
    """Six significant figures on the DECISION rows, and the reason is a floor boundary.

    At four, coverage 0.0199953 rounds to 0.02 and then PASSES a 0.02 floor that detect.py
    failed it on, so the viewer called one run differently from the pipeline. Rounding that
    can move a value across a threshold is not a size optimisation, it is a wrong answer.
    The cloud stays coarse (3 s.f.) because nothing is decided from it.
    """
    try:
        return float(f"{float(x):.6g}")
    except (TypeError, ValueError):
        return None


def main() -> None:
    rng = random.Random(SEED)
    runs, cur, dec, cloud = {}, None, [], []

    def flush():
        if cur is None:
            return
        keep = cloud if len(cloud) <= CLOUD_N else rng.sample(cloud, CLOUD_N)
        runs[cur["r"]] = {**cur, "d": dec[:], "k": keep, "kn": len(cloud)}

    for row in csv.DictReader(open(SRC), delimiter="\t"):
        try:
            red, cov = float(row["redundancy"]), float(row["coverage"])
        except ValueError:
            continue
        if red <= 0 or cov <= 0:
            continue
        if row["run"] != (cur or {}).get("r"):
            flush()
            dec, cloud = [], []
            cur = {"r": row["run"], "p": row["BioProject"], "c": row.get("crop", ""),
                   "h": row.get("host", ""), "loc": row.get("location", ""),
                   "y": row.get("year", ""), "t": row.get("tissue", ""),
                   "dp": row.get("declared_pathogens", ""),
                   "np": int(row.get("n_pathogens") or 0),
                   "nu": int(row.get("n_undeclared") or 0)}
        try:
            dist = float(row["distinct"])
        except ValueError:
            dist = 0.0
        if row["call"] == "yes" or dist >= MIN_DISTINCT:
            dec.append([
                row["species"],
                TIERS.index(row["tier"]) if row["tier"] in TIERS else 5,
                TROPHIC.index(row["trophic"]) if row.get("trophic") in TROPHIC else 0,
                EVID.index(row["evidence"]) if row.get("evidence") in EVID else 0,
                1 if row.get("declared") == "declared" else 0,
                r6(cov), r6(row["saturation"]), r6(red), int(dist),
                int(float(row["reads"] or 0)),
                # CONGENER-ONLY, carried so the viewer's re-call reproduces detect.py exactly.
                # detect.py applies this AFTER the three floors -- a taxon whose only
                # membership route is "shares a genus with a pathogen" is not reportable -- and
                # a re-call knowing only the floors called 8 runs differently. Found by diffing
                # the viewer against the pipeline, not by reading the code.
                1 if "congener" in (row.get("reason") or "") else 0,
            ])
        else:
            cloud.append([float(f"{red:.3g}"), float(f"{cov:.3g}")])
    flush()

    order = sorted(runs.values(), key=lambda x: (-x["np"], x["r"]))
    payload = {
        "meta": {"cov_floor": COV_FLOOR, "sat_floor": SAT_FLOOR,
                 "min_distinct": MIN_DISTINCT, "cloud_n": CLOUD_N,
                 "tiers": TIERS, "trophic": TROPHIC, "evid": EVID,
                 "source": SRC.name, "n_runs": len(order)},
        "runs": order,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with atomic_write(OUT) as fh:
        json.dump(payload, fh, separators=(",", ":"))

    # ── the self-contained file ────────────────────────────────────────────────────────────
    # A file:// page cannot fetch: Chrome blocks it as a cross-origin request from a null
    # origin, which is what happens when you double-click the template. So the payload is
    # inlined into a second file that opens straight from disk with nothing running.
    # `</script>` cannot appear in JSON-encoded data, but escape the slash anyway rather than
    # rely on that -- a species name is attacker-free here, but the habit is cheap.
    if TPL.exists():
        blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
        html = TPL.read_text().replace(
            "<script>",
            '<script id="payload" type="application/json">' + blob + "</script>\n<script>", 1)
        with atomic_write(SELF) as fh:
            fh.write(html)
        print(f"-> {SELF.relative_to(MODULE.parent)}: "
              f"{SELF.stat().st_size / 1e6:.1f} MB — OPEN THIS, straight from disk")

    mb = OUT.stat().st_size / 1e6
    nd = sum(len(x["d"]) for x in order)
    nk = sum(len(x["k"]) for x in order)
    tr = sum(x["kn"] for x in order)
    print(f"-> {OUT.relative_to(MODULE.parent)}: {mb:.1f} MB, {len(order):,} runs")
    print(f"   decision rows kept in full : {nd:,}")
    print(f"   cloud points kept          : {nk:,} of {tr:,} "
          f"({nk / tr * 100:.0f}%, capped at {CLOUD_N}/run)")


if __name__ == "__main__":
    main()
