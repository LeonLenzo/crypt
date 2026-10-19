#!/usr/bin/env python3
"""fetch_depth.py — how much sequence each run actually holds, from ENA.

Caches one TSV per BioProject under bronze, then writes a single joined table. Run it before
`resolve.py`, which joins the result into `run_list.tsv`.

## Why this is a step and not a column someone looks up later

`runs.tsv` carries `spots`, and a spot is not a base. The cohort's libraries run from 51 nt
single-end to 300 nt paired, so spots differ from bases by a factor that varies 12-fold across
projects, and every question about detection sensitivity is about BASES. Measured over the
7,567 fetch-eligible runs, project median depth spans **300x**:

    PRJDB7234      1,557 runs  SINGLE   51 nt   median 0.057 Gb   Rice
    PRJNA726716      104 runs  SINGLE   51 nt   median 0.109 Gb   Rice
    PRJEB39201       538 runs  PAIRED  202 nt   median 4.533 Gb   Wheat
    PRJNA1181034      12 runs  PAIRED  300 nt   median 17.12 Gb   Wheat

1,571 runs (21% of the cohort) hold under 0.2 Gb. That is not a curiosity, it is a confound on
two of the project's central comparisons, and it is visible only once bases are in the table:

  * **Crop.** PRJDB7234 and PRJNA726716 are both rice and are 1,661 of the 2,784 rice runs, so
    60% of rice is ultra-shallow. A low rice co-infection rate would have depth as a sufficient
    explanation, indistinguishable from biology.
  * **Evidence standard.** The pathology stratum medians 4.67 Gb against 1.46 Gb for
    unselected. The stratum expected to show a SECOND pathogen is also the deeper one.

So depth belongs beside the run, in the plan, before anything is fetched or classified. What to
DO about it is a 05_detect question (a per-run detection limit, or strata, or a minimum depth
for rate claims) and is deliberately not decided here; this step only makes the number visible.

## Mechanics

ONE REQUEST PER BIOPROJECT, never per run. The same lesson as 01_literature's portal fetch:
per-run requests were ~150x slower and returned less. 41 requests cover 7,567 runs.

EBI shares a 2 req/s limiter across EuropePMC and ENA, so requests are spaced and retried.
Responses are cached per project, because this table is regenerated whenever the cohort moves
and re-fetching 41 projects to learn nothing new is just rude to EBI.

    python 02_classify/03_cohort/fetch_depth.py
    python 02_classify/03_cohort/fetch_depth.py --refresh   # ignore the cache
"""

import argparse
import csv
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import BRONZE, COHORT_SAMPLES, SILVER  # noqa: E402

CACHE = BRONZE / "ena_depth"
DEPTH = SILVER / "run_depth.tsv"

FIELDS = ("run_accession,library_layout,library_strategy,library_selection,"
          "read_count,base_count,instrument_platform,instrument_model,fastq_ftp,fastq_bytes")
PORTAL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
          "?accession={acc}&result=read_run&fields=" + FIELDS + "&format=tsv")

COLS = ["Run", "BioProject", "layout", "strategy", "selection", "reads", "bases", "read_len",
        "platform", "model", "n_fastq", "fastq_bytes", "fastq_urls"]

PAUSE = 0.6          # EBI shares a 2 req/s limiter; stay well under it
RETRIES = 3


def fetch_project(proj: str, refresh: bool) -> list[dict]:
    cached = CACHE / f"{proj}.tsv"
    if cached.exists() and not refresh:
        return list(csv.DictReader(open(cached), delimiter="\t"))
    url = PORTAL.format(acc=proj)
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "crypt/02_classify"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                body = fh.read().decode()
            break
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
            time.sleep(3 * (attempt + 1))
    else:
        print(f"  {proj}: FAILED after {RETRIES} attempts ({last})", file=sys.stderr)
        return []
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = cached.with_suffix(".tmp")
    tmp.write_text(body)
    tmp.replace(cached)
    time.sleep(PAUSE)
    return list(csv.DictReader(body.splitlines(), delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="ignore the per-project cache")
    args = ap.parse_args()

    cohort = list(csv.DictReader(open(COHORT_SAMPLES), delimiter="\t"))
    want = defaultdict(set)
    for r in cohort:
        want[r["BioProject"]].add(r["Run"])

    out, missing = [], []
    for i, proj in enumerate(sorted(want), 1):
        rows = fetch_project(proj, args.refresh)
        byrun = {r["run_accession"]: r for r in rows}
        got = 0
        for run in sorted(want[proj]):
            r = byrun.get(run)
            if not r:
                missing.append((proj, run))
                continue
            try:
                reads, bases = int(r["read_count"] or 0), int(r["base_count"] or 0)
            except ValueError:
                reads, bases = 0, 0
            urls = [u.strip() for u in (r.get("fastq_ftp") or "").split(";") if u.strip()]
            out.append({
                "Run": run, "BioProject": proj,
                "layout": r.get("library_layout", ""),
                "strategy": r.get("library_strategy", ""),
                "selection": r.get("library_selection", ""),
                "reads": reads, "bases": bases,
                "read_len": f"{bases / reads:.1f}" if reads else "",
                "platform": r.get("instrument_platform", ""),
                "model": r.get("instrument_model", ""),
                "n_fastq": len(urls),
                "fastq_bytes": r.get("fastq_bytes", ""),
                # Kept so fetch_reads.py needs no second round-trip to ENA.
                "fastq_urls": ";".join(urls),
            })
            got += 1
        print(f"  [{i:2d}/{len(want)}] {proj:14s} {got}/{len(want[proj])} runs", flush=True)

    SILVER.mkdir(parents=True, exist_ok=True)
    tmp = DEPTH.with_suffix(".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    tmp.replace(DEPTH)
    print(f"\n-> {DEPTH.name}  ({len(out)} runs)")
    if missing:
        print(f"{len(missing)} cohort runs absent from ENA: {missing[:5]}")
        print("  ENA mirrors SRA, so this means the run is suppressed or very new. "
              "They stay in the cohort and simply have no depth recorded.")
    noseq = [r for r in out if not r["bases"]]
    if noseq:
        print(f"{len(noseq)} runs report ZERO bases: "
              f"{[r['Run'] for r in noseq][:5]} -- withdrawn or empty submissions.")


if __name__ == "__main__":
    main()
