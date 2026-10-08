#!/usr/bin/env python3
"""The STAT-frame assessment queue: which classified projects still need a verdict.

Written 2026-10-08. The STAT frame (01_stat) screened 10,995 runs and 03_kraken classified
a subset of them; those runs are the ones whose co-infection calls are already paid for, so
a `setting` verdict on their project is worth more than a verdict on an unclassified one.
This lists the projects that carry classified runs and have no setting yet, newest-evidence
first, so the expensive tail can be worked in a defined order rather than rediscovered.

Three facts per project, because they decide how to read it:
  classified  - runs already classified by 03_kraken (what a verdict buys)
  group       - host group from _hostgroup (cereal first; `cereal?` is pathogen-inferred)
  text        - whether a full text is on disk, so `sources.py cached DOI` can open it
                without a network call

A project is "pending" when the literature frame's silver `runs.tsv` carries no `setting`
for any of its runs and the root `_exclusions.EXCLUDED` registry does not hold it. The test
has to be on the DERIVED table, not on curation.tsv: the eight resolver-curated projects
(Ada21's 539 classified runs among them) get their setting from `joins.tsv` or a resolver
writing runs.tsv directly, and never acquire a `curation.tsv` rule. Testing curation.tsv
alone put those back in the queue as unassessed.

Run:  python 01a_Literature/report/stat_queue.py            # pending cereal projects
      python 01a_Literature/report/stat_queue.py --group all --done
"""

from __future__ import annotations

import argparse, collections, csv, json, sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _exclusions import is_excluded, reason_for
from _hostgroup import ORDER, project_group
from _layout import ROOT, SEEDS, SILVER, TEXT_CACHE

STAT_RUNS  = ROOT / "01_stat/03_filter/data/runs.tsv"
CLASSIFIED = ROOT / "03_kraken/05_filter/data/_classified.tsv"


def tsv(path: Path) -> list[dict]:
    if not path.exists():
        sys.exit(f"missing: {path}")
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="cereal",
                    help="cereal (incl. cereal?), grass, other, or all")
    ap.add_argument("--done", action="store_true", help="list the settled ones instead")
    ap.add_argument("--min", type=int, default=1, help="minimum classified runs")
    a = ap.parse_args()

    # Which BioProject each classified run belongs to, and its submitted organism.
    bp_of, organism = {}, {}
    for r in tsv(STAT_RUNS):
        bp_of[r["Run"]] = r["BioProject"]
        organism.setdefault(r["BioProject"], r.get("library_organism") or r.get("host") or "")

    # Distinct RUNS, not classification rows: _classified.tsv carries one row per
    # species call, so 24,467 rows cover only 3,141 runs.
    runs = collections.defaultdict(set)
    hosts = collections.defaultdict(set)
    for r in tsv(CLASSIFIED):
        bp = bp_of.get(r["run"])
        if bp:
            runs[bp].add(r["run"])
            hosts[bp].add(r.get("host", ""))

    settled = {r["BioProject"] for r in tsv(SILVER / "runs.tsv") if r.get("setting")}

    # A DOI we hold full text for means the read costs nothing but tokens.
    have_text = set()
    if TEXT_CACHE.exists():
        with TEXT_CACHE.open() as fh:
            for line in fh:
                doi, _, blob = line.partition("\t")
                if blob.strip():
                    have_text.add(doi.strip().lower())
    papers = collections.defaultdict(set)
    ap_tsv = SILVER / "accession_papers.tsv"
    if ap_tsv.exists():
        for r in tsv(ap_tsv):
            bp, doi = r.get("accession", ""), (r.get("doi") or "").strip().lower()
            if bp and doi:
                papers[bp].add(doi)

    rows = []
    for bp, rs in runs.items():
        n = len(rs)
        if n < a.min:
            continue
        if ((bp in settled) or is_excluded(bp)) != a.done:
            continue
        g, _why = project_group(organism.get(bp, ""), *sorted(hosts[bp]))
        if a.group != "all":
            want = {"cereal": {"cereal", "cereal?"}}.get(a.group, {a.group})
            if g not in want:
                continue
        dois = sorted(papers[bp])
        cached = [d for d in dois if d in have_text]
        rows.append((ORDER.get(g, 9), -n, bp, g, n,
                     organism.get(bp, "")[:38], len(dois), len(cached)))

    rows.sort()
    label = "settled" if a.done else "pending"
    print(f"{len(rows)} {label} projects  "
          f"{sum(r[4] for r in rows)} classified runs  group={a.group}\n")
    print(f"{'BioProject':<14}{'group':<9}{'cls':>5}  {'doi':>3} {'txt':>3}  organism")
    for _, _, bp, g, n, org, nd, nc in rows:
        print(f"{bp:<14}{g:<9}{n:>5}  {nd:>3} {nc:>3}  {org}")


if __name__ == "__main__":
    main()
