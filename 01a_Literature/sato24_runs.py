#!/usr/bin/env python3
"""Resolve the Sato 2024 Arabidopsis field cluster, and deduplicate it.

Eight BioProjects hold one experiment: 199 A. thaliana accessions plus Col-0, in randomised
blocks, in outdoor gardens at Zurich and Otsu, in 2017 and 2018. 4,796 runs. The cluster has
no paper of its own, so the evidence comes from two:

  Sato et al. 2024, Nat Commun, 10.1038/s41467-024-52374-7 (PMC11458863)
      The field experiment: both sites with coordinates, the potting and transfer, the design,
      the herbivore surveys. Contains NO RNA-Seq and deposits no sequence.
  Tomita et al. v3, bioRxiv 10.1101/2025.05.29.656841
      The sequencing: leaf harvested three weeks post-transplantation within two hours of
      solar noon, Lasy-Seq v1.1, both platforms. Names all eight accessions and states
      plainly that it REUSED the data rather than generating it.

## The duplication, which is the point of this script

Every library was sequenced twice. Tom26 v3: "These libraries were sequenced by Illumina
HiSeq 2500 and HiSeq X platforms to generate 50-bp and 150-bp reads, from the 3' end,
respectively." The archive registered the two runs under two BioSamples, so:

    LibraryName  overlaps 100% between the platform-paired projects
    BioSample    overlaps   0%

Deduplicating on BioSample therefore counts every leaf twice, and `biosample_representative`
downstream would not catch it. The real unit is LibraryName; each of the 2,398 appears in
exactly two projects, once per platform.

## Which run to keep

Ranking on `spots` picks the HiSeq2500 run for 159 of the 2,398 pairs, which is an artefact
of ignoring read length. What matters for k-mer classification is how many k-mers a library
yields, and at k=31 a 150-bp read gives 120 while a 50-bp read gives 20. Ranked on
classifiable k-mers (reads x (readlen - k + 1)) the HiSeqX run wins all 2,398 pairs
unanimously, and carries about 10x the k-mer content. So the rule is simply: keep HiSeqX.

The discarded runs are flagged `library_representative = no`, never deleted. They are real
sequencing of real samples and remain available as a technical-replicate check.

Run:  python 01a_Literature/sato24_runs.py
"""

from __future__ import annotations

import collections, csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUNS = HERE / "data" / "runs.tsv"
STUDY = HERE / "studies" / "doi_10.1038_s41467-024-52374-7"
OUT = STUDY / "sato24_runs.csv"

PROJECTS = {
    # BioProject:      (site key, year, platform) - from the NCBI project titles, each of which
    # names its site, year and instrument explicitly.
    "PRJNA1055060": ("zurich", "2017", "2500"),
    "PRJNA1055734": ("zurich", "2017", "X"),
    "PRJNA1055104": ("zurich", "2018", "2500"),
    "PRJNA1056126": ("zurich", "2018", "X"),
    "PRJNA1055317": ("otsu", "2017", "2500"),
    "PRJNA1055736": ("otsu", "2017", "X"),
    "PRJNA1055424": ("otsu", "2018", "2500"),
    "PRJNA1055939": ("otsu", "2018", "X"),
}

# Sato 2024, Methods "Field setting". Coordinates and institution are the paper's own.
SITES = {
    "zurich": "Switzerland: Zurich, University of Zurich-Irchel campus outdoor garden",
    "otsu": "Japan: Shiga, Otsu, Center for Ecological Research, Kyoto University",
}

READ_LEN = {"2500": 50, "X": 150}   # Tom26 v3, respectively
K = 31

# Potted plants raised indoors for 1.5 months, then outdoors for three weeks under the
# natural herbivore community. Leon's call 2026-10-06 that this counts as field, carrying the
# nuance in the label as Hashida's transferred seedlings do.
SETTING = "field (potted, outdoor common garden)"
SELECTION = "unselected"


def kmers(spots: int, plat: str) -> int:
    return spots * max(0, READ_LEN[plat] - K + 1)


def main() -> None:
    runs = [r for r in csv.DictReader(RUNS.open(), delimiter="\t")
            if r["BioProject"] in PROJECTS]
    if not runs:
        sys.exit("no runs from the eight projects found in runs.tsv")
    print(f"{len(runs)} runs across {len({r['BioProject'] for r in runs})} projects",
          file=sys.stderr)

    by_lib = collections.defaultdict(list)
    for r in runs:
        site, year, plat = PROJECTS[r["BioProject"]]
        r["_site"], r["_year"], r["_plat"] = site, year, plat
        r["_kmers"] = kmers(int(r["spots"] or 0), plat)
        by_lib[(site, year, r["LibraryName"])].append(r)

    sizes = collections.Counter(len(v) for v in by_lib.values())
    if set(sizes) != {2}:
        sys.exit(f"REFUSED: expected every library to appear exactly twice, got sizes {dict(sizes)}. "
                 f"The platform pairing is not what the methods describe.")
    print(f"  {len(by_lib)} libraries, each appearing exactly twice", file=sys.stderr)

    rows, chosen_plat = [], collections.Counter()
    for (site, year, lib), pair in by_lib.items():
        best = max(pair, key=lambda r: r["_kmers"])
        chosen_plat[best["_plat"]] += 1
        for r in pair:
            rows.append(dict(
                Run=r["Run"], LibraryName=lib, Platform=r["_plat"],
                Site=site, Year=year,
                Setting=SETTING, Location=SITES[site], Tissue="leaf",
                SamplingSelection=SELECTION,
                LibraryRepresentative="yes" if r["Run"] == best["Run"] else "no",
                Spots=r["spots"], ClassifiableKmers=r["_kmers"],
            ))

    # The methods say HiSeqX is the deeper arm; if the data ever stopped agreeing, the
    # deduplication would be silently picking a different thing than the docstring claims.
    if chosen_plat.get("2500"):
        print(f"  NOTE: {chosen_plat['2500']} libraries chose the HiSeq2500 run on k-mer "
              f"content, against the expected HiSeqX", file=sys.stderr)
    else:
        print(f"  all {chosen_plat['X']} picks are HiSeqX, as the read lengths predict",
              file=sys.stderr)

    rows.sort(key=lambda r: r["Run"])
    STUDY.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    keep = sum(1 for r in rows if r["LibraryRepresentative"] == "yes")
    print(f"\nwrote {OUT.relative_to(ROOT)}: {len(rows)} runs, {keep} representative",
          file=sys.stderr)
    for site in sorted(SITES):
        for year in ("2017", "2018"):
            n = sum(1 for r in rows if r["Site"] == site and r["Year"] == year
                    and r["LibraryRepresentative"] == "yes")
            print(f"  {site:<8} {year}  {n:>5} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
