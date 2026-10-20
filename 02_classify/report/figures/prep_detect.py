#!/usr/bin/env python3
"""prep_detect.py — what 05_detect decided, and on what grounds.

Three tables for three panels:

    detect_surface.tsv   every candidate binned in (redundancy, coverage) space, with the
                         call. The DECISION SURFACE itself, not an illustration of it.
    detect_persample.tsv pathogens per sample, by crop and by sampling stratum
    detect_species.tsv   which species are actually being called, labelled by guild

Binned rather than one row per candidate because there are 300,061 of them and a scatter at
that density draws a silhouette rather than a distribution. Bins are in log space on both
axes, and the per-bin call fraction is carried so a boundary region shows as a gradient
instead of being flattened to whichever call happens to be more common there.

    python 02_classify/report/figures/prep_detect.py
"""

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import CLASSIFIED, COINFECTION, DETECTIONS  # noqa: E402
from _util import atomic_write  # noqa: E402

SURFACE = HERE / "detect_surface.tsv"
PERSAMPLE = HERE / "detect_persample.tsv"
SPECIES = HERE / "detect_species.tsv"
CAPTION = HERE / "detect_caption.txt"

NBIN = 70


def main() -> None:
    rows = list(csv.DictReader(open(CLASSIFIED), delimiter="\t"))
    det = list(csv.DictReader(open(DETECTIONS), delimiter="\t"))
    summ = list(csv.DictReader(open(COINFECTION), delimiter="\t"))
    print(f"{len(rows):,} candidates, {len(det):,} calls, {len(summ):,} runs")

    # ── a. the decision surface ─────────────────────────────────────────────────────────────
    pts = []
    for r in rows:
        try:
            x, y = float(r["redundancy"]), float(r["coverage"])
        except ValueError:
            continue
        if x <= 0 or y <= 0:
            continue
        pts.append((math.log10(x), math.log10(y), r["call"], r["reason"]))
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    xlo, xhi = min(xs), max(xs)
    ylo, yhi = min(ys), max(ys)
    bins = defaultdict(lambda: {"n": 0, "yes": 0, "piled": 0})
    for lx, ly, call, why in pts:
        i = min(NBIN - 1, int((lx - xlo) / (xhi - xlo) * NBIN))
        j = min(NBIN - 1, int((ly - ylo) / (yhi - ylo) * NBIN))
        b = bins[(i, j)]
        b["n"] += 1
        b["yes"] += call == "yes"
        b["piled"] += why == "reads piled, not spread"
    with atomic_write(SURFACE) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["x", "y", "n", "n_called", "frac_called", "n_piled"])
        for (i, j), b in sorted(bins.items()):
            x = xlo + (i + 0.5) * (xhi - xlo) / NBIN
            y = ylo + (j + 0.5) * (yhi - ylo) / NBIN
            w.writerow([f"{x:.5f}", f"{y:.5f}", b["n"], b["yes"],
                        f"{b['yes']/b['n']:.4f}", b["piled"]])

    # ── b. per sample ───────────────────────────────────────────────────────────────────────
    with atomic_write(PERSAMPLE) as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["run", "crop", "sampling_selection", "n_fungi",
                                       "n_pathogens", "n_undeclared", "coinfected"])
        w.writeheader()
        for s in summ:
            w.writerow({k: s[k] for k in w.fieldnames})

    # ── c. species ──────────────────────────────────────────────────────────────────────────
    # Guild collapsed to ONE label per species for drawing. A taxon can carry several FUNGuild
    # strings ("Endophyte-Plant Pathogen"), and precedence is pathogen first because that is
    # the claim being made; the full string stays in _classified.tsv.
    def guild_of(g: str) -> str:
        if not g:
            return "unlabelled"
        if "Plant Pathogen" in g:
            return "plant pathogen"
        if "Endophyte" in g:
            return "endophyte"
        if "Saprotroph" in g:
            return "saprotroph"
        if "Animal Pathogen" in g:
            return "animal pathogen"
        return "other"

    per = defaultdict(lambda: {"runs": set(), "guild": "", "cov": [], "sat": []})
    for r in det:
        if r["kingdom"] != "Fungi":
            continue
        e = per[r["species"]]
        e["runs"].add(r["run"])
        e["guild"] = guild_of(r["guild"])
        try:
            e["cov"].append(float(r["coverage"]))
            e["sat"].append(float(r["saturation"]) if r["saturation"] else 0.0)
        except ValueError:
            pass
    with atomic_write(SPECIES) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "guild", "n_runs", "median_coverage", "median_saturation"])
        for sp, e in sorted(per.items(), key=lambda kv: -len(kv[1]["runs"])):
            cov = sorted(e["cov"])
            sat = sorted(e["sat"])
            w.writerow([sp, e["guild"], len(e["runs"]),
                        f"{cov[len(cov)//2]:.5f}" if cov else "",
                        f"{sat[len(sat)//2]:.5f}" if sat else ""])

    # ── caption ─────────────────────────────────────────────────────────────────────────────
    rej = Counter(r["reason"] for r in rows if r["call"] == "no")
    co = sum(1 for s in summ if s["coinfected"] == "yes")
    L, A = [], None
    A = L.append
    A("What 05_detect decided")
    A("Caption for report/figures/detect.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_detect.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"(a) THE DECISION SURFACE. All {len(rows):,} species-rank candidates across "
      f"{len(summ):,} runs,")
    A("binned in log space by how much of a taxon the reads could have covered (redundancy)")
    A("against how much they did (coverage), both normalised by that taxon's own db_v4")
    A("minimizer total so taxa of any genome size share one scale. Shading is the fraction of")
    A("each bin that was CALLED, so the boundary shows as a gradient rather than a line.")
    A("")
    A("Two boundaries are drawn. The horizontal is the abundance floor (coverage >= 1%): is")
    A("the organism's DNA there at all. The curve is the shape test (saturation >= 0.05, i.e.")
    A("coverage >= 0.05 x (1-exp(-redundancy))): are the reads spread across the taxon's k-mer")
    A("space or piled onto a few loci. They answer different questions and are not")
    A("interchangeable -- that distinction is carried over from the db_v3 filter, where score")
    A("alone gave wheat 80% co-infected by admitting trace detections and the coverage floor")
    A("brought it to 47%.")
    A("")
    A(f"  candidates                        {len(rows):8,d}")
    for why, n in rej.most_common():
        A(f"  rejected, {why:30s} {n:8,d}")
    A(f"  CALLED                            {len(det):8,d}")
    A("")
    A("The third floor cannot be drawn here: an absolute 10,000 distinct minimizers is a")
    A("different line for every taxon once normalised, so it does not project onto these axes.")
    A("It is not optional -- 1% of db_v4's smallest species is 29 k-mers, reachable by")
    A("cross-mapping alone, and 20 species have ceilings under 100,000.")
    A("")
    A("(b) PATHOGENS PER SAMPLE, by crop and by sampling stratum. Plant-pathogen status is a")
    A("LABEL from holobase applied AFTER detection, never a filter before it: the db_v3 filter")
    A("skipped any taxon outside PHI-base's 632, which under db_v4 would discard every")
    A("resident the database was expanded to hold.")
    A("")
    A(f"  >=2 plant pathogens: {co:,} of {len(summ):,} runs ({co/len(summ):.0%})")
    A("")
    A("(c) WHAT IS ACTUALLY BEING CALLED, by number of runs, coloured by FUNGuild guild. A")
    A("species carrying several guild strings is drawn under the strongest claim (pathogen")
    A("before endophyte before saprotroph); the full string stays in _classified.tsv.")
    A("")
    A("WHAT THIS IS NOT. These are not prevalence estimates. Every run classified so far is")
    A("STAT-retained, because the reads on disk were downloaded through the STAT run table,")
    A("which selected for pathogen presence: projects where STAT kept <50% of runs give 83%")
    A("co-infection against 58% where it kept >=90%. The thresholds are also a provisional")
    A("operating point, as the db_v3 weights were -- the sensitivity grid spans 3,396 to")
    A("17,755 calls across plausible floors, and calibrating the saturation floor needs the")
    A("axenic controls.")
    A("")
    A("SHORT CAPTION")
    A("")
    A(f"(a) All {len(rows):,} candidate detections in redundancy-coverage space, both axes")
    A("normalised by the taxon's database k-mer total; shading is the fraction called. The")
    A("horizontal boundary is the abundance floor, the curve is the shape test. (b) Plant")
    A("pathogens per sample by crop and sampling stratum, labelled after detection. (c) The")
    A("species being called, by run count and guild.")
    A("")
    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {SURFACE.name} ({len(bins):,} bins)")
    print(f"-> {PERSAMPLE.name}  -> {SPECIES.name} ({len(per)} species)  -> {CAPTION.name}")


if __name__ == "__main__":
    main()
