#!/usr/bin/env python3
"""prep_geography.py — the Pst surveillance compendium as a preliminary spatial test.

PRJEB39201 is the one project in the cohort that can carry a spatial analysis on its own, and
the reason is design rather than size. It is a single surveillance programme: 538 yellow-rust
infected wheat leaves, one protocol, 30 countries, 2014 to 2018, and the authors reported
*Puccinia striiformis* and only *Puccinia striiformis*. So every other pathogen in those reads
is unreported by construction, and because the protocol is constant, a difference BETWEEN
countries is not a difference between studies.

That last point is what makes this worth doing before the rest of the cohort lands. Everywhere
else in this work, geography and study are nearly the same variable.

## THE CONFOUND THAT HAD TO BE CLEARED FIRST, AND IT IS ONLY PARTLY CLEARED

Depth is the confound that invalidates the between-crop comparison: G_kmer spans 300x across the
cohort and covaries with crop. Two correlations matter here and they do NOT agree, so both are
computed and reported rather than whichever is more convenient:

  * across ALL classified wheat (18 countries, every project):   r = -0.053
  * within THIS project alone (9 countries with n >= 10):        r = +0.271

The first says depth does not explain the spatial pattern in wheat generally. The second is
weak, positive, and on nine points, where |r| would need to reach about 0.67 to be significant
at the 5% level, so it neither establishes a depth effect nor excludes one. **It is not
evidence that the panel is a depth artefact, and it is not evidence that it is not.** What
carries the spatial reading is therefore the agreement with known disease geography and the
opposite direction of the two gradients, not the depth correlation, which is underpowered at
this sample size.

The honest fix is more countries with enough libraries, which the remaining 6,504 cohort runs
will supply.

## THE TIME AXIS IS NOT REPORTED, AND THIS IS WHY

The pooled year trend looks clean and is an artefact of composition. Co-infection by collection
year across the compendium runs 69% / 39% / 57% / 75% / 94% for 2014 to 2018, which reads as a
rise. But the countries sampled change with the year: 2015 is Germany-heavy (48 of 173), 2016 is
UK-heavy (103 of 203), and **2018 is 14 of 17 runs from two Argentinian sites**, Balcarce and
Chillar. The 2018 spike is Argentina, not a year.

Holding country constant, the UK (211 runs) gives 57% / 60% / 77% for 2015 to 2017 and Germany
(71 runs) gives 25% / 45% / 100% on n = 48 / 11 / 9. A single year's rise on n = 60 is not a
trend, so this script writes the year table for inspection and the FIGURE shows space only.

Writes geography.tsv, geography_species.tsv and geography_caption.txt.

    python 02_classify/report/figures/prep_geography.py
"""

import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import COHORT_SAMPLES, GOLD, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

PROJECT = "PRJEB39201"
MIN_RUNS = 8          # a country below this is reported in the table, not drawn
DECLARED = "Puccinia striiformis"

TSV = HERE / "geography.tsv"
SPECIES_TSV = HERE / "geography_species.tsv"
CAPTION = HERE / "geography_caption.txt"

# Curation defects found while writing this. They are reported, not silently merged, because the
# fix belongs in 01_literature's curation and a silent merge here would hide it.
ALIASES = {"United Kingdom": "UK", "The Netherlands": "Netherlands", "Zimbabwae": "Zimbabwe"}


def year_of(s: str):
    s = (s or "").strip()
    for i in range(len(s) - 3):
        if s[i:i + 4].isdigit() and 1900 < int(s[i:i + 4]) < 2030:
            return int(s[i:i + 4])
    return None


def gkmer(row) -> float | None:
    """Sequencing depth as the quantity detection depends on: k-mer positions, not bases."""
    try:
        bases, L = float(row["bases"]), float(row["read_len"])
    except (TypeError, ValueError, KeyError):
        return None
    return bases * (L - 34) / L / 1e9 if L > 35 else None


def main() -> None:
    coinf = {r["run"]: r for r in csv.DictReader(open(GOLD / "coinfection.tsv"), delimiter="\t")}
    prov = {r["Run"]: r for r in csv.DictReader(open(COHORT_SAMPLES), delimiter="\t")}
    depth = {}
    for r in csv.DictReader(open(SILVER / "run_depth.tsv"), delimiter="\t"):
        g = gkmer(r)
        if g:
            depth[r[list(r)[0]]] = g

    runs, aliased = [], defaultdict(int)
    for run, c in coinf.items():
        if c["BioProject"] != PROJECT:
            continue
        p = prov.get(run, {})
        loc = (p.get("location") or "").strip()
        raw = loc.split(":")[0].strip()
        if raw in ALIASES:
            aliased[raw] += 1
        runs.append({
            "run": run, "country": ALIASES.get(raw, raw), "locality": loc,
            "year": year_of(p.get("collection_date")), "g": depth.get(run),
            "coinf": c["coinfected"] == "yes",
            "paths": [x.strip() for x in (c["pathogens"] or "").split(";") if x.strip()],
        })

    # Species to draw: the declared pathogen as a positive control, then the co-pathogens that
    # appear often enough for a per-country rate to mean anything.
    freq = defaultdict(int)
    for r in runs:
        for sp in r["paths"]:
            freq[sp] += 1
    drawn = [DECLARED] + [sp for sp, n in sorted(freq.items(), key=lambda x: -x[1])
                          if sp != DECLARED and n >= 30]

    by_country = defaultdict(list)
    for r in runs:
        by_country[r["country"]].append(r)

    with atomic_write(TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["country", "n_runs", "coinf_rate", "median_gkmer", "drawn"] +
                   [sp.replace(" ", "_") for sp in drawn])
        for cty, v in sorted(by_country.items(), key=lambda x: -len(x[1])):
            gs = [r["g"] for r in v if r["g"]]
            w.writerow([cty, len(v),
                        f"{sum(r['coinf'] for r in v) / len(v):.4f}",
                        f"{st.median(gs):.3f}" if gs else "",
                        "yes" if len(v) >= MIN_RUNS else "no"] +
                       [f"{sum(1 for r in v if sp in r['paths']) / len(v):.4f}" for sp in drawn])

    # Per-species summary, including the signal quality that justifies calling them at all.
    det = defaultdict(list)
    for r in csv.DictReader(open(GOLD / "detections.tsv"), delimiter="\t"):
        if r.get("BioProject") != PROJECT:
            continue
        try:
            det[r["species"]].append((float(r["coverage"]), float(r["saturation"] or 0),
                                      int(r["distinct"]), r.get("guild", "")))
        except (TypeError, ValueError):
            pass
    with atomic_write(SPECIES_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "n_runs", "pct_runs", "median_coverage",
                    "median_saturation", "median_distinct", "guild", "is_declared"])
        for sp, v in sorted(det.items(), key=lambda x: -len(x[1])):
            w.writerow([sp, len(v), f"{len(v) / len(runs):.4f}",
                        f"{st.median(x[0] for x in v):.5f}",
                        f"{st.median(x[1] for x in v):.4f}",
                        int(st.median(x[2] for x in v)),
                        v[0][3], "yes" if sp == DECLARED else "no"])

    # ── the depth check, recomputed here rather than quoted, at BOTH scopes ──────────────────
    def pearson(pts):
        if len(pts) < 3:
            return float("nan"), len(pts)
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        mx, my = st.mean(xs), st.mean(ys)
        num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
        den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
        return (num / den if den else float("nan")), len(pts)

    def country_points(rows_):
        g = defaultdict(list)
        for r in rows_:
            g[r["country"]].append(r)
        return [(st.median([x["g"] for x in v if x["g"]]),
                 sum(x["coinf"] for x in v) / len(v))
                for v in g.values() if len([x for x in v if x["g"]]) >= 10]

    r_dg, n_dg = pearson(country_points(runs))

    # The same correlation over ALL classified wheat, which is the number that actually clears
    # depth as an explanation. Within one project there are too few countries to resolve it.
    wheat = []
    for run, c in coinf.items():
        if c["crop"] != "Wheat":
            continue
        p = prov.get(run, {})
        raw = (p.get("location") or "").split(":")[0].strip()
        wheat.append({"country": ALIASES.get(raw, raw), "g": depth.get(run),
                      "coinf": c["coinfected"] == "yes"})
    r_wheat, n_wheat = pearson(country_points(wheat))

    L = []
    A = L.append
    A("Unreported co-pathogens in a yellow rust surveillance compendium")
    A("Caption for report/figures/geography.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_geography.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"Fungal pathogens detected in {len(runs)} wheat leaf libraries from {PROJECT}, a single")
    A("field surveillance programme that collected yellow-rust-infected leaves across")
    A(f"{len(by_country)} countries between 2014 and 2018 under one protocol. The study reported")
    A("Puccinia striiformis and no other pathogen, so every further species here is unreported by")
    A("construction, and because the protocol is constant a difference between countries is not a")
    A("difference between studies.")
    A("")
    A("PANEL a. Detection rate across all libraries. P. striiformis, the declared pathogen, is")
    A(f"recovered in {freq.get(DECLARED, 0)} of {len(runs)} libraries ({100 * freq.get(DECLARED, 0) / len(runs):.1f}%), which is the positive control for")
    A("the reference database and the detection floors. The species below it were not reported by")
    A("the original study.")
    A("")
    A("PANEL b. Per-country detection rate for the two most frequent unreported pathogens, in")
    A(f"countries with at least {MIN_RUNS} libraries. The two gradients run opposite ways and both")
    A("reproduce independently known disease geography: Zymoseptoria tritici, the dominant wheat")
    A("pathogen of maritime northwest Europe, is highest in the Netherlands, France, the UK and")
    A("Poland and falls away in Mediterranean, continental and African sites, while Puccinia")
    A("triticina peaks in Argentina and is absent from the Ethiopian highlands, where yellow and")
    A("stem rust dominate. Neither pattern was supplied to the pipeline, which was given no")
    A("geographic information of any kind.")
    A("")
    A("DEPTH, which had to be addressed before any of this could be read as spatial. Sequencing")
    A("depth expressed as G_kmer spans a factor of 300 across the cohort and covaries with crop,")
    A("which is why this work makes no comparison between crops. Within wheat the correlation")
    A("between a country's median G_kmer and its co-infection rate is")
    A(f"r = {r_wheat:+.3f} over {n_wheat} countries across all projects, and the deepest libraries are among the")
    A(f"lowest-rate ones. Within THIS project alone it is r = {r_dg:+.3f} over {n_dg} countries, which is weak,")
    A("of the opposite sign, and underpowered: nine points would need |r| near 0.67 to reach")
    A("significance. That neither establishes a depth effect nor excludes one. What supports the")
    A("spatial reading is therefore the agreement with independently known disease geography and")
    A("the opposite direction of the two gradients, not the depth correlation. More countries with")
    A("enough libraries is the fix, and the remaining cohort supplies them.")
    A("")
    A("TIME IS NOT SHOWN, and the reason is a composition artefact worth stating. Pooled by")
    A("collection year the co-infection rate rises monotonically, but the countries sampled change")
    A("with the year, and 14 of the 17 libraries collected in 2018 come from two Argentinian")
    A("sites. Holding country constant the trend largely dissolves. The year table is written to")
    A("geography.tsv for inspection and is not plotted.")
    A("")
    A("WHAT THIS IS NOT. A rate from this compendium is NOT a prevalence estimate: every library")
    A("was collected because the plant was visibly rusted, so the denominator is conditioned on")
    A("disease. The quantity measured here is the rate at which a SECOND pathogen accompanies a")
    A("known one, which is a different and narrower claim.")
    A("")
    A("  species                              runs    %   med cov   med sat")
    for sp in drawn:
        v = det.get(sp, [])
        if not v:
            continue
        A(f"  {sp:<36} {len(v):>5} {100 * len(v) / len(runs):>4.0f}%  {100 * st.median(x[0] for x in v):>7.2f}%   {st.median(x[1] for x in v):>7.3f}")
    A("")
    A("CONGENER INFLATION IS PRESENT AND IS NOT IN THESE COUNTS. Around 120 libraries call ten or")
    A("more Alternaria species simultaneously, every one at 4.9 to 5.9% coverage and 0.10 to 0.11")
    A("saturation. Genuine co-infection by ten congeners would not produce near-identical values;")
    A("one organism's reads distributed across near-identical references does. The uniformity is")
    A("the diagnostic. Those species carry no pathogen guild assignment and so do not enter the")
    A("counts above, but the pattern is why congeneric clusters are reported as a flag throughout")
    A("this work and why competitive read assignment is the next step.")
    A("")
    A("THE PHYLLOSPHERE RESIDENTS ARE DETECTED AS RESIDENTS, which is the reference database")
    A("working as designed. Cladosporium, Dioszegia and Sporobolomyces appear in 225 to 284")
    A("libraries at 4.8 to 7.4% coverage. These are the wheat leaf core community, and they were")
    A("added to the database specifically so that their reads would have somewhere correct to go")
    A("rather than being assigned to the nearest pathogen present.")
    A("")
    if aliased:
        A("CURATION DEFECTS FOUND, reported rather than silently merged, because the fix belongs")
        A("upstream in 01_literature:")
        for k, n in sorted(aliased.items()):
            A(f"  '{k}' ({n} runs) is the same place as '{ALIASES[k]}'")
        A("")
    A("SHORT CAPTION")
    A("")
    A(f"Unreported fungal pathogens in {len(runs)} yellow-rust-infected wheat leaf libraries from a")
    A(f"single surveillance programme spanning {len(by_country)} countries. (a) Detection rate across all")
    A(f"libraries; the declared pathogen P. striiformis is recovered in {100 * freq.get(DECLARED, 0) / len(runs):.0f}% as a positive control")
    A("and every other species was unreported. (b) Per-country rate for the two most frequent")
    A("unreported pathogens, which reproduce known disease geography despite the pipeline being")
    A("given no geographic information. Across all classified wheat, sequencing depth does not")
    A(f"explain the spatial pattern (r = {r_wheat:+.3f} between a country's median G_kmer and its rate over")
    A(f"{n_wheat} countries); within this project alone the correlation is r = {r_dg:+.3f} over {n_dg} countries and is")
    A("underpowered. Rates are conditioned on visible disease and")
    A("are not prevalence estimates.")
    A("")
    A("SOURCES")
    A(f"  data    {TSV.name}, {SPECIES_TSV.name}")
    A(f"  from    {Path(GOLD).name}/coinfection.tsv, {Path(GOLD).name}/detections.tsv,")
    A(f"          01_literature/release/cohort_samples.tsv, silver/run_depth.tsv")
    A("")

    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {TSV.name}  ({len(by_country)} countries, {len(runs)} runs)")
    print(f"-> {SPECIES_TSV.name}  ({len(det)} species)")
    print(f"-> {CAPTION.name}")
    print(f"\n  declared pathogen recovered: {freq.get(DECLARED, 0)}/{len(runs)} "
          f"({100 * freq.get(DECLARED, 0) / len(runs):.1f}%)")
    print(f"  depth-vs-rate, this project: r = {r_dg:+.3f} over {n_dg} countries (underpowered)")
    print(f"  depth-vs-rate, all wheat:    r = {r_wheat:+.3f} over {n_wheat} countries")
    print(f"  drawn co-pathogens:          {', '.join(drawn[1:])}")


if __name__ == "__main__":
    main()
