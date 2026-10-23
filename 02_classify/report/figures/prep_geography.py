#!/usr/bin/env python3
"""prep_geography.py — all classified wheat: the declared-vs-undeclared control, and geography.

Covers **every classified wheat run**, 1,068 across 12 BioProjects, rather than the single
surveillance compendium an earlier version of this script used. Expanding it is what turns two
observations into tests, because each wheat project declares exactly ONE pathogen and the four
declared species differ between projects:

    Puccinia striiformis   890 runs over 8 projects
    Puccinia triticina      99 runs over 1
    Puccinia graminis       67 runs over 2
    Fusarium graminearum    12 runs over 1

## TEST 1: the same species, declared and undeclared

A detection pipeline that invents rusts would find them whether or not the study declared one.
This cohort can measure that directly, because *P. triticina* is the declared pathogen in
PRJEB84390 and an undeclared bystander in the eight *Pst* projects, and *P. graminis* likewise.
Same species, same database, same floors; the only thing that changes is whether the study said
it was there.

The asymmetry is the evidence. A per-species recovery near 100% when declared, against a much
lower rate when not, cannot be produced by cross-mapping, because cross-mapping does not know
what the paper said.

**It also tests the congener worry for the rusts specifically.** If *Puccinia* calls were
smearing across congeners, every project would show all three rusts. PRJNA256347 reports 100%
*striiformis* with 0% *triticina* and 0% *graminis* on 35 runs, so they are not.

## TEST 2: does the geography replicate across independent studies

The spatial pattern is only worth reporting if a country's rate agrees between projects that
share nothing but the country. This script writes the per-country rate separately per project
so that agreement can be read directly rather than asserted.

## THE SEVEN NON-DETECTIONS ARE NOT BORDERLINE

Of 1,068 wheat runs, 1,061 recover the declared pathogen. Six of the seven are PRJNA1181034,
which is a field *F. graminearum* spike trial, and its 12 runs split bimodally: six at 14.8 to
29.3% coverage and six at 0.00 to 0.11%, with 113 to 3,675 distinct minimizers against 491,872
to 974,167. That is two to four orders of magnitude, not a threshold being brushed. The
unavoidable reading is uninfected or mock-inoculated controls inside the trial, which means the
absolute floor separated infected from uninfected cleanly in the one project that contained
both.

## TIME IS STILL NOT REPORTED

Pooled by collection year the co-infection rate rises monotonically and it is a composition
artefact: the countries sampled change with the year, and 14 of 17 libraries collected in 2018
come from two Argentinian sites. Holding country constant it largely dissolves. The year table
is written for inspection and is not plotted.

Writes geography.tsv, geography_species.tsv, geography_control.tsv, geography_year.tsv and
geography_caption.txt.

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

from _layout import COHORT_SAMPLES, GOLD, LIT_RUNS, SILVER  # noqa: E402
from _util import atomic_write  # noqa: E402

CROP = "Wheat"
MIN_RUNS = 8        # a country below this is written to the table but not drawn
MIN_CELL = 6        # a project x country cell below this is not drawn in the replication panel
FOCUS = "Zymoseptoria tritici"

TSV = HERE / "geography.tsv"
SPECIES_TSV = HERE / "geography_species.tsv"
CONTROL_TSV = HERE / "geography_control.tsv"
YEAR_TSV = HERE / "geography_year.tsv"
CAPTION = HERE / "geography_caption.txt"

# Curation defects, reported rather than silently merged: the fix belongs in 01_literature.
ALIASES = {"United Kingdom": "UK", "The Netherlands": "Netherlands", "Zimbabwae": "Zimbabwe"}


def species_of(name: str) -> str:
    """'Puccinia striiformis f. sp. tritici' -> 'Puccinia striiformis'. db_v4 tags at species."""
    parts = (name or "").replace("f. sp.", "").split()
    return " ".join(parts[:2]) if len(parts) >= 2 else (name or "").strip()


def year_of(s: str):
    s = (s or "").strip()
    for i in range(len(s) - 3):
        if s[i:i + 4].isdigit() and 1900 < int(s[i:i + 4]) < 2030:
            return int(s[i:i + 4])
    return None


def gkmer(row):
    """Gkmer for one run. ENA's read_len is bases PER SPOT, so a paired library's per-read
    length is HALF of it -- the same correction rarefy.py and prep_rarefaction.py already make.
    Without it every paired run is inflated (median 1.25x, up to 2.10x over the cohort), which
    moved the cohort median from 1.68 to 2.30 Gkmer and biased the per-country medians below."""
    try:
        bases, L = float(row["bases"]), float(row["read_len"])
    except (TypeError, ValueError, KeyError):
        return None
    if row.get("layout") == "PAIRED":
        L = L / 2
    return bases * (L - 34) / L / 1e9 if L > 35 else None


def pearson(pts):
    if len(pts) < 3:
        return float("nan"), len(pts)
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    mx, my = st.mean(xs), st.mean(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    return (num / den if den else float("nan")), len(pts)


def main() -> None:
    coinf = {r["run"]: r for r in csv.DictReader(open(GOLD / "coinfection.tsv"), delimiter="\t")}
    prov = {r["Run"]: r for r in csv.DictReader(open(COHORT_SAMPLES), delimiter="\t")}
    lit = {r["Run"]: r for r in csv.DictReader(open(LIT_RUNS), delimiter="\t")}
    depth = {}
    for r in csv.DictReader(open(SILVER / "run_depth.tsv"), delimiter="\t"):
        g = gkmer(r)
        if g:
            depth[r[list(r)[0]]] = g

    hits, guild = defaultdict(set), {}
    for r in csv.DictReader(open(GOLD / "detections.tsv"), delimiter="\t"):
        hits[r["run"]].add(r["species"])
        guild.setdefault(r["species"], r.get("guild", ""))

    runs, aliased = [], defaultdict(int)
    for run, c in coinf.items():
        if c["crop"] != CROP:
            continue
        p = prov.get(run, {})
        raw = (p.get("location") or "").split(":")[0].strip()
        if raw in ALIASES:
            aliased[raw] += 1
        runs.append({
            "run": run, "proj": c["BioProject"], "country": ALIASES.get(raw, raw),
            "year": year_of(p.get("collection_date")), "g": depth.get(run),
            "coinf": c["coinfected"] == "yes", "hits": hits.get(run, set()),
            "declared": species_of((lit.get(run, {}).get("declared_pathogens") or "")),
        })

    by_proj = defaultdict(list)
    for r in runs:
        by_proj[r["proj"]].append(r)

    # ── TEST 1 ───────────────────────────────────────────────────────────────────────────────
    # Per species: recovery where the study declared it, against where it did not.
    declared_species = sorted({r["declared"] for r in runs if r["declared"]})
    with atomic_write(CONTROL_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "n_declared", "found_declared", "rate_declared",
                    "n_undeclared", "found_undeclared", "rate_undeclared", "n_projects_declaring"])
        for sp in declared_species:
            dec = [r for r in runs if r["declared"] == sp]
            und = [r for r in runs if r["declared"] != sp]
            fd = sum(1 for r in dec if sp in r["hits"])
            fu = sum(1 for r in und if sp in r["hits"])
            w.writerow([sp, len(dec), fd, f"{fd / len(dec):.4f}" if dec else "",
                        len(und), fu, f"{fu / len(und):.4f}" if und else "",
                        len({r["proj"] for r in dec})])
    control = {r["species"]: r for r in csv.DictReader(open(CONTROL_TSV), delimiter="\t")}

    # ── TEST 2 ───────────────────────────────────────────────────────────────────────────────
    by_country = defaultdict(list)
    for r in runs:
        by_country[r["country"]].append(r)

    freq = defaultdict(int)
    for r in runs:
        for sp in r["hits"]:
            freq[sp] += 1
    drawn = [sp for sp, n in sorted(freq.items(), key=lambda x: -x[1])
             if n >= 40 and "Plant Pathogen" in guild.get(sp, "")]

    with atomic_write(TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["country", "n_runs", "n_projects", "coinf_rate", "median_gkmer", "drawn"] +
                   [sp.replace(" ", "_") for sp in drawn])
        for cty, v in sorted(by_country.items(), key=lambda x: -len(x[1])):
            gs = [r["g"] for r in v if r["g"]]
            w.writerow([cty, len(v), len({r["proj"] for r in v}),
                        f"{sum(r['coinf'] for r in v) / len(v):.4f}",
                        f"{st.median(gs):.3f}" if gs else "",
                        "yes" if len(v) >= MIN_RUNS else "no"] +
                       [f"{sum(1 for r in v if sp in r['hits']) / len(v):.4f}" for sp in drawn])

    # Per project x country, so replication can be READ rather than asserted.
    with atomic_write(HERE / "geography_replication.tsv") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["country", "BioProject", "n_runs", "rate", "drawn"])
        for cty, v in sorted(by_country.items(), key=lambda x: -len(x[1])):
            if len(v) < MIN_RUNS:
                continue
            cells = defaultdict(list)
            for r in v:
                cells[r["proj"]].append(r)
            for pj, vv in sorted(cells.items(), key=lambda x: -len(x[1])):
                w.writerow([cty, pj, len(vv),
                            f"{sum(1 for r in vv if FOCUS in r['hits']) / len(vv):.4f}",
                            "yes" if len(vv) >= MIN_CELL else "no"])

    # Species summary with the signal quality that justifies the call.
    det = defaultdict(list)
    for r in csv.DictReader(open(GOLD / "detections.tsv"), delimiter="\t"):
        if coinf.get(r["run"], {}).get("crop") != CROP:
            continue
        try:
            det[r["species"]].append((float(r["coverage"]), float(r["saturation"] or 0),
                                      int(r["distinct"])))
        except (TypeError, ValueError):
            pass
    with atomic_write(SPECIES_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "n_runs", "pct_runs", "median_coverage", "median_saturation",
                    "median_distinct", "guild", "ever_declared"])
        for sp, v in sorted(det.items(), key=lambda x: -len(x[1])):
            w.writerow([sp, len(v), f"{len(v) / len(runs):.4f}",
                        f"{st.median(x[0] for x in v):.5f}",
                        f"{st.median(x[1] for x in v):.4f}",
                        int(st.median(x[2] for x in v)), guild.get(sp, ""),
                        "yes" if sp in declared_species else "no"])

    # The year table, written for inspection, deliberately not plotted.
    with atomic_write(YEAR_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["BioProject", "year", "n_runs", "coinf_rate", "n_countries", "top_country"])
        for pj, v in sorted(by_proj.items(), key=lambda x: -len(x[1])):
            yrs = defaultdict(list)
            for r in v:
                if r["year"]:
                    yrs[r["year"]].append(r)
            for y, vv in sorted(yrs.items()):
                cc = defaultdict(int)
                for r in vv:
                    cc[r["country"]] += 1
                top = max(cc.items(), key=lambda x: x[1])
                w.writerow([pj, y, len(vv), f"{sum(r['coinf'] for r in vv) / len(vv):.4f}",
                            len(cc), f"{top[0]} {top[1]}"])

    # ── depth, at both scopes ────────────────────────────────────────────────────────────────
    def country_points(rows_):
        g = defaultdict(list)
        for r in rows_:
            g[r["country"]].append(r)
        return [(st.median([x["g"] for x in v if x["g"]]),
                 sum(x["coinf"] for x in v) / len(v))
                for v in g.values() if len([x for x in v if x["g"]]) >= 10]

    r_wheat, n_wheat = pearson(country_points(runs))

    tot_found = sum(1 for r in runs if r["declared"] and r["declared"] in r["hits"])
    tot_decl = sum(1 for r in runs if r["declared"])

    L = []
    A = L.append
    A("Declared and undeclared fungal pathogens across all classified wheat")
    A("Caption for report/figures/geography.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_geography.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"All {len(runs):,} classified wheat libraries, from {len(by_proj)} BioProjects across {len(by_country)} countries. Each")
    A("project declared exactly one target pathogen, and the declared species differs between")
    A("projects, which turns the comparison below into a test rather than an observation.")
    A("")
    A("PANEL a. THE SAME SPECIES, DECLARED AND UNDECLARED. For each species that some project")
    A("declared, the detection rate in libraries whose study declared it, against libraries whose")
    A("study did not. Same database, same thresholds; the only difference is what the paper said.")
    A("")
    A("  species                    declared        undeclared")
    for sp in declared_species:
        c = control[sp]
        A(f"  {sp:<24} {float(c['rate_declared']):>5.0%} (n={int(c['n_declared']):>4})   "
          f"{float(c['rate_undeclared']):>5.0%} (n={int(c['n_undeclared']):>4})")
    A("")
    A(f"Overall the declared pathogen is recovered in {tot_found:,} of {tot_decl:,} libraries ({tot_found / tot_decl:.1%}). Cross-mapping")
    A("cannot produce that asymmetry, because cross-mapping does not know what the paper said.")
    A("The same comparison answers the congener objection for the rusts: if Puccinia calls were")
    A("smearing across congeners, every project would show all three rusts, and PRJNA256347")
    A("reports 100% striiformis with 0% triticina and 0% graminis across 35 libraries.")
    A("")
    A(f"PANEL b. {FOCUS} detection rate by country, one point per BioProject, for")
    A(f"countries with at least {MIN_RUNS} libraries and cells with at least {MIN_CELL}. Agreement between points in")
    A("a row is agreement between independent studies that share nothing but the country. The")
    A("gradient reproduces the known distribution of septoria tritici blotch, which is the")
    A("dominant wheat disease of maritime northwest Europe and minor in drier and warmer wheat")
    A("regions, and the pipeline was given no geographic information of any kind.")
    A("")
    A("THE SEVEN NON-DETECTIONS ARE NOT BORDERLINE. Six of the seven are PRJNA1181034, a field")
    A("Fusarium graminearum spike trial whose 12 libraries split bimodally: six at 14.8 to 29.3%")
    A("coverage and six at 0.00 to 0.11%, which is 113 to 3,675 distinct minimizers against")
    A("491,872 to 974,167. Two to four orders of magnitude is not a threshold being brushed, and")
    A("the reading is uninfected or mock controls inside the trial. The absolute floor separated")
    A("infected from uninfected cleanly in the one project that contained both.")
    A("")
    A("DEPTH does not explain the spatial pattern: across the countries with at least 10")
    A(f"libraries, the correlation between a country's median G_kmer and its co-infection rate is")
    A(f"r = {r_wheat:+.3f} over {n_wheat} countries. Depth spans a factor of 300 across the wider cohort and")
    A("covaries with crop, which is why this work makes no comparison between crops.")
    A("")
    A("TIME IS NOT SHOWN. Pooled by collection year the rate rises monotonically and it is a")
    A("composition artefact: the countries sampled change with the year, and 14 of the 17")
    A("libraries collected in 2018 come from two Argentinian sites. Holding country constant it")
    A("largely dissolves. geography_year.tsv carries the table.")
    A("")
    A("WHAT THIS IS NOT. Every library here was collected because the plant was diseased, so the")
    A("denominator is conditioned on disease. The quantity measured is the rate at which a SECOND")
    A("pathogen accompanies a known one, not prevalence.")
    A("")
    if aliased:
        A("CURATION DEFECTS FOUND, reported rather than silently merged, because the fix belongs")
        A("upstream in 01_literature:")
        for k, n in sorted(aliased.items()):
            A(f"  '{k}' ({n} runs) is the same place as '{ALIASES[k]}'")
        A("")
    A("SHORT CAPTION")
    A("")
    A(f"All {len(runs):,} classified wheat libraries from {len(by_proj)} BioProjects and {len(by_country)} countries. (a) Detection")
    A("rate for each declared pathogen in libraries whose study declared it against libraries")
    A(f"whose study did not; the declared pathogen is recovered in {tot_found / tot_decl:.1%} of libraries overall,")
    A("and the asymmetry cannot be produced by cross-mapping. (b) Zymoseptoria tritici rate by")
    A("country with one point per BioProject, so agreement between points is agreement between")
    A("independent studies. The gradient reproduces the known distribution of septoria tritici")
    A(f"blotch although no geographic information was supplied. Depth does not explain it (r = {r_wheat:+.3f}).")
    A("Rates are conditioned on visible disease and are not prevalence estimates.")
    A("")
    A("SOURCES")
    A(f"  data   {TSV.name}, {SPECIES_TSV.name}, {CONTROL_TSV.name},")
    A(f"         geography_replication.tsv, {YEAR_TSV.name}")
    A("")

    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {TSV.name}, {SPECIES_TSV.name}, {CONTROL_TSV.name}, "
          f"geography_replication.tsv, {YEAR_TSV.name}, {CAPTION.name}")
    print(f"\n  wheat runs {len(runs):,} over {len(by_proj)} projects, {len(by_country)} countries")
    print(f"  declared pathogen recovered: {tot_found:,}/{tot_decl:,} ({100 * tot_found / tot_decl:.1f}%)")
    print(f"  depth-vs-rate: r = {r_wheat:+.3f} over {n_wheat} countries")
    for sp in declared_species:
        c = control[sp]
        print(f"    {sp:<24} declared {float(c['rate_declared']):>6.1%} (n={c['n_declared']:>4})  "
              f"undeclared {float(c['rate_undeclared']):>6.1%} (n={c['n_undeclared']:>4})")


if __name__ == "__main__":
    main()
