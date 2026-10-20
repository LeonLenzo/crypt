#!/usr/bin/env python3
"""prep_pathogen_map.py — where the UNDECLARED pathogens turn up, in space.

One small-multiple map per secondary pathogen over all classified wheat. The question is
Leon's: not how often a second pathogen appears, but *where*.

## WHAT COUNTS AS A POINT, AND WHY IT IS NOT JUST "DETECTED"

A detection only speaks to this project's question if the study that produced the library did
not report it. Every wheat BioProject declares exactly one target, and the declared species
differs between projects, so the same species can be declared in one study and an undeclared
bystander in the next. A locality is therefore drawn for a species when at least one library
from that locality detected it **in a project that did not declare it**.

That distinction does real work. *Puccinia triticina* is the declared target of PRJEB84390 and
an undeclared bystander in the eight *Pst* projects; drawing both together would show the
sampling design rather than the biology.

## THE MAP LOCATES; IT DOES NOT MEASURE

Point density on this map follows SAMPLING EFFORT, because a point is a locality where a
species was found and more libraries make a find more likely. The per-country RATES in
geography.tsv have a denominator and this does not. Read the map for where a species turns up
and the rate table for how often, never the map for how often.

The trap is concrete and nearly caught me. *Microdochium nivale* looks purely northern on the
panel, which invites "snow mould, a cold-climate pathogen". But the median absolute latitude of
ALL sampled localities is 49.6 degrees and *M. nivale*'s is 51.6, so it is not distinguishable
from where people sampled. The script now writes that comparison so the eye is checked against
the baseline rather than trusted.

The two species that DO sit away from the baseline are *Cladosporium cladosporioides* (34.1)
and *Epicoccum nigrum* (36.1), both markedly lower-latitude, and both are phyllosphere
saprotrophs that FUNGuild also labels plant pathogens. Neither gets a panel.

## THE DENOMINATOR IS DRAWN

Every panel carries all {n} sampled localities as open grey points underneath. A map of
detections alone shows where people sampled, which for wheat is overwhelmingly northwest
Europe, and a reader cannot tell a real absence from an unsampled one. The grey layer is what
makes an empty region in a panel mean something.

## THE COORDINATES ARE NOT ALL EQUAL

88% of classified wheat libraries carry coordinates. Those at `country` or `region` precision
name a country and nothing finer, and are deliberately NOT given a centroid here: a centroid is
a statement about a country drawn at a place nobody sampled, and on a detection map it would
invent a field site. The house crop map draws them as diamonds because its question is "which
countries are represented"; this one is about location, so they are dropped and counted in the
caption instead.

Writes pathogen_map.tsv (one row per species x locality), pathogen_map_sites.tsv (the
denominator layer) and pathogen_map_caption.txt.

    python 02_classify/report/figures/prep_pathogen_map.py
"""

import csv
import sys
from math import comb
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import COHORT_SAMPLES, GOLD, LIT_RUNS  # noqa: E402
from _util import atomic_write  # noqa: E402

CROP = "Wheat"
MIN_SITES = 8          # a species needs this many distinct localities to get a panel
# Precisions that name a real place. 'country' and 'region' are excluded: see the docstring.
GOOD_PRECISION = {"locality", "locality-nofacility",
                  "corrected-high", "corrected-medium", "corrected-low"}

# The Europe window. 67% of localities and 69% of libraries sit inside it, which is why the
# world panels give the signal about a fifth of their area and a zoom is worth drawing.
EUROPE = {"lon": (-11.0, 32.0), "lat": (35.0, 60.0)}

# Maritime northwest Europe against the continental and southern remainder. The split is by
# climate, which is what wheat foliar pathogens actually respond to, and it is declared here as
# a hypothesis rather than drawn from the data after looking at it.
MARITIME = {"UK", "Ireland", "Netherlands", "Belgium", "Denmark"}

TSV = HERE / "pathogen_map.tsv"
EUROPE_TSV = HERE / "pathogen_map_europe.tsv"
SITES_TSV = HERE / "pathogen_map_sites.tsv"
CAPTION = HERE / "pathogen_map_caption.txt"

ALIASES = {"United Kingdom": "UK", "The Netherlands": "Netherlands", "Zimbabwae": "Zimbabwe"}


def species_of(name: str) -> str:
    parts = (name or "").replace("f. sp.", "").split()
    return " ".join(parts[:2]) if len(parts) >= 2 else (name or "").strip()


def main() -> None:
    coinf = {r["run"]: r for r in csv.DictReader(open(GOLD / "coinfection.tsv"), delimiter="\t")}
    prov = {r["Run"]: r for r in csv.DictReader(open(COHORT_SAMPLES), delimiter="\t")}
    lit = {r["Run"]: r for r in csv.DictReader(open(LIT_RUNS), delimiter="\t")}

    hits, guild = defaultdict(set), {}
    for r in csv.DictReader(open(GOLD / "detections.tsv"), delimiter="\t"):
        hits[r["run"]].add(r["species"])
        guild.setdefault(r["species"], r.get("guild", ""))

    runs, dropped = [], defaultdict(int)
    for run, c in coinf.items():
        if c["crop"] != CROP:
            continue
        p = prov.get(run, {})
        prec = (p.get("location_precision") or "").strip()
        lat, lon = (p.get("lat") or "").strip(), (p.get("lon") or "").strip()
        if not (lat and lon):
            dropped["no coordinates"] += 1
            continue
        if prec not in GOOD_PRECISION:
            dropped[f"precision '{prec or 'blank'}'"] += 1
            continue
        raw = (p.get("location") or "").split(":")[0].strip()
        runs.append({
            "run": run, "proj": c["BioProject"], "loc": (p.get("location") or "").strip(),
            "country": ALIASES.get(raw, raw), "lat": float(lat), "lon": float(lon),
            "hits": hits.get(run, set()),
            "declared": species_of(lit.get(run, {}).get("declared_pathogens") or ""),
        })

    # The denominator layer: every locality that was sampled at all.
    sites = defaultdict(lambda: {"n": 0, "lat": None, "lon": None, "country": "", "projs": set()})
    for r in runs:
        s = sites[r["loc"]]
        s["n"] += 1
        s["lat"], s["lon"], s["country"] = r["lat"], r["lon"], r["country"]
        s["projs"].add(r["proj"])
    with atomic_write(SITES_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["locality", "country", "lat", "lon", "n_libraries", "n_projects"])
        for loc, s in sorted(sites.items()):
            w.writerow([loc, s["country"], f"{s['lat']:.5f}", f"{s['lon']:.5f}",
                        s["n"], len(s["projs"])])

    # Candidate species: plant pathogens seen often enough to carry a map.
    freq = defaultdict(int)
    for r in runs:
        for sp in r["hits"]:
            if sp != r["declared"] and "Plant Pathogen" in guild.get(sp, ""):
                freq[sp] += 1

    rows = []
    for sp in sorted(freq, key=lambda s: -freq[s]):
        agg = defaultdict(lambda: {"n": 0, "det": 0})
        for r in runs:
            if r["declared"] == sp:          # declared here: not this project's question
                continue
            a = agg[r["loc"]]
            a["n"] += 1
            if sp in r["hits"]:
                a["det"] += 1
        n_sites = sum(1 for a in agg.values() if a["det"])
        if n_sites < MIN_SITES:
            continue
        for loc, a in agg.items():
            if not a["det"]:
                continue
            s = sites[loc]
            rows.append([sp, loc, s["country"], f"{s['lat']:.5f}", f"{s['lon']:.5f}",
                         a["n"], a["det"], f"{a['det'] / a['n']:.4f}"])

    drawn = []
    for sp in sorted(freq, key=lambda s: -freq[s]):
        n = sum(1 for r in rows if r[0] == sp)
        if n:
            drawn.append((sp, n, freq[sp]))

    with atomic_write(TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "locality", "country", "lat", "lon",
                    "n_libraries", "n_detected", "rate"])
        # Panel order is detection frequency, so the eye meets the strongest signal first.
        order = {sp: i for i, (sp, _, _) in enumerate(drawn)}
        for row in sorted(rows, key=lambda x: (order[x[0]], -int(x[6]))):
            w.writerow(row)

    # ── the Europe rate, which the map itself cannot show ────────────────────────────────────
    # A panel shows WHERE a species was found; this is HOW OFTEN, with a denominator, inside
    # the zoom window. Written so the figure can print it rather than implying it from point
    # density, which follows sampling effort.
    def in_europe(r):
        return (EUROPE["lon"][0] <= r["lon"] <= EUROPE["lon"][1]
                and EUROPE["lat"][0] <= r["lat"] <= EUROPE["lat"][1])

    eu_runs = [r for r in runs if in_europe(r)]
    with atomic_write(EUROPE_TSV) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "n_libraries", "n_detected", "rate",
                    "n_localities", "n_localities_detected"])
        for sp, _, _ in drawn:
            elig = [r for r in eu_runs if r["declared"] != sp]
            det_ = [r for r in elig if sp in r["hits"]]
            locs_all = {r["loc"] for r in elig}
            locs_det = {r["loc"] for r in det_}
            w.writerow([sp, len(elig), len(det_),
                        f"{len(det_) / len(elig):.4f}" if elig else "",
                        len(locs_all), len(locs_det)])

    # ── the maritime/continental test ────────────────────────────────────────────────────────
    # The EYE reads several gradients off the Europe panels. Most do not survive a denominator,
    # and this is the check that says which. Fisher's exact rather than chi-square because tan
    # spot's maritime cell is a single library.
    def fisher2(a, b, c, d):
        n, r1, r2, c1 = a + b + c + d, a + b, c + d, a + c
        def pr(x):
            return comb(r1, x) * comb(r2, c1 - x) / comb(n, c1)
        p0, tot = pr(a), 0.0
        for x in range(max(0, c1 - r2), min(r1, c1) + 1):
            px = pr(x)
            if px <= p0 * (1 + 1e-9):
                tot += px
        return min(1.0, tot)

    split = []
    for sp, _, _ in drawn:
        cell = {}
        for grp in ("maritime", "continental"):
            elig = [r for r in eu_runs
                    if r["declared"] != sp
                    and ((r["country"] in MARITIME) == (grp == "maritime"))]
            cell[grp] = (sum(1 for r in elig if sp in r["hits"]), len(elig))
        (ad, an), (cd, cn) = cell["maritime"], cell["continental"]
        if not (an and cn):
            continue
        ra, rc = ad / an, cd / cn
        split.append((sp, ad, an, ra, cd, cn, rc,
                      (ra / rc) if rc else float("inf"),
                      fisher2(ad, an - ad, cd, cn - cd)))
    bonf = 0.05 / len(split) if split else float("nan")
    with atomic_write(HERE / "pathogen_map_climate.tsv") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "maritime_detected", "maritime_n", "maritime_rate",
                    "continental_detected", "continental_n", "continental_rate",
                    "ratio", "fisher_p", "survives_bonferroni"])
        for sp, ad, an, ra, cd, cn, rc, rt, pv in split:
            w.writerow([sp, ad, an, f"{ra:.4f}", cd, cn, f"{rc:.4f}",
                        f"{rt:.3f}", f"{pv:.3e}", "yes" if pv < bonf else "no"])

    lons = [s["lon"] for s in sites.values()]
    lats = [s["lat"] for s in sites.values()]

    # The latitude baseline. A species' detection latitudes mean nothing on their own; what
    # matters is whether they differ from where anybody sampled.
    import statistics as _st
    base_lat = _st.median(abs(x) for x in lats)
    lat_rows = []
    for sp, n_sites, _ in drawn:
        ls = [float(r[3]) for r in rows if r[0] == sp]
        lat_rows.append((sp, n_sites, _st.median(abs(x) for x in ls)))

    L = []
    A = L.append
    A("Where the undeclared pathogens are")
    A("Caption for report/figures/pathogen_map.{png,pdf}")
    A("")
    A("GENERATED by report/figures/prep_pathogen_map.py -- do not hand-edit.")
    A("")
    A("")
    A("FULL CAPTION")
    A("")
    A(f"Geographic distribution of fungal pathogens detected in {len(runs):,} classified wheat libraries")
    A(f"that carry field coordinates, across {len(sites)} localities in {len({s['country'] for s in sites.values()})} countries. Each panel is one")
    A("species, and a locality is marked only where that species was detected in a library whose")
    A("study DID NOT declare it, so every coloured point is an unreported detection.")
    A("")
    A("That restriction matters rather than being a formality. Puccinia triticina is the declared")
    A("target of one project and an undeclared bystander in eight others; mapping both together")
    A("would draw the sampling design instead of the biology.")
    A("")
    A("OPEN GREY POINTS ARE EVERY SAMPLED LOCALITY and are repeated in all panels. A detection")
    A("map without its denominator shows where people sampled, which for wheat is overwhelmingly")
    A("northwest Europe, and an empty region could mean nobody looked. The grey layer is what")
    A("makes an absence readable.")
    A("")
    A("  species                          localities    libraries")
    for sp, n_sites, n_lib in drawn:
        A(f"  {sp:<34} {n_sites:>7}      {n_lib:>7}")
    A("")
    A("THE MAP LOCATES, IT DOES NOT MEASURE. Point density follows sampling effort: a point is a")
    A("locality where a species was found, and more libraries make a find more likely. For how")
    A("OFTEN a species occurs, which has a denominator, see the per-country rates in")
    A("geography.tsv. The latitude of detections is the clearest illustration, since the median")
    A(f"absolute latitude of every sampled locality is {base_lat:.1f} degrees:")
    A("")
    A("  species                          localities   median |lat|   vs sampling baseline")
    for sp, n_sites, ml in lat_rows:
        d = ml - base_lat
        note = "indistinguishable" if abs(d) < 4 else ("LOWER latitude" if d < 0 else "higher latitude")
        A(f"  {sp:<34} {n_sites:>7}   {ml:>10.1f}     {d:>+5.1f}  {note}")
    A("")
    A("Microdochium nivale looks purely northern on its panel, which invites reading snow mould")
    A("as a cold-climate signal. It is not distinguishable from the sampling baseline. The two")
    A("species that DO sit away from it are Cladosporium cladosporioides and Epicoccum nigrum,")
    A("both lower-latitude and both phyllosphere saprotrophs that FUNGuild also labels plant")
    A("pathogens; neither is drawn.")
    A("")
    A(f"THE EUROPE PANEL SET covers {EUROPE['lon'][0]:.0f} to {EUROPE['lon'][1]:.0f} longitude and {EUROPE['lat'][0]:.0f} to {EUROPE['lat'][1]:.0f} latitude, which holds")
    A(f"{len(eu_runs)} of the {len(runs)} drawn libraries. Each panel prints its European detection rate, which")
    A("is the quantity with a denominator; the points show where, the printed rate shows how")
    A("often:")
    A("")
    A("  species                          libraries   detected      rate")
    for row in csv.DictReader(open(EUROPE_TSV), delimiter="\t"):
        A(f"  {row['species']:<34} {int(row['n_libraries']):>7}   {int(row['n_detected']):>8}   {float(row['rate']):>7.1%}")
    A("")
    A("ONE OF THE SIX APPARENT GRADIENTS SURVIVES A DENOMINATOR. Maritime northwest Europe (UK,")
    A("Ireland, Netherlands, Belgium, Denmark) against the continental and southern remainder,")
    A(f"Fisher's exact, Bonferroni over {len(split)} tests so significance needs p < {bonf:.4f}:")
    A("")
    A("  species                          maritime   continental   ratio    Fisher p")
    for sp, ad, an, ra, cd, cn, rc, rt, pv in sorted(split, key=lambda x: x[8]):
        mark = "  <-- survives" if pv < bonf else ""
        A(f"  {sp:<34} {ra:>6.1%}      {rc:>6.1%}   {rt:>5.2f}   {pv:>9.2e}{mark}")
    A("")
    A("BOTH SURVIVORS ARE Pyrenophora, AND THEY ARE NOT EACH OTHER. A genus-level coincidence")
    A("in this pipeline is first a congener-smearing suspect, so it was checked: of the 60")
    A("libraries calling either species, only 3 call both, 40 are tan spot alone and 17 net")
    A("blotch alone. Smearing produces near-complete overlap at near-identical coverage, which is")
    A("what the Alternaria cluster does; these are independent, and where each is detected alone")
    A("the median coverages differ (7.0% against 2.7%). Two stubble-borne Pleosporales sharing a")
    A("continental distribution is coherent rather than suspicious. P. teres is primarily a")
    A("BARLEY pathogen and its presence on wheat, while recorded, is worth flagging.")
    A("")
    A("Pyrenophora tritici-repentis is detected in 1 of 350 maritime libraries and 15 of 292")
    A("continental ones, an eighteenfold difference, and tan spot is indeed a disease of")
    A("reduced-tillage continental wheat rather than of the maritime northwest. Everything else")
    A("fails, INCLUDING the two gradients the panels make most obvious: Zymoseptoria tritici")
    A("reads as strongly maritime and is 42.3% against 34.3% at p = 0.04, which does not survive")
    A("correction, and Microdochium nivale looks almost exclusively British and is 9.1% against")
    A("7.2% at p = 0.39. Both are")
    A("mostly the sampling: 354 of the 648 European libraries are maritime, and 317 of those are")
    A("from the UK alone.")
    A("")
    A("POINT SIZE is the number of libraries at that locality in which the species was detected,")
    A("not a rate. Most localities contributed a single library, so a rate there would be 0 or")
    A("100% and would overstate what one sample can say.")
    A("")
    A(f"COORDINATE HANDLING. {sum(dropped.values())} of {sum(dropped.values()) + len(runs):,} classified wheat libraries are not drawn:")
    for k, n in sorted(dropped.items(), key=lambda x: -x[1]):
        A(f"  {n:>4}  {k}")
    A("Libraries whose location resolves only to a country or region are excluded rather than")
    A("placed at a centroid. On a map about where disease was found, a centroid invents a field")
    A("site at a place nobody sampled.")
    A("")
    A("WHAT THIS IS NOT. Every library was collected because the plant was diseased, so these")
    A("are not prevalence surfaces. A locality with no coloured point was sampled and the species")
    A("was not detected there, which under a disease-conditioned design is weaker than absence.")
    A("")
    A("SHORT CAPTION")
    A("")
    A(f"Undeclared fungal pathogens across {len(sites)} wheat sampling localities in {len({s['country'] for s in sites.values()})} countries.")
    A("Each panel is one species; coloured points mark localities where it was detected in a")
    A("library whose study did not report it, sized by the number of such libraries. Open grey")
    A("points are every sampled locality and are repeated in each panel, so an empty region can")
    A("be read as sampled-and-not-found rather than unsampled. Detections are conditioned on")
    A("visible disease and are not prevalence estimates.")
    A("")
    A("SOURCES")
    A(f"  data  {TSV.name}, {SITES_TSV.name}")
    A("")

    with atomic_write(CAPTION) as fh:
        fh.write("\n".join(L))

    print(f"-> {TSV.name}  ({len(rows)} species x locality rows)")
    print(f"-> {EUROPE_TSV.name}, pathogen_map_climate.tsv")
    print(f"-> {SITES_TSV.name}  ({len(sites)} localities)")
    print(f"-> {CAPTION.name}")
    print(f"\n  drawn libraries {len(runs):,}, localities {len(sites)}, "
          f"lon {min(lons):.0f} to {max(lons):.0f}, lat {min(lats):.0f} to {max(lats):.0f}")
    print(f"  excluded: {dict(dropped)}")
    print(f"  latitude baseline (all sampled localities): {base_lat:.1f} deg")
    print(f"  maritime vs continental, Bonferroni p < {bonf:.4f}:")
    for sp, ad, an, ra, cd, cn, rc, rt, pv in sorted(split, key=lambda x: x[8])[:4]:
        print(f"    {sp:<34} {ra:>6.1%} vs {rc:>6.1%}  ratio {rt:>5.2f}  p={pv:.2e}"
              f"{'  SURVIVES' if pv < bonf else ''}")
    for sp, n_sites, n_lib in drawn:
        print(f"    {sp:<34} {n_sites:>4} localities, {n_lib:>4} libraries")


if __name__ == "__main__":
    main()
