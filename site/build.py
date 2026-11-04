#!/usr/bin/env python3
"""build.py -- pack the release layer into the JSON the static site reads.

    python site/build.py

WHAT THIS IS. The site is static: no server, no database, no API. Every view is a projection
of files that already exist in `01_literature/release/` and `02_classify/data/gold/`, and this
script is the only thing that writes `site/data/`. That mirrors `build_release.py`: the site is
a projection of the seeds, so to change a number on the page you change the seed it comes from
and regenerate. NEVER hand-edit anything under site/data/.

WHY STATIC. The reader-facing data is small -- 505 sites, 142 called species, 2,862 panelled
runs, 191 evidence quotations -- and the whole payload is ~12 MB, which a browser holds without
complaint. A database buys nothing until submissions exist, and a submission tier is a separate
decision (see docs/web.md).

THE ONE SOURCE THAT IS NOT IN THE REPO. `report/figures/detect_runs/detect_runs_all_proposed.tsv`
is 195 MB and gitignored. It is the only table carrying `tier`, `trophic`, `evidence` and
`declared` together, which the run panel and the species view both need, so the site cannot be
rebuilt from a fresh clone alone. Regenerate it first with detect_run_single_proposed.R's
upstream step; the script says so and exits rather than silently emitting a thinner site.

FLOORS ARE COPIED, NOT REDERIVED. COV_FLOOR, SAT_FLOOR and MIN_DISTINCT must track
05_detect/detect.py. The site re-applies them client-side so a reader can move them, and a
moved floor is labelled WHAT-IF on the page. That is the same contract prep_browser.py set and
the same reason: a number read off a moved slider must not be mistakable for a pipeline result.

ROUNDING. Decision rows go to 6 significant figures, the rejected cloud to 3. At 4 a coverage
of 0.0199953 rounds to 0.02 and passes a 0.02 floor it should fail, so the rows a floor is
applied to are never rounded coarsely. Inherited from prep_browser.py, which found it.
"""

import collections
import csv
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # phd/01-review/
sys.path.insert(0, str(ROOT))

from _util import atomic_write           # noqa: E402

LIT = ROOT / "01_literature" / "release"
GOLD = ROOT / "02_classify" / "data" / "gold"
SILVER = ROOT / "02_classify" / "data" / "silver"
PANEL_SRC = ROOT / "02_classify" / "report" / "figures" / "detect_runs" / "detect_runs_all_proposed.tsv"
OUT = HERE / "data"

COV_FLOOR = 0.02        # 05_detect/detect.py
SAT_FLOOR = 0.05
MIN_DISTINCT = 10_000
CLOUD_N = 140           # rejected points kept per run, for the shape of the cloud
SEED = 20261030

TIERS = ["host", "pathogen_crop", "pathogen_other",
         "resident_crop", "resident_other", "no_record"]
TROPHIC = ["", "biotroph", "necrotroph", "weak", "opportunist"]
EVID = ["", "A", "B", "C"]


def tsv(path):
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def r6(v):
    try:
        return float(f"{float(v):.6g}")
    except (TypeError, ValueError):
        return None


def num(v, cast=float, default=None):
    try:
        return cast(float(v))
    except (TypeError, ValueError):
        return default


def write(name, payload):
    OUT.mkdir(parents=True, exist_ok=True)
    with atomic_write(OUT / name) as fh:
        json.dump(payload, fh, separators=(",", ":"), allow_nan=False)
    print(f"  {name:<16} {(OUT / name).stat().st_size / 1e6:6.2f} MB")


# ── the cohort, which every other view joins to ──────────────────────────────────────────────
def load_cohort():
    """One row per run: the curated host, tissue, setting, location and year, plus lat/lon.

    `location_precision` is carried through to the site and shown, because it is the difference
    between a paddock and a university postcode. `country` is split off the location string's
    leading token -- the archive's own convention, "Japan: Takatsuki" -- and is used only for
    grouping, never for placement; placement uses lat/lon or nothing.
    """
    cohort = {}
    for r in tsv(LIT / "cohort_samples.tsv"):
        loc = r["location"]
        cohort[r["Run"]] = {
            "p": r["BioProject"], "crop": r["crop"], "host": r["host"],
            "set": r["setting"], "sel": r["sampling_selection"], "tis": r["tissue"],
            "loc": loc, "country": loc.split(":")[0].strip() if loc else "",
            "year": (r["collection_date"] or "")[:4],
            "lat": num(r["lat"]), "lon": num(r["lon"]),
            "prec": r["location_precision"],
            "rules": {"setting": r["setting_rule"], "selection": r["sampling_selection_rule"],
                      "tissue": r["tissue_rule"], "location": r["location_rule"],
                      "year": r["collection_date_rule"]},
        }
    return cohort


# ── evidence: the quotation behind every curated claim ───────────────────────────────────────
def build_evidence():
    """`rule` -> the claim it asserts, its quotation, and where in which paper it was read.

    Keyed by rule id because that is what `cohort_samples.tsv` carries beside each curated
    value, so the site can answer "why does this run say field?" by a dictionary lookup rather
    than a search. One rule can back several fields, hence a list per key.
    """
    by_rule = collections.defaultdict(list)
    by_project = collections.defaultdict(list)
    for r in tsv(LIT / "evidence.tsv"):
        item = {"field": r["field"], "n": num(r["n_runs"], int, 0), "doi": r["doi"],
                "locus": r["locus"], "quote": r["quote"], "bp": r["bioproject"],
                "id": r["evidence_id"]}
        by_rule[r["rule"]].append(item)
        by_project[r["bioproject"]].append(item)

    projects = {}
    for r in tsv(LIT / "cohort_bioprojects.tsv"):
        projects[r.get("bioproject") or r.get("BioProject")] = r
    return {"by_rule": by_rule, "by_project": by_project, "projects": projects}


# ── the panel rows, per run ──────────────────────────────────────────────────────────────────
def build_runs(cohort):
    """Decision-relevant rows in full; the rejected cloud downsampled to a shape.

    A row is decision-relevant if it was called OR it cleared the absolute floor: those are the
    rows where the call was a decision, so the viewer can re-run the floors on them. Everything
    below the absolute floor carries no identity worth shipping and goes to the cloud as (x, y)
    only. Downsampling is stated on the page rather than hidden, which is the same rule
    prep_browser.py set.

    `congener` is carried because detect.py applies it AFTER the three floors -- a taxon whose
    only membership route is sharing a genus with a pathogen is not reportable -- and a re-call
    knowing only the floors called 8 runs differently. Found by diffing the viewer against the
    pipeline, not by reading the code.
    """
    rng = random.Random(SEED)
    runs, cur, dec, cloud, mass = {}, None, [], [], collections.Counter()

    def flush():
        if cur is None:
            return
        keep = cloud if len(cloud) <= CLOUD_N else rng.sample(cloud, CLOUD_N)
        c = cohort.get(cur["r"], {})
        runs[cur["r"]] = {
            **cur,
            "crop": c.get("crop", ""), "host": c.get("host", ""), "loc": c.get("loc", ""),
            "year": c.get("year", ""), "tis": c.get("tis", ""), "set": c.get("set", ""),
            "sel": c.get("sel", ""), "prec": c.get("prec", ""),
            "lat": c.get("lat"), "lon": c.get("lon"), "rules": c.get("rules", {}),
            # read mass by tier: what the library IS made of, before any identity question.
            "mass": [mass[t] for t in TIERS],
            "d": dec[:], "k": keep, "kn": len(cloud),
        }

    for row in tsv(PANEL_SRC):
        red, cov = num(row["redundancy"]), num(row["coverage"])
        if red is None or cov is None or red <= 0 or cov <= 0:
            continue
        if row["run"] != (cur or {}).get("r"):
            flush()
            dec, cloud, mass = [], [], collections.Counter()
            cur = {"r": row["run"], "p": row["BioProject"],
                   "dp": row.get("declared_pathogens", ""),
                   "np": num(row.get("n_pathogens"), int, 0),
                   "nu": num(row.get("n_undeclared"), int, 0)}
        tier = row["tier"] if row["tier"] in TIERS else "no_record"
        mass[tier] += num(row["reads"], int, 0)
        dist = num(row["distinct"], float, 0.0)
        if row["call"] == "yes" or dist >= MIN_DISTINCT:
            dec.append([
                row["species"],
                TIERS.index(tier),
                TROPHIC.index(row["trophic"]) if row.get("trophic") in TROPHIC else 0,
                EVID.index(row["evidence"]) if row.get("evidence") in EVID else 0,
                1 if row.get("declared") == "declared" else 0,
                r6(cov), r6(row["saturation"]), r6(red), int(dist),
                num(row["reads"], int, 0),
                1 if "congener" in (row.get("reason") or "") else 0,
            ])
        else:
            cloud.append([float(f"{red:.3g}"), float(f"{cov:.3g}")])
    flush()

    # Most pathogens first, so the list opens on the runs with something to look at.
    return sorted(runs.values(), key=lambda x: (-x["np"], x["r"]))


# ── the species view: the same detections, pivoted ───────────────────────────────────────────
def build_species(cohort):
    """For one species: where it was found, on what, declared or not, and how strongly.

    This is the pivot `detections.tsv` does not offer. The table is one row per (run, species),
    so the run view reads it down and the species view reads it across; nothing is recomputed,
    only regrouped.

    `declared` comes from the panel source, not from here -- detections.tsv has no such column
    and reconstructing it by matching species names against `declared_pathogens` strings is how
    a near-miss silently becomes a discovery.
    """
    declared = {}
    for row in tsv(PANEL_SRC):
        if row["call"] == "yes":
            declared[(row["run"], row["species"])] = row.get("declared") == "declared"

    clade = {}
    for row in tsv(GOLD / "detections.tsv"):
        clade[row["species"]] = num(row["db_clade"], int)

    agg = {}
    for row in tsv(GOLD / "detections.tsv"):
        sp = row["species"]
        s = agg.setdefault(sp, {
            "sp": sp, "taxid": row["taxid"], "kingdom": row["kingdom"],
            "guild": row["guild"], "route": row["route"], "pathogen": row["pathogen"] == "yes",
            "clade": clade.get(sp), "runs": [], "cov": [],
            "crop": collections.Counter(), "country": collections.Counter(),
            "proj": set(), "ndec": 0, "nundec": 0, "year": collections.Counter(),
        })
        c = cohort.get(row["run"], {})
        s["runs"].append(row["run"])
        s["proj"].add(row["BioProject"])
        s["crop"][c.get("crop", "unknown")] += 1
        if c.get("country"):
            s["country"][c["country"]] += 1
        if c.get("year"):
            s["year"][c["year"]] += 1
        cov = num(row["coverage"])
        if cov is not None:
            s["cov"].append(cov)
        if declared.get((row["run"], sp)):
            s["ndec"] += 1
        else:
            s["nundec"] += 1

    out = []
    for s in agg.values():
        cov = sorted(s["cov"])
        def q(p):
            return r6(cov[min(len(cov) - 1, int(p * len(cov)))]) if cov else None
        out.append({
            "sp": s["sp"], "taxid": s["taxid"], "kingdom": s["kingdom"],
            "guild": s["guild"], "route": s["route"], "pathogen": s["pathogen"],
            "clade": s["clade"], "n_runs": len(s["runs"]), "n_proj": len(s["proj"]),
            "n_declared": s["ndec"], "n_undeclared": s["nundec"],
            "cov_q": [q(0.0), q(0.25), q(0.5), q(0.75), q(0.99)],
            "crop": s["crop"].most_common(),
            "country": s["country"].most_common(),
            "year": sorted(s["year"].items()),
            # Capped: a species in 1,500 runs does not need all of them in the payload, and the
            # run view is one click away. Sorted so the cap is deterministic, not whichever
            # rows the reader happened to get.
            "runs": sorted(s["runs"])[:400],
        })
    return sorted(out, key=lambda x: -x["n_runs"])


# ── the map: runs aggregated to the place they were collected ────────────────────────────────
def build_sites(cohort, species):
    """One marker per distinct (location, lat, lon), carrying what was found there.

    Aggregated on the LOCATION STRING as well as the coordinate, because two sites can round to
    the same coordinate at 4 dp and still be different places in the record. The string is what
    the curation asserts; the coordinate is a placement of it.

    Runs with no coordinate are counted into an `unplaced` total and reported on the page. They
    are 149 of 7,688 and dropping them silently would make the map claim completeness it does
    not have.
    """
    by_run_sp = collections.defaultdict(list)
    for row in tsv(GOLD / "detections.tsv"):
        by_run_sp[row["run"]].append(row["species"])

    coinf = {r["run"]: r for r in tsv(GOLD / "coinfection.tsv")}

    sites, unplaced = {}, collections.Counter()
    for run, c in cohort.items():
        if c["lat"] is None or c["lon"] is None:
            unplaced[c["crop"]] += 1
            continue
        key = (c["loc"], round(c["lat"], 4), round(c["lon"], 4))
        s = sites.setdefault(key, {
            "loc": c["loc"], "lat": round(c["lat"], 4), "lon": round(c["lon"], 4),
            "country": c["country"], "prec": c["prec"],
            "n_runs": 0, "n_classified": 0, "proj": set(),
            "crop": collections.Counter(), "year": collections.Counter(),
            "sp": collections.Counter(), "n_coinf": 0, "n_undec": 0,
        })
        s["n_runs"] += 1
        s["proj"].add(c["p"])
        s["crop"][c["crop"]] += 1
        if c["year"]:
            s["year"][c["year"]] += 1
        if run in by_run_sp:
            s["n_classified"] += 1
            s["sp"].update(by_run_sp[run])
        cf = coinf.get(run)
        if cf:
            if cf.get("coinfected") in ("yes", "1", "True", "true"):
                s["n_coinf"] += 1
            s["n_undec"] += num(cf.get("n_undeclared"), int, 0)

    out = []
    for s in sites.values():
        out.append({**{k: v for k, v in s.items() if k not in ("proj", "crop", "year", "sp")},
                    "n_proj": len(s["proj"]),
                    "crop": s["crop"].most_common(),
                    "year": sorted(s["year"].items()),
                    "sp": s["sp"].most_common(12)})
    return sorted(out, key=lambda x: -x["n_runs"]), dict(unplaced)


def preview_page():
    """Emit `preview.html`: index.html with the document shell stripped.

    The Artifacts host wraps whatever it is given in its own `<!doctype><head><body>`, so a page
    that brings its own would nest one document inside another. index.html stays the real thing
    -- a complete file anyone can open from a Zenodo snapshot with no server -- and this is a
    projection of it, so the two cannot drift. GENERATED AND GITIGNORED; edit index.html.
    """
    src = (HERE / "index.html").read_text()
    body = src.split("<body", 1)[-1] if "<body" in src else src
    out, skip = [], ("<!doctype", "<meta charset", "<meta name=\"viewport\"")
    for line in body.splitlines():
        t = line.strip().lower()
        if any(t.startswith(k) for k in skip):
            continue
        out.append(line)
    (HERE / "preview.html").write_text("\n".join(out).lstrip("\n"))
    print(f"  preview.html     {(HERE / 'preview.html').stat().st_size / 1e6:6.2f} MB"
          "   (generated from index.html; do not edit)")


def main():
    if not PANEL_SRC.exists():
        sys.exit(f"missing {PANEL_SRC.relative_to(ROOT)}\n"
                 "It is 195 MB and gitignored. Regenerate it before building the site; a site\n"
                 "built without it would be missing tier, trophic, evidence and declared.")

    print("reading the release layer")
    cohort = load_cohort()
    ev = build_evidence()
    species = build_species(cohort)
    sites, unplaced = build_sites(cohort, species)
    runs = build_runs(cohort)

    funnel = [dict(r) for r in tsv(LIT / "funnel.tsv")]

    print("writing site/data/")
    write("meta.json", {
        "floors": {"cov": COV_FLOOR, "sat": SAT_FLOOR, "distinct": MIN_DISTINCT},
        "tiers": TIERS, "trophic": TROPHIC, "evid": EVID, "cloud_n": CLOUD_N,
        "counts": {"cohort_runs": len(cohort), "panelled_runs": len(runs),
                   "species": len(species), "sites": len(sites),
                   "projects": len(ev["projects"]), "evidence": sum(
                       len(v) for v in ev["by_project"].values())},
        "unplaced": unplaced, "funnel": funnel,
    })
    write("sites.json", sites)
    write("species.json", species)
    write("evidence.json", {"by_rule": ev["by_rule"], "by_project": ev["by_project"]})
    write("runs.json", {"runs": runs})
    preview_page()
    print(f"\n{len(runs)} runs, {len(species)} species, {len(sites)} sites, "
          f"{len(cohort)} cohort runs ({sum(unplaced.values())} unplaced)")


if __name__ == "__main__":
    main()
