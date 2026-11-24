#!/usr/bin/env python3
"""Decide which Kraken2 assignments in one arm are real organisms.

    python release/03_classify/scripts/detect.py --arm db1p
    python release/03_classify/scripts/detect.py --arm db2 --cov-floor 0.03
    python release/03_classify/scripts/detect.py --arm db1p --sensitivity

Reads  data/raw/reports/<arm>/<BioProject>/<Run>.report.tsv   the classification
       data/inputs/thresholds.tsv                             the operating point
       ../02_reference/data/outputs/db_inspect_<arm>.tsv       the coverage denominator
       ../02_reference/data/outputs/supp_data_s3_reference.tsv the labels, one row per taxon
       ../02_reference/data/curation/extra_cds_db2.tsv         which taxa are predicted CDS
       ../01_cohort/data/outputs/supp_data_s2_cohort.tsv       host, library, declaration
Writes data/raw/candidates_<arm>.tsv    every species row with its metrics and its verdict
       data/outputs/detections_<arm>.tsv  the calls

## The call

Four tests, all of which must hold, in this order, every one written POSITIVELY:

    distinct   >= distinct_floor      is there enough signal to divide at all
    coverage   >= coverage_floor      is the DNA there
    saturation >= saturation_floor    is it abundant RELATIVE TO WHAT THIS DEPTH ALLOWED
    fsp_frac   <= max_fsp             is one forma specialis carrying the call

## saturation is a depth-normalised ABUNDANCE floor, not a shape test

It was specified as a shape test -- are the reads spread over the taxon's k-mer space or
piled on a few fragments -- and it cannot be one, which was measured on db2's 430,820
candidates rather than argued:

  * Every row saturation alone rejects has redundancy >= 0.534. Every row the coverage
    floor alone rejects has redundancy <= 0.480. The two partition at r = 0.51 with no
    overlap, because that is where 0.05*(1-exp(-r)) crosses 0.02.
  * Above r = 3 the median saturation/coverage ratio is 1.007, and above 10 it is 1.000.
    1-exp(-r) has gone to 1 and the statistic IS coverage.

So saturation only ever acts where the shape information has vanished, and in the regime
where it would carry shape information the coverage floor has already bound. Taken
together the two are ONE floor, coverage >= max(coverage_floor, saturation_floor*(1-e^-r)).

That is still worth having, and it is not what 04_align does. Read correctly it asks for
5% of what the depth ALLOWED: shallow libraries get a proportionally lower bar because the
taxon could not have reached 5% whatever was there, and at full depth it becomes a flat 5%.
04_align tests the genuine shape question -- where in a transcript the reads fell -- which
needs positions this module does not have. The two do not overlap.

A missing value must FAIL the test. `coverage` is None whenever the report has no minimizer
columns or the taxon is absent from the inspect table, and None compares False against
everything, so `if coverage < floor: continue` ADMITS those rows. That inflated a cereal
species count from 32 to 44 before it was caught, on a table where 13.8% of rows carried NaN.
Every test here reads `value is not None and value >= floor`.

The order matters. The absolute floor is first because a relative floor alone is meaningless at
the bottom of the ceiling distribution. Below `min_redundancy` saturation is silent and the
call rests on the absolute and coverage floors, which is recorded per row in `judgeable`.

## Three rules about labels, and all three are load-bearing

**Pathogen status is a label applied AFTER detection, never a filter applied before it.**
`on_wheat`, `is_plant_pathogen` and `is_wheat_pathogen` are columns on S3 and are read off a
call. The study exists to find unreported host associations, so filtering the candidates to
what is already known on wheat would make the finding impossible to make.

**Predicted CDS is a competitor only and can never be a positive call.** The 135 gene-predicted
genomes are in db2 to absorb reads that would otherwise be forced onto a congener. A detection
resting on them would be an artefact of our own build, so they are called `competitor_only`
whatever their metrics, and their metrics are still written to the candidate table because what
they absorbed is the purpose test's measurement. Refusing to run without that list is
deliberate: a db2 pass that cannot identify the predicted taxa would silently report them.

**An uncurated declaration is UNKNOWN, never undeclared.** 5 of the 31 cohort BioProjects
(88 runs) still have no curated declaration, and counting those as "nobody reported it" would
manufacture the chapter's headline. `declared` is one of `declared`, `undeclared` or `unknown`,
and the unknowns are reported separately rather than folded in.

## Congener clusters are a flag, never a filter

`genus_calls` and `genus_cov_cv` are written on every call so a cluster is visible. One run
calls fourteen *Alternaria* species at coverage 5.6 to 7.3%: fourteen genuinely co-infecting
congeners would differ wildly in abundance, and the uniformity is the diagnostic. But every one
of the 96 runs a collapse rule would flip is *Puccinia*, mostly *striiformis* with *triticina*,
which differ 46.9% against 6.9% coverage and are separable by exact k-mers across a 16 point
ANI gap. Collapsing those discards a real co-infection. No per-detection threshold catches the
difference, because each call is individually well supported. `04_align` resolves it.
"""

import argparse
import collections
import csv
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _reports import (SPECIES_RANK, absolute_floor_ok, coverage, kmers_per_spot,  # noqa: E402
                      load_inspect, parse_report, redundancy, saturation)

HERE = Path(__file__).resolve().parent
MODULE = HERE.parent
RELEASE = MODULE.parent
THRESHOLDS = MODULE / "data" / "inputs" / "thresholds.tsv"
REPORTS = MODULE / "data" / "raw" / "reports"
COHORT = RELEASE / "01_cohort" / "data" / "outputs" / "supp_data_s2_cohort.tsv"
REFERENCE = RELEASE / "02_reference" / "data" / "outputs" / "supp_data_s3_reference.tsv"
EXTRA_CDS = RELEASE / "02_reference" / "data" / "curation" / "extra_cds_db2.tsv"
INSPECT_DIR = RELEASE / "02_reference" / "data" / "outputs"

ARMS = ("db1p", "db2")

CAND_COLS = ["run", "BioProject", "arm", "taxid", "taxon", "rank", "reads", "minimizers",
             "distinct", "ceiling", "coverage", "redundancy", "saturation", "fsp_frac",
             "judgeable", "call", "reason"]
DET_COLS = ["run", "BioProject", "host", "setting", "arm", "taxid", "taxon", "genus",
            "kingdom", "taxon_class", "reads", "distinct", "ceiling", "coverage",
            "redundancy", "saturation", "fsp_frac", "judgeable", "on_wheat",
            "is_plant_pathogen", "is_wheat_pathogen", "declared", "declared_rank",
            "declared_pathogens", "declared_status", "genus_calls", "genus_cov_cv"]


def read(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_thresholds(path: Path, arm: str) -> dict:
    """The operating point for this arm, from the authored table. `arm` beats `all`."""
    out = {}
    for r in read(path):
        if r["arm"] in ("all", arm):
            if r["arm"] == "all" and r["parameter"] in out:
                continue
            out[r["parameter"]] = float(r["value"])
            if r["arm"] == arm:
                out[f"{r['parameter']}__status"] = r["status"]
    missing = {"distinct_floor", "coverage_floor", "saturation_floor", "min_redundancy",
               "max_fsp"} - set(out)
    if missing:
        sys.exit(f"{path.name} is missing {', '.join(sorted(missing))} for arm {arm}")
    return out


def load_cohort() -> dict:
    """{run: the cohort row, plus kmers_per_spot from its own library}."""
    out = {}
    for r in read(COHORT):
        r["kmers_per_spot"] = kmers_per_spot(r.get("base_count"), r.get("read_count"),
                                             r.get("library_layout", ""))
        out[r["run"]] = r
    return out


def load_reference() -> dict:
    """{taxid: the S3 row}. The labels, read off a call and never applied to the candidates."""
    return {int(r["ncbi_taxid"]): r for r in read(REFERENCE) if r["ncbi_taxid"].isdigit()}


def load_predicted(arm: str) -> set:
    """Species taxids whose CDS is gene-predicted. Competitors only, never a positive call.

    Read from `extra_cds_db2.tsv` and not from the build plan, because the plan deliberately
    holds only CDS that exists to download: predicted CDS enters the build through
    `build_db.py --extra-cds`, and that table is the one record of which taxa it covered.

    Refuses rather than returning an empty set for db2: a pass that cannot name the predicted
    taxa would report them as detections, which is the one failure this rule exists to stop.
    """
    taxids = set()
    if EXTRA_CDS.is_file():
        taxids = {int(r["taxid"]) for r in read(EXTRA_CDS) if r["taxid"].isdigit()}
    if arm == "db2" and not taxids:
        sys.exit(f"no usable {EXTRA_CDS.name}, so the predicted species cannot be told from "
                 "the rest of db2.\n"
                 "Predicted CDS is a competitor only and must never support a positive call, "
                 "so this refuses rather than guessing.\n"
                 "That table is written by 02_reference/scripts/plan_extra_cds.py, on the "
                 "machine holding research/predict.")
    return taxids


def normalise(name: str) -> str:
    """A declared or database name reduced to genus + specific epithet, lowercased.

    A declaration is routinely a forma specialis (`Puccinia striiformis f. sp. tritici`) while
    the database node is the species, so the comparison is made at species rank. The genus
    alone is kept separately: a genus match is reported as such and never as a species match.
    """
    toks = [t for t in name.replace("[", "").replace("]", "").split() if t]
    return " ".join(toks[:2]).lower()


GENUS_ONLY = {"sp.", "spp.", "sp", "spp"}


def declared_sets(value: str) -> tuple[set, set]:
    """(species, genera) a study declared, from the curated `; `-separated value.

    **A GENUS IS RECORDED ONLY WHERE THE DECLARATION NAMES NO SPECIES.** This was wrong until
    2026-11-23 and it silently erased the chapter's own result. The genus was added for every
    declaration, including full binomials, so `declaration()`'s genus fallback fired for a
    congener of the declared species: a study declaring *Puccinia striiformis* f. sp.
    *tritici* also marked every *P. triticina* and *P. graminis* call as declared. The tell
    was that not one *Puccinia* call in the whole cohort came out undeclared, while
    *Zymoseptoria* and *Blumeria*, whose genera no study declares, behaved correctly.

    A study that names a species has not reported its congeners. Wheat rust surveys declare
    exactly one rust each, so the three *Puccinia* on wheat are precisely the case where the
    distinction decides the result.

    The genus branch is kept because a declaration CAN be genus-only ("Fusarium spp."), in
    which case any species of that genus is genuinely reported. `main` prints how many such
    declarations exist so the branch cannot quietly stop being exercised.
    """
    species, genera = set(), set()
    for part in value.split(";"):
        part = part.strip()
        if not part:
            continue
        toks = part.split()
        if len(toks) >= 2 and toks[1].lower() not in GENUS_ONLY:
            species.add(normalise(part))
        else:
            genera.add(toks[0].lower())
    return species, genera


def declaration(crow: dict, taxon: str, genus: str, is_host: bool) -> tuple[str, str]:
    """Was this taxon reported by the study that generated the run.

    `n/a` for the plant host, which no study declares as its pathogen and which would
    otherwise land in the unreported count. `unknown` wherever the declaration has not been
    curated: an uncurated project cannot support a claim that a detection went unreported, and
    folding those in would manufacture the headline.
    """
    if is_host:
        return "n/a", ""
    if crow.get("declared_status", "not_curated") == "not_curated":
        return "unknown", ""
    species, genera = declared_sets(crow.get("declared_pathogens", ""))
    if normalise(taxon) in species:
        return "declared", "species"
    if genus.lower() in genera:
        return "declared", "genus"
    return "undeclared", ""


def fsp_fractions(path: Path) -> dict:
    """{species taxid: the largest share of its reads sitting under one sub-species node}.

    Inert in both arms, which the run reports rather than hiding: neither contains a
    sub-species node, so every value is 0. It is computed anyway because the smearing it
    catches is a property of the build, and a build that admits one strain-level assembly makes
    it live again with no other code change.
    """
    out, totals, sub = {}, {}, collections.defaultdict(list)
    for rec in parse_report(path):
        if rec["rank"] == SPECIES_RANK:
            totals[rec["taxid"]] = rec["reads_clade"]
        elif rec["rank"].startswith(SPECIES_RANK) and len(rec["rank"]) > 1:
            parent = next((tid for rank, tid, _ in reversed(rec["lineage"])
                           if rank == SPECIES_RANK), None)
            if parent is not None:
                sub[parent].append(rec["reads_clade"])
    for taxid, total in totals.items():
        out[taxid] = (max(sub[taxid]) / total) if (total and sub.get(taxid)) else 0.0
    return out


def call_one(cov, dist, sat, fsp, judgeable, th, is_host=False) -> tuple[str, str]:
    """The four tests, in order, each written positively. Returns (call, reason).

    **THE HOST IS RANKED, NOT THRESHOLDED**, which 02_classify settled and this function
    was ignoring. Every plant library has a host by construction, so "is there a host" is
    not the question; "which plant" is. The abundance floors are therefore not applied to a
    host row -- they cost 77 of them, including the three uninoculated host-only controls
    where the right answer is known and is wheat -- but dropping the floors ALONE is wrong
    in the other direction: it admitted every plant clearing the absolute floor and called
    Hordeum vulgare in 1,141 runs of a wheat cohort, at coverage down to 0.0002. That is the
    one-directional cross-mapping between cereal genomes this project already documented.

    So the floors are replaced for the host, not removed: `host_rank` keeps the single
    highest-read plant per run. Ranking by reads agrees with the curated host 91.2% of the
    time, which is the basis 02_classify established.
    """
    if not absolute_floor_ok(dist, th["distinct_floor"]):
        return "rejected", "distinct"
    if is_host:
        # The host's call is decided per RUN by rank, not per row by threshold; `host_rank`
        # below picks the winner. Returning "detected" here only means "eligible".
        return "detected", "host, ranked per run"
    if not (cov is not None and cov >= th["coverage_floor"]):
        return "rejected", "coverage"
    if judgeable and not (sat is not None and sat >= th["saturation_floor"]):
        return "rejected", "saturation"
    if not (fsp is not None and fsp <= th["max_fsp"]):
        return "rejected", "forma_specialis"
    return "detected", "all floors" if judgeable else "floors, saturation silent"


def taxon_class(ref_row: dict, predicted: bool) -> str:
    if predicted:
        return "predicted_competitor"
    if not ref_row:
        return "not_in_reference"
    if (ref_row.get("kingdom") or "") not in ("Fungi", ""):
        return "host"
    if ref_row.get("is_wheat_pathogen") in ("yes", "1", "true"):
        return "wheat_pathogen"
    if ref_row.get("is_plant_pathogen") in ("yes", "1", "true"):
        return "plant_pathogen"
    return "resident"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", required=True, choices=ARMS)
    ap.add_argument("--reports-dir", help="default data/raw/reports/<arm>")
    ap.add_argument("--cov-floor", type=float, help="override the arm's coverage floor")
    ap.add_argument("--sat-floor", type=float, help="override the saturation floor")
    ap.add_argument("--distinct-floor", type=int, help="override the absolute floor")
    ap.add_argument("--sensitivity", action="store_true",
                    help="report the call count across a range of coverage floors, write nothing")
    args = ap.parse_args()

    arm = args.arm
    reports = Path(args.reports_dir) if args.reports_dir else REPORTS / arm
    inspect_path = INSPECT_DIR / f"db_inspect_{arm}.tsv"
    th = load_thresholds(THRESHOLDS, arm)
    if args.cov_floor is not None:
        th["coverage_floor"] = args.cov_floor
    if args.sat_floor is not None:
        th["saturation_floor"] = args.sat_floor
    if args.distinct_floor is not None:
        th["distinct_floor"] = float(args.distinct_floor)

    inspect = load_inspect(inspect_path)
    cohort = load_cohort()
    ref = load_reference()
    predicted = load_predicted(arm)

    print(f"arm {arm}")
    print(f"  floors: distinct >= {th['distinct_floor']:.0f}, "
          f"coverage >= {th['coverage_floor']:.4f} "
          f"({th.get('coverage_floor__status', 'provisional')}), "
          f"saturation >= {th['saturation_floor']:.3f}, fsp <= {th['max_fsp']:.2f}")
    sub_nodes = sum(1 for e in inspect.values()
                    if e["rank"].startswith("S") and len(e["rank"]) > 1)
    print(f"  f. sp. feature: {'LIVE' if sub_nodes else 'INERT'} "
          f"({sub_nodes} sub-species nodes in {inspect_path.name})")
    print(f"  predicted competitors: {len(predicted)} species, never callable")

    if not reports.is_dir():
        sys.exit(f"no reports at {reports}; run classify.py --arm {arm} first")

    cands, dets = [], []
    no_minimizers = seen_runs = missing_cohort = 0
    for proj_dir in sorted(p for p in reports.iterdir() if p.is_dir()):
        for report in sorted(proj_dir.glob("*.report.tsv")):
            run = report.name.split(".")[0]
            crow = cohort.get(run)
            if not crow:
                missing_cohort += 1
                continue
            seen_runs += 1
            kps = crow["kmers_per_spot"]
            fsp = fsp_fractions(report)
            host_rows = []              # (reads, row) for every plant clearing the floors
            for rec in parse_report(report):
                if rec["rank"] != SPECIES_RANK or not rec["reads_clade"]:
                    continue
                taxid, reads, dist = rec["taxid"], rec["reads_clade"], rec["distinct"]
                if dist is None:
                    no_minimizers += 1
                cov = coverage(dist, taxid, inspect)
                red = redundancy(reads, taxid, inspect, kps) if kps else None
                sat = saturation(dist, reads, taxid, inspect, kps) if kps else None
                judgeable = red is not None and red >= th["min_redundancy"]
                f = fsp.get(taxid, 0.0)
                klass = taxon_class(ref.get(taxid, {}), taxid in predicted)
                call, reason = call_one(cov, dist, sat, f, judgeable, th,
                                        is_host=klass == "host")
                if taxid in predicted:
                    call, reason = "competitor_only", "predicted CDS, never a positive call"
                e = inspect.get(taxid)
                row = {"run": run, "BioProject": crow["BioProject"], "arm": arm,
                       "taxid": taxid, "taxon": rec["name"], "rank": rec["rank"],
                       "reads": reads, "minimizers": rec["minimizers"], "distinct": dist,
                       "ceiling": e["clade"] if e else "",
                       "coverage": f"{cov:.6g}" if cov is not None else "",
                       "redundancy": f"{red:.6g}" if red is not None else "",
                       "saturation": f"{sat:.6g}" if sat is not None else "",
                       "fsp_frac": f"{f:.4g}", "judgeable": int(judgeable),
                       "call": call, "reason": reason}
                cands.append(row)
                if klass == "host":
                    host_rows.append((reads, row))
                    continue            # the host's verdict is settled after the run
                if call == "detected":
                    rr = ref.get(taxid, {})
                    status = crow.get("declared_status", "not_curated")
                    gen = rr.get("genus") or rec["name"].split()[0]
                    dec, dec_rank = declaration(crow, rec["name"], gen, klass == "host")
                    dets.append({**{k: v for k, v in row.items() if k in DET_COLS},
                                 "host": crow.get("host", ""), "setting": crow.get("setting", ""),
                                 "genus": gen, "kingdom": rr.get("kingdom", ""),
                                 "taxon_class": klass,
                                 "on_wheat": rr.get("on_wheat", ""),
                                 "is_plant_pathogen": rr.get("is_plant_pathogen", ""),
                                 "is_wheat_pathogen": rr.get("is_wheat_pathogen", ""),
                                 "declared": dec, "declared_rank": dec_rank,
                                 "declared_pathogens": crow.get("declared_pathogens", ""),
                                 "declared_status": status})

            # ── the host, ranked, one per run ─────────────────────────────────
            # By READS, which is what 02_classify established agrees with the curated host
            # 91.2% of the time. Every other plant that cleared the absolute floor is
            # written to the candidate table as rejected, with the reason, so the runner-up
            # is visible rather than silently absent: in a wheat cohort that runner-up is
            # Hordeum vulgare in most runs and it is cross-mapping, not barley.
            if host_rows:
                host_rows.sort(key=lambda hr: -hr[0])
                for n, (_, hrow) in enumerate(host_rows):
                    if n:
                        hrow["call"], hrow["reason"] = "rejected", "host, outranked"
                        continue
                    rr = ref.get(hrow["taxid"], {})
                    dets.append({**{k: v for k, v in hrow.items() if k in DET_COLS},
                                 "host": crow.get("host", ""),
                                 "setting": crow.get("setting", ""),
                                 "genus": rr.get("genus") or hrow["taxon"].split()[0],
                                 "kingdom": rr.get("kingdom", ""), "taxon_class": "host",
                                 "on_wheat": rr.get("on_wheat", ""),
                                 "is_plant_pathogen": rr.get("is_plant_pathogen", ""),
                                 "is_wheat_pathogen": rr.get("is_wheat_pathogen", ""),
                                 "declared": "n/a", "declared_rank": "",
                                 "declared_pathogens": crow.get("declared_pathogens", ""),
                                 "declared_status": crow.get("declared_status",
                                                             "not_curated")})

    if args.sensitivity:
        print("\nsensitivity to the coverage floor, nothing written:")
        print(f"  {'floor':>7s} {'calls':>8s} {'runs':>6s} {'undeclared':>11s}"
              "   (undeclared counts CALLS with a curated declaration that did "
              "not name them)")
        for floor in (0.005, 0.01, 0.02, 0.03, 0.05, 0.10):
            kept = [c for c in cands if c["call"] != "competitor_only"
                    and call_one(float(c["coverage"]) if c["coverage"] else None,
                                 c["distinct"],
                                 float(c["saturation"]) if c["saturation"] else None,
                                 float(c["fsp_frac"]), bool(c["judgeable"]),
                                 {**th, "coverage_floor": floor},
                                 # candidates_<arm>.tsv carries no taxon_class column, so
                                 # this must be derived, not looked up
                                 is_host=taxon_class(ref.get(int(c["taxid"]), {}),
                                                     False) == "host")[0] == "detected"]
            und = [c for c in kept
                   if declaration(cohort[c["run"]], c["taxon"],
                                  ref.get(int(c["taxid"]), {}).get("genus")
                                  or c["taxon"].split()[0],
                                  taxon_class(ref.get(int(c["taxid"]), {}), False) == "host"
                                  )[0] == "undeclared"]
            print(f"  {floor:7.3f} {len(kept):8d} {len({c['run'] for c in kept}):6d} "
                  f"{len(und):11d}")
        print("\nThe floor is anchored to the off-host null, not chosen off this table. "
              "See calibrate_floors.py.")
        return

    # genus cluster flags, computed over the CALLS in each run: a count and the spread of
    # coverage within the genus. A tight cluster of congeners at near-identical coverage is one
    # organism's reads spread over several references; it is flagged and never filtered.
    by_run_genus = collections.defaultdict(list)
    for d in dets:
        by_run_genus[(d["run"], d["genus"])].append(float(d["coverage"]))
    for d in dets:
        covs = by_run_genus[(d["run"], d["genus"])]
        d["genus_calls"] = len(covs)
        d["genus_cov_cv"] = (f"{statistics.pstdev(covs)/statistics.fmean(covs):.4g}"
                             if len(covs) > 1 and statistics.fmean(covs) else "")

    cand_path = MODULE / "data" / "raw" / f"candidates_{arm}.tsv"
    det_path = MODULE / "data" / "outputs" / f"detections_{arm}.tsv"
    for path, cols, rows in ((cand_path, CAND_COLS, cands), (det_path, DET_COLS, dets)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    tally = collections.Counter(c["call"] for c in cands)
    rej = collections.Counter(c["reason"] for c in cands if c["call"] == "rejected")
    cls = collections.Counter(d["taxon_class"] for d in dets)
    dec = collections.Counter(d["declared"] for d in dets)
    print(f"\n{seen_runs} runs, {len(cands):,} candidates")
    if missing_cohort:
        print(f"  {missing_cohort} reports for runs not in the cohort, skipped")
    if no_minimizers:
        print(f"  WARNING {no_minimizers:,} rows carry no distinct-minimizer column: those "
              "reports were written without --report-minimizer-data and cannot be scored")
    print(f"  calls {tally['detected']:,}  rejected {tally['rejected']:,} "
          f"({', '.join(f'{k} {v:,}' for k, v in rej.most_common())})")
    print(f"  competitor_only {tally['competitor_only']:,}")
    print(f"  by class:   {', '.join(f'{k} {v:,}' for k, v in cls.most_common())}")
    print(f"  by declaration: {', '.join(f'{k} {v:,}' for k, v in dec.most_common())}")
    print(f"  runs with a fungal call: "
          f"{len({d['run'] for d in dets if d['taxon_class'] != 'host'}):,} of {seen_runs:,}")
    print(f"\n-> {cand_path.relative_to(MODULE)}  {len(cands):,} rows")
    print(f"-> {det_path.relative_to(MODULE)}  {len(dets):,} calls")
    print("\nNo rate out of this step is a prevalence. The cohort cannot support one.")


if __name__ == "__main__":
    main()
