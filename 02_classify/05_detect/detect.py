#!/usr/bin/env python3
"""detect.py — which Kraken2 assignments are real organisms, and which are artefacts.

    python 02_classify/05_detect/detect.py
    python 02_classify/05_detect/detect.py --sat-floor 0.05 --cov-floor 0.01
    python 02_classify/05_detect/detect.py --sensitivity

Writes `silver/_classified.tsv` (every candidate, with its metrics) and
`gold/detections.tsv` + `gold/coinfection.tsv` (the calls, and the per-sample summary).

## What carries across from 05_filter, and what does not

CARRIED, because they are knowledge rather than code:

  * **The `unbiased` f. sp. feature, WHICH IS INERT UNDER db_v4 AND SAYS SO AT STARTUP.**
    A detection whose reads pile onto ONE sub-taxon (S1/S2 under the species) is LCA smearing,
    not a genome-wide infection, and saturation cannot see it -- the reads genuinely do spread
    across the k-mer space of the sliver they came from. *Puccinia graminis* on wheat forced it
    in db_v3: two overlapping clouds, one rising and one flat, split by forma specialis.
    Measured 2026-10-20, **db_v4 contains no sub-species nodes whatsoever** (1,547 species,
    zero S1/S2, against db_v3's 222), so `fsp_frac` is 0 on all 300,061 candidates and the test
    never fires. It is RETAINED rather than deleted because the smearing it catches is a
    property of the database build, not of the method: admit one strain-level assembly and it
    matters again. A feature that silently does nothing is how knowledge rots, so the run
    reports whether it is live.
  * **Per-host application.** The same organism is a real infection on its own host and a
    cross-mapping artefact on another; *Sclerotinia* is real on canola and noise on wheat,
    *P. striiformis* real on wheat and noise on maize. Host comes from the CURATED value in
    cohort_samples.tsv, not from a read-based call: the read-based call has a documented
    one-directional bias in rust surveys.
  * **Two thresholds doing two jobs, never conflated.** The abundance floor asks whether the
    organism's DNA is there; the shape test asks whether the signal has the right structure.

REPLACED: `tightness`, `rising` and `on_line`. All three were features of a FITTED cloud --
a per-host-per-species log-log regression over >=15 points spanning >=0.7 decades -- and
`saturation` answers the same question analytically from a single report. That matters for
three reasons beyond elegance:

  1. **No NaN.** The fit needed a populated cloud; 3,382 of 24,467 rows (13.8%) in the db_v3
     table had `score = nan` and `cls = unclassified` because no curve could be fitted. NaN
     compares False against everything, so `if score < FLOOR: continue` ADMITS them -- that
     inflated a cereal species count from 32 to 44 before it was caught, and 393 of those rows
     cleared the 1% coverage floor. An analytic null has nothing to fail to fit.
  2. **No refitting.** Adding samples re-fits every cloud and can move an existing call.
  3. **An explicit null** instead of five hand-weighted features. The old weights were a
     stated provisional operating point precisely because automated fits failed: fitting to
     PHI-base labels let coverage dominate and crushed low-abundance reals, while dropping
     coverage over-rejected obvious reals.

Agreement was measured before replacing anything: Spearman(score, saturation) = 0.872 over
9,442 (run, species) pairs, holding within redundancy decades (+0.798, +0.867, +0.821).

## The criterion

For each species-rank taxon in each run, with M = that taxon's clade minimizers in db_v4:

    redundancy = reads * kmers_per_spot / M     how many times over the reads could cover it
    coverage   = distinct / M                   how much of it was actually seen
    saturation = coverage / (1 - exp(-redundancy))

A detection is called when ALL of:

    coverage   >= COV_FLOOR       the abundance floor, 1%: is the DNA there
    distinct   >= MIN_DISTINCT    the absolute floor, and NOT optional -- 1% of db_v4's
                                  smallest species is 29 k-mers, reachable by cross-mapping
                                  alone, and 20 species have ceilings under 100,000
    saturation >= SAT_FLOOR       the shape test: are the reads spread or piled
    fsp_frac   <= MAX_FSP         not dominated by one forma specialis

SATURATION MUST BE READ AGAINST REDUNDANCY, NEVER AS A SCALAR. For small x,
1 - exp(-x) ~= x, so coverage ~= redundancy mechanically and every taxon scores near 1
whether it is real or three cross-mapped loci. Median saturation falls 0.67 -> 0.23 -> 0.053
-> 0.019 -> 0.002 across redundancy decades, which is the mechanism and not biology. Below
MIN_REDUNDANCY the statistic is silent and the call rests on the floors alone. In practice
nothing clearing a 1% coverage floor is in that zone -- coverage ~= redundancy there -- but
the guard matters the moment the coverage floor is lowered.

## The thresholds are a provisional operating point, exactly as the old weights were

SAT_FLOOR = 0.05 sits below the 25th percentile (0.091) of detections the db_v3 score called
real, and above the 460 calls that clear coverage with piled reads. It is NOT calibrated: real
calls span 0.07-0.56 rather than sitting near 1, because the null assumes uniform sampling and
RNA-seq expression is not uniform -- reads concentrate on highly expressed genes. If expression
skew varies by organism or infection stage then this floor is not constant across taxa, which
the literature search flagged independently. **Calibrating it needs the axenic controls**, which
give the saturation distribution of a known-present organism in a real RNA-seq library. Until
then, run --sensitivity and report the range.

## CONGENER INFLATION, and why the genus count is reported alongside the species count

Found 2026-10-20 on the first full run. ERR1607767 calls FOURTEEN *Alternaria* species, every
one at coverage 5.6-7.3% and saturation 0.070-0.080. Fourteen genuinely co-infecting congeners
would differ wildly in abundance; one organism's reads distributed across fourteen
near-identical references looks exactly like this. The uniformity IS the diagnostic.

Saturation cannot catch it, and that is not a failure of the statistic: the reads genuinely do
spread across each reference, because the references are nearly the same sequence. Nor can any
per-detection threshold, because each call is individually well-supported. It is a property of
the database's congener density, and db_v4 deliberately has a lot of that -- breadth was added
to give the classifier competitors.

Affected genera, by runs carrying >=3 species of one genus: *Alternaria* 260, *Cladosporium*
61, *Bipolaris* 30, *Cercospora* 27, *Fusarium* 22, *Puccinia* 16, *Didymella* 16,
*Curvularia* 15.

So `coinfection.tsv` carries BOTH counts and neither is hidden:

    n_pathogens        species-level. An UPPER bound on co-infection.
    n_pathogen_genera  genus-level. The conservative reading, and the one to quote until
                       alignment resolves congeners.

Measured here: 701/1,177 runs (60%) co-infected by species, 605/1,177 (51%) by genus.

Resolving it properly is `04_align`'s job -- competitive EM assigns a read to one reference
instead of counting it against all plausible ones -- and this gives that module a defined
target set rather than "a stratified sample". Note it may not fully succeed: *Puccinia* species
are separable by exact k-mers (a 16-point ANI gap) but formae speciales are not, and whether
*Alternaria* sections separate is unmeasured.

## Pathogen status is a LABEL, applied after detection. Never a filter before it.

The db_v3 detect.py skipped any taxon outside PHI-base's 632 fungal taxa, which under db_v4
would discard every resident the database was expanded to hold. Here everything is detected
first and labelled second, so a resident, an endophyte and a pathogen are all measurable and
the pathogen question is asked of the output rather than of the input.
"""

import argparse
import csv
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _layout import (CLASSIFIED, COHORT_SAMPLES, COINFECTION, DB_INSPECT, DETECTIONS,
                     GOLD, HOLOBASE, LIT_RUNS, REPORTS, SILVER)  # noqa: E402
from _reports import load_inspect, parse_report  # noqa: E402

K = 35
COV_FLOOR = 0.01        # the published abundance floor
MIN_DISTINCT = 10_000   # the absolute companion; see the docstring
SAT_FLOOR = 0.05        # provisional; calibrate with the axenic controls
MAX_FSP = 0.80          # a detection >80% from one forma specialis is LCA smearing
MIN_REDUNDANCY = 1e-3   # below this, saturation is silent

CLS_COLS = ["run", "BioProject", "biosample", "host", "species", "taxid", "kingdom", "guild",
            "reads", "distinct", "db_clade", "redundancy", "coverage", "saturation",
            "fsp_frac", "judgeable", "call", "reason"]


def labels(db: Path) -> dict:
    """ncbi_taxid -> (kingdom, guilds). Applied AFTER detection, never before."""
    out = {}
    if not db.exists():
        print(f"  no holobase at {db}; guild labels will be blank", file=sys.stderr)
        return out
    con = sqlite3.connect(db)
    for tid, king, guild in con.execute(
            "SELECT t.ncbi_taxid, t.kingdom, g.guild FROM taxon t "
            "LEFT JOIN taxon_guild g ON g.taxon_id = t.taxon_id "
            "WHERE t.ncbi_taxid IS NOT NULL"):
        if tid is None:
            continue
        e = out.setdefault(int(tid), [king or "", set()])
        if guild:
            e[1].add(guild)
    return out


def subtaxa(inspect_path: Path) -> dict:
    """species taxid -> [S1/S2 taxids beneath it], for the f. sp. feature.

    Read from the inspect table's own indentation order rather than from a taxonomy dump,
    because that is where the ranks db_v4 actually tagged are recorded.
    """
    out, cur = defaultdict(list), None
    for line in open(inspect_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        rank = f[3].strip()
        try:
            tid = int(f[4])
        except ValueError:
            continue
        if rank == "S":
            cur = tid
            out.setdefault(cur, [])
        elif rank.startswith("S") and cur is not None:
            out[cur].append(tid)
    return dict(out)


def kmers_per_spot(read_len: float, layout: str) -> float:
    """(L-K+1) summed over a spot's mates. ENA's read_len is bases PER SPOT."""
    L = read_len / 2 if layout == "PAIRED" else read_len
    if L <= K:
        return 0.0
    per = L - K + 1
    return per * 2 if layout == "PAIRED" else per


def call_one(cov, dist, sat, fsp, judgeable, cov_floor, sat_floor) -> tuple[str, str]:
    """All tests POSITIVE. `sat` can be None, and `None < floor` raises while
    `not (None >= floor)` would silently admit it -- the NaN trap, one layer down."""
    if not (dist is not None and dist >= MIN_DISTINCT):
        return "no", "below absolute floor"
    if not (cov is not None and cov >= cov_floor):
        return "no", "below coverage floor"
    if judgeable and not (sat is not None and sat >= sat_floor):
        return "no", "reads piled, not spread"
    if fsp is not None and fsp > MAX_FSP:
        return "no", "dominated by one forma specialis"
    return "yes", ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cov-floor", type=float, default=COV_FLOOR)
    ap.add_argument("--sat-floor", type=float, default=SAT_FLOOR)
    ap.add_argument("--reports", default=str(REPORTS))
    ap.add_argument("--inspect", default=str(DB_INSPECT))
    ap.add_argument("--sensitivity", action="store_true",
                    help="report the call count across a grid of both floors and stop")
    args = ap.parse_args()

    insp = load_inspect(Path(args.inspect))
    subs = subtaxa(Path(args.inspect))
    n_sub = sum(len(v) for v in subs.values())
    print(f"  f. sp. feature: {'LIVE, ' + str(n_sub) + ' sub-taxa' if n_sub else 'INERT'} "
          f"-- this database has {n_sub} S1/S2 nodes under {len(subs)} species")
    lab = labels(HOLOBASE)
    # TWO tables, and both are needed. cohort_samples.tsv is the membership authority and
    # carries crop/host/setting, but it is keyed on Run and has NO BioSample and NO
    # declared_pathogens columns -- those live in runs.tsv (cols 22-23). Reading the declared
    # set off the wrong table yields a silently empty `dec`, which would make every detection
    # look undeclared.
    cohort = {r["Run"]: r for r in csv.DictReader(open(COHORT_SAMPLES), delimiter="\t")}
    litruns = {r["Run"]: r for r in csv.DictReader(open(LIT_RUNS), delimiter="\t")}
    depth = {r["Run"]: r for r in csv.DictReader(open(SILVER / "run_depth.tsv"),
                                                 delimiter="\t")}
    root = Path(args.reports)

    rows = []
    for pdir in sorted(d for d in root.iterdir() if d.is_dir()):
        for rep in sorted(pdir.glob("*.report.tsv")):
            run = rep.name.split(".")[0]
            d = depth.get(run)
            if not d:
                continue
            kps = kmers_per_spot(float(d["read_len"] or 0), d["layout"])
            if not kps:
                continue
            recs = list(parse_report(rep))
            sub_reads = {r["taxid"]: r["reads_clade"] for r in recs}
            meta = cohort.get(run, {})
            for rec in recs:
                if rec["rank"] != "S" or not rec["reads_clade"] or rec["distinct"] is None:
                    continue
                e = insp.get(rec["taxid"])
                if not e or not e["clade"]:
                    continue
                m = e["clade"]
                red = rec["reads_clade"] * kps / m
                cov = rec["distinct"] / m
                exp = m * (1.0 - math.exp(-min(700.0, red)))
                sat = rec["distinct"] / exp if exp else None
                sr = sum(sub_reads.get(s, 0) for s in subs.get(rec["taxid"], []))
                fsp = sr / rec["reads_clade"] if rec["reads_clade"] else 0.0
                judge = red >= MIN_REDUNDANCY
                kind, guilds = lab.get(rec["taxid"], ["", set()])
                call, why = call_one(cov, rec["distinct"], sat, fsp, judge,
                                     args.cov_floor, args.sat_floor)
                rows.append({
                    "run": run, "BioProject": pdir.name,
                    "biosample": litruns.get(run, {}).get("BioSample", ""),
                    "host": meta.get("host", ""), "species": rec["name"],
                    "taxid": rec["taxid"], "kingdom": kind,
                    "guild": "|".join(sorted(guilds)),
                    "reads": rec["reads_clade"], "distinct": rec["distinct"],
                    "db_clade": m, "redundancy": f"{red:.6g}", "coverage": f"{cov:.6g}",
                    "saturation": f"{sat:.6g}" if sat is not None else "",
                    "fsp_frac": f"{fsp:.4f}", "judgeable": "yes" if judge else "no",
                    "call": call, "reason": why,
                })

    if args.sensitivity:
        print("Calls across a grid of both floors. The thresholds are a PROVISIONAL "
              "operating point;\nreport this range until the axenic controls calibrate "
              "the saturation floor.\n")
        print(f"{'cov floor':>10s} " + " ".join(f"{s:>9.2f}" for s in
                                                (0.0, 0.02, 0.05, 0.10, 0.20)))
        print(f"{'':>10s} " + " ".join(f"{'sat=' + str(s):>9s}" for s in
                                       (0.0, 0.02, 0.05, 0.10, 0.20)))
        for cf in (0.001, 0.005, 0.01, 0.02, 0.05):
            cells = []
            for sf in (0.0, 0.02, 0.05, 0.10, 0.20):
                n = sum(1 for r in rows
                        if call_one(float(r["coverage"]), r["distinct"],
                                    float(r["saturation"]) if r["saturation"] else None,
                                    float(r["fsp_frac"]), r["judgeable"] == "yes",
                                    cf, sf)[0] == "yes")
                cells.append(f"{n:9,d}")
            print(f"{cf:10.3f} " + " ".join(cells))
        return

    SILVER.mkdir(parents=True, exist_ok=True)
    GOLD.mkdir(parents=True, exist_ok=True)
    with open(CLASSIFIED, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CLS_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    det = [r for r in rows if r["call"] == "yes"]
    with open(DETECTIONS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CLS_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(det)

    # ── per-sample summary ──────────────────────────────────────────────────────────────────
    # Pathogen status is read here, from the LABEL, after everything has been detected.
    def is_path(r):
        return any("Plant Pathogen" in g for g in r["guild"].split("|") if g)

    by_run = defaultdict(list)
    for r in det:
        if r["kingdom"] == "Fungi":
            by_run[r["run"]].append(r)
    runs = sorted({r["run"] for r in rows})
    summ = []
    for run in runs:
        rs = by_run.get(run, [])
        meta = cohort.get(run, {})
        declared = litruns.get(run, {}).get("declared_pathogens") or ""
        dec = {" ".join(x.lower().replace("f. sp.", " ").split()[:2])
               for x in declared.split(";") if x.strip()}
        paths = [r for r in rs if is_path(r)]
        undec = [r for r in paths
                 if " ".join(r["species"].lower().split()[:2]) not in dec]
        # Genus is the first token of the species name. Crude, but these are NCBI scientific
        # names from the database's own taxonomy, so the first token IS the genus; parsing it
        # out of a taxonomy dump would add a dependency for no gain.
        genera = {r["species"].split()[0] for r in paths}
        # The congener-cluster signature: >=3 species of one genus in one run.
        clusters = {g for g in genera
                    if sum(1 for r in paths if r["species"].split()[0] == g) >= 3}
        summ.append({
            "run": run, "BioProject": meta.get("BioProject", ""),
            "crop": meta.get("crop", ""), "host": meta.get("host", ""),
            "sampling_selection": meta.get("sampling_selection", ""),
            "n_fungi": len(rs), "n_pathogens": len(paths),
            "n_pathogen_genera": len(genera), "n_undeclared": len(undec),
            "coinfected": "yes" if len(paths) >= 2 else "no",
            "coinfected_genus": "yes" if len(genera) >= 2 else "no",
            "congener_clusters": ";".join(sorted(clusters)),
            "pathogens": ";".join(sorted(r["species"] for r in paths)),
            "undeclared": ";".join(sorted(r["species"] for r in undec)),
        })
    with open(COINFECTION, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summ)

    print(f"-> {CLASSIFIED.name}: {len(rows):,} candidates over {len(runs):,} runs")
    print(f"-> {DETECTIONS.name}: {len(det):,} calls")
    print(f"-> {COINFECTION.name}: {len(summ):,} runs")
    rej = defaultdict(int)
    for r in rows:
        if r["call"] == "no":
            rej[r["reason"]] += 1
    print("\n  rejected because:")
    for why, n in sorted(rej.items(), key=lambda kv: -kv[1]):
        print(f"    {why:36s} {n:8,d}")
    co = sum(1 for s in summ if s["coinfected"] == "yes")
    cog = sum(1 for s in summ if s["coinfected_genus"] == "yes")
    ncl = sum(1 for s in summ if s["congener_clusters"])
    print(f"\n  >=2 plant pathogen SPECIES: {co:,} of {len(summ):,} ({co/len(summ):.1%})"
          f"   <- upper bound")
    print(f"  >=2 plant pathogen GENERA : {cog:,} of {len(summ):,} ({cog/len(summ):.1%})"
          f"   <- conservative, quote this one")
    print(f"  runs with a congener cluster (>=3 species of one genus): {ncl:,}")
    print(f"  >=1 UNDECLARED pathogen: "
          f"{sum(1 for s in summ if s['n_undeclared']):,}")
    print("\n  NOTE: every run classified so far is STAT-retained, so these are NOT "
          "prevalence\n  estimates. See crypt-reads-layout / the 2026-10-20 commit.")


if __name__ == "__main__":
    main()
