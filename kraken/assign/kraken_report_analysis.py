#!/usr/bin/env python3
"""
kraken_report_analysis.py — turn the per-run Kraken2 reports into the co-infection
tables the chapter reports.

What it answers: for each SRA run, which plant pathogens are present above a stated
detection criterion, and which of those the original study never declared. That second
part is the "cryptic" definition, and it is why this cannot be done from the reports
alone: it needs the study's declared pathogens from the metadata module.

Inputs:
  kraken/assign/data/reports/*.txt          per-run Kraken2 reports
  kraken/build/data/db_v3_inspect.tsv       minimizers per taxon IN THE DB (denominator)
  stat/build/data/phibase_db.json           which taxids are pathogens / hosts
  stat/filter/data/runs.tsv                 Run -> BioSample
  metadata/classify/data/samples.tsv        BioSample -> declared pathogen taxids, setting

Outputs (kraken/assign/data/):
  detections.tsv        one row per run x detected pathogen taxon
  run_summary.tsv       one row per run: pathogens found, declared, cryptic
  cryptic_by_setting.tsv aggregate rates, same shape as the STAT table
  threshold_sweep.tsv   how the co-infection rate moves with the floor (--sweep)

On the detection criterion. A flat read floor is not defensible on its own: across the
first db_v3 run the share of runs carrying two or more pathogen species was 98.8% at a
floor of 1 read and still 49.7% at 10,000, so the headline number would be mostly a
choice. Where --report-minimizer-data is present, the primary criterion is instead the
fraction of a taxon's k-mer space actually observed (distinct minimizers seen, over that
taxon's total minimizers in the database). Reads piled on one conserved repeat score low
however many there are. The read floor is kept as a secondary guard, not the basis.

Run from crypt/:
  python kraken/assign/kraken_report_analysis.py --sweep
  python kraken/assign/kraken_report_analysis.py --min-reads 100 --min-kmer-frac 0.01
"""
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
csv.field_size_limit(10 ** 7)          # samples.tsv carries long LLM rationale fields


# ── report parsing ────────────────────────────────────────────────────────────
def parse_report(path):
    """Yield one dict per taxon line, with its lineage.

    Two layouts exist and are told apart per line by column count, never assumed:
        6 cols  pct reads taxReads rank taxid name
        8 cols  pct reads taxReads minimizers distinct rank taxid name   (--report-minimizer-data)
    Indexing blind reads the minimizer count as the rank code and silently drops
    every species.

    Lineage comes from the indentation of the name column, two spaces per level, so
    genus rollup and ancestry need no external taxonomy file.
    """
    stack = []                                   # [(depth, rank, taxid, name)]
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            mini = len(parts) >= 8
            raw_name = parts[7] if mini else parts[5]
            depth = (len(raw_name) - len(raw_name.lstrip(" "))) // 2
            try:
                rec = {
                    "reads_clade":  int(parts[1]),
                    "reads_direct": int(parts[2]),
                    "minimizers":   int(parts[3]) if mini and parts[3].strip().isdigit() else None,
                    "distinct":     int(parts[4]) if mini and parts[4].strip().isdigit() else None,
                    "rank":         (parts[5] if mini else parts[3]).strip(),
                    "taxid":        int((parts[6] if mini else parts[4]).strip()),
                    "name":         raw_name.strip(),
                    "depth":        depth,
                }
            except ValueError:
                continue
            while stack and stack[-1]["depth"] >= depth:
                stack.pop()
            rec["lineage"] = [(s["rank"], s["taxid"], s["name"]) for s in stack]
            stack.append(rec)
            yield rec


def ancestor_at(rec, rank):
    for r, tid, name in rec["lineage"]:
        if r == rank:
            return tid, name
    return None, None


# ── reference data ────────────────────────────────────────────────────────────
def load_reference():
    db = json.loads((ROOT / "stat/build/data/phibase_db.json").read_text())
    pathogens = (set(map(str, db["fungal_to_seed"]))
                 | set(map(str, db["oomycete_to_seed"]))
                 | set(map(str, db["nematode_to_seed"])))
    hosts = set(map(str, db["host_to_seed"]))
    return pathogens, hosts


def load_inspect(path):
    """taxid -> distinct minimizers that taxon holds in the DB (the denominator)."""
    out = {}
    if not path.exists():
        return out
    for rec in parse_report(path):
        # kraken2-inspect writes minimizer counts in the reads columns
        out[rec["taxid"]] = rec["reads_clade"]
    return out


def load_run_to_biosample():
    out = {}
    with open(ROOT / "stat/filter/data/runs.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["Run"].strip()] = row["BioSample"].strip()
    return out


def load_samples():
    out = {}
    with open(ROOT / "metadata/classify/data/samples.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            declared = {t.strip() for t in (row.get("llm_named_pathogens_taxids") or "").replace(";", ",").split(",") if t.strip().isdigit()}
            out[row["BioSample"].strip()] = {
                "BioProject":  row.get("BioProject", ""),
                "declared":    declared,
                "declared_names": row.get("llm_named_pathogens", ""),
                "stress":      row.get("llm_stress", ""),
                "setting":     row.get("llm_study_setting", ""),
                "coinf_intent": row.get("llm_coinfection_intent", ""),
            }
    return out


# ── detection ─────────────────────────────────────────────────────────────────
def detections_for_run(path, pathogens, hosts, db_minimizers, rank,
                       min_reads, min_kmer_frac):
    """Pathogen taxa in one report that clear the criterion, at the given rank.

    Species rank is "S"; Kraken2 also emits S1/S2 for subspecies and formae speciales,
    whose reads are already included in their parent S clade count, so counting only
    "S" avoids double counting a species and its f.sp.
    """
    out = []
    for rec in parse_report(path):
        if rec["rank"] != rank:
            continue
        tid = str(rec["taxid"])
        if tid in hosts or tid not in pathogens:
            continue
        frac = None
        if rec["distinct"] is not None:
            total = db_minimizers.get(rec["taxid"])
            if total:
                frac = rec["distinct"] / total
        if rec["reads_clade"] < min_reads:
            continue
        if frac is not None and frac < min_kmer_frac:
            continue
        g_tid, g_name = ancestor_at(rec, "G")
        out.append({
            # every ancestor, so a study that declared its pathogen at genus (or at
            # f.sp. below the species) still matches its own detection instead of
            # being scored cryptic against itself
            "lineage_taxids": {str(rec["taxid"])} | {str(tid) for _, tid, _ in rec["lineage"]},
            "taxid": rec["taxid"], "name": rec["name"],
            "reads": rec["reads_clade"],
            "distinct_minimizers": rec["distinct"],
            "kmer_frac": round(frac, 6) if frac is not None else "",
            "genus_taxid": g_tid or "", "genus": g_name or "",
        })
    return out


def sweep(report_paths, pathogens, hosts, db_minimizers, out_path):
    """Document how the floor drives the answer, instead of hiding the choice."""
    floors = [1, 10, 50, 100, 500, 1000, 5000, 10000]
    rows = []
    for rank in ("S", "G"):
        counts = {f: [] for f in floors}
        for p in report_paths:
            hits = [r["reads_clade"] for r in parse_report(p)
                    if r["rank"] == rank and str(r["taxid"]) in pathogens
                    and str(r["taxid"]) not in hosts]
            for f in floors:
                counts[f].append(sum(1 for h in hits if h >= f))
        for f in floors:
            v = sorted(counts[f])
            n = len(v) or 1
            rows.append({
                "rank": rank, "min_reads": f,
                "runs": len(v),
                "median_per_run": v[n // 2] if v else 0,
                "runs_ge2": sum(1 for x in v if x >= 2),
                "pct_ge2": round(100 * sum(1 for x in v if x >= 2) / n, 1),
                "runs_zero": sum(1 for x in v if x == 0),
            })
    _write_tsv(out_path, rows)
    return rows


def _write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reports-dir", default="kraken/assign/data/reports")
    ap.add_argument("--inspect", default="kraken/build/data/db_v3_inspect.tsv",
                    help="kraken2-inspect output: minimizers per taxon in the DB")
    ap.add_argument("--out-dir", default="kraken/assign/data")
    ap.add_argument("--min-reads", type=int, default=100,
                    help="secondary guard, not the basis for a call (default: 100)")
    ap.add_argument("--min-kmer-frac", type=float, default=0.01,
                    help="primary criterion: fraction of the taxon's k-mer space seen. "
                         "Ignored for reports without --report-minimizer-data "
                         "(default: 0.01)")
    ap.add_argument("--sweep", action="store_true",
                    help="also write the floor sensitivity table and stop")
    args = ap.parse_args()

    reports_dir = Path(args.reports_dir)
    out_dir = Path(args.out_dir)
    paths = sorted(reports_dir.glob("*.txt"))
    if not paths:
        sys.exit(f"No reports in {reports_dir}")

    pathogens, hosts = load_reference()
    db_minimizers = load_inspect(Path(args.inspect))
    print(f"reports: {len(paths):,} | pathogen taxids: {len(pathogens):,} | "
          f"DB taxa with minimizer counts: {len(db_minimizers):,}")
    if not db_minimizers:
        print("  WARNING: no inspect table, so the k-mer criterion is inactive and "
              "detection falls back to the read floor alone.")

    if args.sweep:
        rows = sweep(paths, pathogens, hosts, db_minimizers,
                     out_dir / "threshold_sweep.tsv")
        print(f"\n{'rank':>5} {'floor':>7} {'median/run':>11} {'runs>=2':>9} {'rate':>7}")
        for r in rows:
            print(f"{r['rank']:>5} {r['min_reads']:>7} {r['median_per_run']:>11} "
                  f"{r['runs_ge2']:>9} {r['pct_ge2']:>6}%")
        print(f"\nwrote {out_dir / 'threshold_sweep.tsv'}")
        return

    run2bs = load_run_to_biosample()
    samples = load_samples()
    detections, run_rows = [], []
    by_setting = defaultdict(lambda: {"runs": 0, "cryptic_runs": 0, "coinf_runs": 0})

    for p in paths:
        run = p.stem
        bs = run2bs.get(run, "")
        meta = samples.get(bs, {})
        declared = meta.get("declared", set())
        for rank, label in (("S", "species"), ("G", "genus")):
            hits = detections_for_run(p, pathogens, hosts, db_minimizers, rank,
                                      args.min_reads, args.min_kmer_frac)
            for h in hits:
                h["is_declared"] = bool(h.pop("lineage_taxids") & declared)
            cryptic = [h for h in hits if not h["is_declared"]]
            for h in hits:
                row = {k: v for k, v in h.items() if k != "is_declared"}
                detections.append({
                    "run": run, "biosample": bs,
                    "bioproject": meta.get("BioProject", ""),
                    "rank": label, **row,
                    "declared": "yes" if h["is_declared"] else "no",
                })
            run_rows.append({
                "run": run, "biosample": bs, "rank": label,
                "n_pathogens": len(hits), "n_cryptic": len(cryptic),
                "coinfection": "yes" if len(hits) >= 2 else "no",
                "declared_taxids": ",".join(sorted(declared)),
                "declared_names": meta.get("declared_names", ""),
                "stress": meta.get("stress", ""),
                "setting": meta.get("setting", ""),
                "coinf_intent": meta.get("coinf_intent", ""),
            })
            if label == "species":
                # A deliberate multi-pathogen experiment is not cryptic by definition,
                # so it is excluded from the rate rather than counted as a finding.
                if meta.get("coinf_intent", "") == "coinf_experiment":
                    continue
                s = by_setting[meta.get("setting", "unknown")]
                s["runs"] += 1
                s["cryptic_runs"] += 1 if cryptic else 0
                s["coinf_runs"] += 1 if len(hits) >= 2 else 0

    _write_tsv(out_dir / "detections.tsv", detections)
    _write_tsv(out_dir / "run_summary.tsv", run_rows)
    setting_rows = [{"setting": k, "runs": v["runs"],
                     "cryptic_runs": v["cryptic_runs"],
                     "pct_cryptic": round(100 * v["cryptic_runs"] / v["runs"], 1) if v["runs"] else 0,
                     "coinf_runs": v["coinf_runs"],
                     "pct_coinf": round(100 * v["coinf_runs"] / v["runs"], 1) if v["runs"] else 0}
                    for k, v in sorted(by_setting.items())]
    _write_tsv(out_dir / "cryptic_by_setting.tsv", setting_rows)

    print(f"\ndetections: {len(detections):,}   runs summarised: {len(paths):,}")
    print(f"{'setting':<18} {'runs':>6} {'cryptic':>8} {'%':>6} {'co-inf':>7} {'%':>6}")
    for r in setting_rows:
        print(f"{r['setting'][:18]:<18} {r['runs']:>6} {r['cryptic_runs']:>8} "
              f"{r['pct_cryptic']:>5}% {r['coinf_runs']:>7} {r['pct_coinf']:>5}%")
    for f in ("detections.tsv", "run_summary.tsv", "cryptic_by_setting.tsv"):
        print(f"  wrote {out_dir / f}")


if __name__ == "__main__":
    main()
