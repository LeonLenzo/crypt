#!/usr/bin/env python3
"""
filter_refs.py — apply BUSCO completeness thresholds to ref_candidates.tsv,
producing the final ref_screen.tsv used by build.py.

Rules:
  - For each accession: include if BUSCO complete_pct >= kingdom threshold.
  - For each taxid where no assembly passes: include the best available
    (highest complete_pct, or status=no_cds if nothing else) and flag it
    as busco_below_threshold so build.py / the analyst can see it.
  - Accessions with status=no_cds or busco_error are excluded unless they
    are the only option for a taxid.

Thresholds (can be overridden via CLI):
  --fungi-threshold    50.0   (low to accommodate rusts on fungi_odb10)
  --oomycete-threshold 65.0

Run from crypt/:
    python kraken/filter_refs.py [--fungi-threshold N] [--oomycete-threshold N]

Output: kraken/ref_screen.tsv
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

CANDIDATES_TSV = Path("kraken/ref_candidates.tsv")
SCORES_TSV     = Path("kraken/busco_scores.tsv")
OUT_TSV        = Path("kraken/ref_screen.tsv")

OUT_COLS = [
    "taxid", "organism_name", "kingdom", "source",
    "accession", "assembly_level", "release_date", "has_annotation",
    "country", "protein_coding_genes", "scaffold_n50_kb", "total_length_mb",
    "fasta_type", "busco_lineage", "busco_complete_pct", "busco_status",
    "selection_rank", "selection_reason",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates-tsv",    default=str(CANDIDATES_TSV))
    ap.add_argument("--scores-tsv",        default=str(SCORES_TSV))
    ap.add_argument("--out-tsv",           default=str(OUT_TSV))
    ap.add_argument("--fungi-threshold",   type=float, default=50.0)
    ap.add_argument("--oomycete-threshold", type=float, default=65.0)
    args = ap.parse_args()

    thresholds = {
        "fungal":   args.fungi_threshold,
        "oomycete": args.oomycete_threshold,
    }

    # Load BUSCO scores keyed by accession
    scores: dict[str, dict] = {}
    scores_path = Path(args.scores_tsv)
    if not scores_path.exists():
        print(f"ERROR: {scores_path} not found — run busco_screen.py first")
        raise SystemExit(1)
    with open(scores_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row["accession"]
            # Keep last entry per accession (append-only cache)
            scores[acc] = row
    print(f"BUSCO scores loaded: {len(scores):,} accessions")

    # Load candidates
    candidates: list[dict] = []
    with open(Path(args.candidates_tsv), newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            candidates.append(row)
    print(f"Candidates:          {len(candidates):,}")

    # Group candidates by taxid, annotate with BUSCO score
    by_taxid: dict[str, list[dict]] = defaultdict(list)
    n_no_score = 0
    for cand in candidates:
        acc     = cand["accession"]
        kingdom = cand["kingdom"]
        sc      = scores.get(acc, {})
        status  = sc.get("status", "not_run")
        try:
            pct = float(sc.get("complete_pct", "")) if sc.get("complete_pct") else None
        except ValueError:
            pct = None

        threshold = thresholds.get(kingdom, 50.0)
        if status == "pass" and pct is not None and pct >= threshold:
            busco_status = "pass"
        elif status in ("no_cds", "busco_error", "not_run"):
            busco_status = status
            if status == "not_run":
                n_no_score += 1
        else:
            busco_status = "fail"

        cand["busco_complete_pct"] = f"{pct:.1f}" if pct is not None else ""
        cand["busco_status"]       = busco_status
        by_taxid[cand["taxid"]].append(cand)

    if n_no_score:
        print(f"Warning: {n_no_score:,} candidates have no BUSCO score (not yet run)")

    # Select assemblies per taxid
    selected: list[dict] = []
    n_pass      = 0
    n_fallback  = 0
    n_no_option = 0

    for taxid, cands in by_taxid.items():
        passing = [c for c in cands if c["busco_status"] == "pass"]

        if passing:
            # Re-rank within the passing set (preserve original diversity order)
            for rank, c in enumerate(passing, 1):
                c["selection_rank"] = rank
            selected.extend(passing)
            n_pass += len(passing)
        else:
            # Fallback: best available — prefer highest complete_pct, then any cds
            def fallback_key(c):
                try:
                    pct = float(c["busco_complete_pct"]) if c["busco_complete_pct"] else -1
                except ValueError:
                    pct = -1
                has_cds = 1 if c["fasta_type"] == "cds" else 0
                return (has_cds, pct)

            best = sorted(cands, key=fallback_key, reverse=True)[0]
            best["selection_reason"] = (
                best.get("selection_reason", "") + "; busco_below_threshold"
            ).lstrip("; ")
            best["busco_status"] = "below_threshold_fallback"
            best["selection_rank"] = 1
            selected.append(best)
            if best["busco_complete_pct"]:
                n_fallback += 1
            else:
                n_no_option += 1

    # Sort for readability
    selected.sort(key=lambda r: (
        r["kingdom"], r["source"], r["organism_name"].lower(),
        int(r["selection_rank"])
    ))

    # Write output
    out_path = Path(args.out_tsv)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(selected)

    # Summary
    n_cds    = sum(1 for r in selected if r["fasta_type"] == "cds")
    n_genome = sum(1 for r in selected if r["fasta_type"] == "genome")
    n_fungi  = sum(1 for r in selected if r["kingdom"] == "fungal")
    n_oomy   = sum(1 for r in selected if r["kingdom"] == "oomycete")

    print(f"\n── Filter summary ───────────────────────────────────────────────")
    print(f"  Total selected:         {len(selected):>6,}")
    print(f"    BUSCO pass:           {n_pass:>6,}")
    print(f"    Fallback (scored):    {n_fallback:>6,}  (below threshold, only option)")
    print(f"    Fallback (no score):  {n_no_option:>6,}  (busco_error / no_cds)")
    print(f"  CDS FASTA:              {n_cds:>6,}")
    print(f"  Genomic FASTA:          {n_genome:>6,}")
    print(f"  Fungal:                 {n_fungi:>6,}")
    print(f"  Oomycete:               {n_oomy:>6,}")
    print(f"  Thresholds:  fungi ≥ {args.fungi_threshold}%,  "
          f"oomycete ≥ {args.oomycete_threshold}%")
    print(f"\nOutput: {out_path}")
    print(f"Next:   python kraken/build.py --genomes-dir <scratch>/cds_v2 "
          f"--db-dir <scratch>/db_v3")


if __name__ == "__main__":
    main()
