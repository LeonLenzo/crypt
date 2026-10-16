#!/usr/bin/env python3
"""compare.py — score predicted CDS against the k-mer space that actually matters.

Final step of the predict evaluation. Reads predicted CDS, compares it to ground truth
(Arm A) or to db_v4's existing k-mer space (Arm B), and writes `data/gold/predict_benchmark.tsv`.

## Why containment and not F1

Vipsania reports Fungi locus-level F1 0.660. That asks "is this gene model correct", which is
the wrong question for a Kraken2 CDS reference: we do not need correct gene models, we need
CDS-like sequence from the right organism. A misplaced exon boundary still yields real
sequence from that species, and its k-mers still point at the right taxon.

The probe made the gap concrete. On *P. nodorum* GCA_016801405.1, vipsania pretrained produced
**13,634 CDS against the true 17,764 (77%), but 19.29 Mb of sequence against the true
20.55 Mb (94%)**. By gene count it looks mediocre; by sequence bulk it looks excellent. Neither
number says whether it is the SAME sequence, which is the only thing a classifier cares about.
That is what this script measures.

Two directions, two distinct failure modes:

  * **recall = |T ∩ P| / |T|** — the fraction of TRUE CDS k-mer space the prediction recovers.
    Too low and the species cannot be detected even when it is genuinely present.
  * **precision = |T ∩ P| / |P|** — the fraction of PREDICTED k-mer space that is real CDS.
    Invented sequence inflates the taxon's minimizer count, and that count is the DENOMINATOR
    every coverage figure in `05_detect` divides by, so over-prediction deflates the species'
    own coverage and weakens it as a competitor.

Note this is narrower than `kraken_annotate.py`'s claim that "an invented gene corrupts every
species' apportionment". K-mers shared between taxa resolve to their LCA whether or not they
were invented, so what over-prediction corrupts is the predicted species' own denominator, not
its neighbours'. That is a real cost and a measurable one, but it is not contagious.

## Arm B has no truth, so it measures novelty instead

For an unannotated species nothing can be scored for correctness. What decides a db_v5 is
**how much of the predicted k-mer space is NEW to db_v4**: a prediction 95% contained in
species the database already holds adds no competitor value however accurate it is.
`--db-sketches` points at the db_v4 sketch cache and the script reports
`novelty = 1 - |P ∩ DB| / |P|`.

## Method

FracMinHash, k=31, scaled=200, imported from `utilities/saturation/saturation.py` rather than
reimplemented. One sketching implementation, already validated against a pure-Python brute
force over soft-masked sequence, ambiguity codes, chunk boundaries, short contigs and
palindromic k-mers. It also means these numbers are directly comparable to the existing
saturation and depth-bias tables, which were computed at exactly these parameters.

    python 02_classify/01_refs/predict/compare.py --tool vipsania_pretrained \\
        --predicted-dir data/bronze/predicted/vipsania_pretrained
"""

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np

MODULE = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import (CDS_PATHOGEN, GOLD, LOGS, PREDICTED, PREDICT_METRICS,  # noqa: E402
                     PREDICT_SAMPLE, SKETCHES)
from _util import _Tee, atomic_write, link_latest, make_log_dir  # noqa: E402

K, SCALED = 31, 200
DERIVED_SUFFIXES = ("_clean.fna", "_combined.fna", ".tagged.fna")


def _load_sketcher():
    """Import sketch_fasta from the saturation utility; one implementation, already validated."""
    path = MODULE / "utilities" / "saturation" / "saturation.py"
    spec = importlib.util.spec_from_file_location("_sat", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.sketch_fasta


def truth_cds(acc: str, cds_dir: Path):
    """The NCBI CDS for an accession — the ground truth for Arm A."""
    d = cds_dir / acc
    if not d.is_dir():
        return None
    c = [f for f in d.glob("**/*.fna")
         if f.is_file() and not f.name.endswith(DERIVED_SUFFIXES)
         and f.stat().st_size >= 100_000]
    return max(c, key=lambda f: f.stat().st_size) if c else None


def fasta_stats(path: Path) -> tuple:
    """(n_sequences, n_bases) in one pass, without holding the file in memory."""
    n_seq = n_base = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n_seq += 1
            else:
                n_base += len(line.strip())
    return n_seq, n_base


def containment(a: np.ndarray, b: np.ndarray) -> float:
    """|a ∩ b| / |a|. Both sketches sample the SAME hash subspace, so this estimates the
    containment of the full k-mer sets."""
    if a.size == 0:
        return float("nan")
    return float(np.intersect1d(a, b, assume_unique=True).size) / float(a.size)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tool", required=True,
                    help="label for this predictor run, e.g. vipsania_pretrained, "
                         "vipsania_finetuned, miniprot")
    ap.add_argument("--predicted-dir", default=None,
                    help=f"tree of <accession>/coding.fa (default {PREDICTED}/<tool>)")
    ap.add_argument("--coding-name", default="coding.fa",
                    help="filename of the predicted CDS inside each accession dir")
    ap.add_argument("--sample", default=None, help=f"default {PREDICT_SAMPLE}")
    ap.add_argument("--cds-dir", default=None, help=f"truth CDS tree (default {CDS_PATHOGEN})")
    ap.add_argument("--db-sketches", default=None,
                    help=f"db_v4 sketch cache for the Arm B novelty measure "
                         f"(default {SKETCHES}); omit to skip novelty")
    ap.add_argument("--single", nargs=2, metavar=("ACCESSION", "CODING_FA"),
                    help="score ONE accession directly, for a probe result that is not yet "
                         "laid out as a tree")
    ap.add_argument("--out", default=None, help=f"default {PREDICT_METRICS}")
    args = ap.parse_args()

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / f"predict_compare_{args.tool}.log")
    link_latest(LOGS, log_dir)
    try:
        sketch_fasta = _load_sketcher()
        cds_dir = Path(args.cds_dir) if args.cds_dir else CDS_PATHOGEN
        pdir = Path(args.predicted_dir) if args.predicted_dir else PREDICTED / args.tool

        if args.single:
            acc, coding = args.single
            jobs = [{"accession": acc, "arm": "annotated", "organism_name": "(single)",
                     "assembly_level": "", "focal": "", "_coding": Path(coding)}]
        else:
            sample_path = Path(args.sample) if args.sample else PREDICT_SAMPLE
            if not sample_path.exists():
                sys.exit(f"no sample at {sample_path} — run sample.py first")
            jobs = []
            for r in csv.DictReader(open(sample_path, newline=""), delimiter="\t"):
                r = dict(r)
                r["_coding"] = pdir / r["accession"] / args.coding_name
                jobs.append(r)

        db_sketch = None
        if args.db_sketches:
            sdir = Path(args.db_sketches)
            files = sorted(sdir.glob("*.npy"))
            if files:
                print(f"loading db_v4 k-mer space from {len(files):,} cached sketches "
                      f"({sdir})", flush=True)
                parts = [np.load(f) for f in files]
                db_sketch = np.unique(np.concatenate(parts)) if parts else None
                print(f"  db_v4 union: {db_sketch.size:,} hashes "
                      f"(~{db_sketch.size * SCALED / 1e9:.2f} B k-mers)")
            else:
                print(f"no sketches under {sdir}; novelty will be blank")

        rows, skipped = [], 0
        for j in jobs:
            coding = j["_coding"]
            if not coding.is_file() or coding.stat().st_size == 0:
                skipped += 1
                continue
            p_seq, p_base = fasta_stats(coding)
            P = sketch_fasta(coding, k=K, scaled=SCALED)
            row = {
                "tool": args.tool, "arm": j["arm"], "accession": j["accession"],
                "organism_name": j["organism_name"],
                "assembly_level": j.get("assembly_level", ""), "focal": j.get("focal", ""),
                "pred_seqs": p_seq, "pred_bases": p_base,
                "pred_kmers_est": int(P.size * SCALED),
                "true_seqs": "", "true_bases": "", "true_kmers_est": "",
                "recall": "", "precision": "", "jaccard": "", "novelty": "",
            }
            truth = truth_cds(j["accession"], cds_dir) if j["arm"] == "annotated" else None
            if truth:
                t_seq, t_base = fasta_stats(truth)
                T = sketch_fasta(truth, k=K, scaled=SCALED)
                inter = np.intersect1d(T, P, assume_unique=True).size
                union = T.size + P.size - inter
                row.update({
                    "true_seqs": t_seq, "true_bases": t_base,
                    "true_kmers_est": int(T.size * SCALED),
                    "recall": f"{inter / T.size:.4f}" if T.size else "",
                    "precision": f"{inter / P.size:.4f}" if P.size else "",
                    "jaccard": f"{inter / union:.4f}" if union else "",
                })
            if db_sketch is not None and P.size:
                row["novelty"] = f"{1.0 - containment(P, db_sketch):.4f}"
            rows.append(row)
            tag = (f"recall {row['recall']} precision {row['precision']}"
                   if row["recall"] else f"novelty {row['novelty'] or '-'}")
            print(f"  {j['accession']:18s} {j['organism_name'][:28]:28s} "
                  f"pred {p_seq:6,d} seq / {p_base/1e6:5.2f} Mb   {tag}", flush=True)

        if not rows:
            sys.exit(f"no predicted CDS found under {pdir} "
                     f"(looked for <accession>/{args.coding_name}); {skipped} missing")

        out = Path(args.out) if args.out else PREDICT_METRICS
        cols = list(rows[0].keys())
        GOLD.mkdir(parents=True, exist_ok=True)
        existing = []
        if out.exists():          # accumulate across tools rather than overwrite
            with open(out, newline="") as fh:
                existing = [r for r in csv.DictReader(fh, delimiter="\t")
                            if r.get("tool") != args.tool]
        with atomic_write(out) as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(existing + rows)
        print(f"\nmetrics → {out}   ({len(rows)} scored, {skipped} missing)")

        scored = [r for r in rows if r["recall"]]
        if scored:
            def med(key):
                v = sorted(float(r[key]) for r in scored if r[key])
                return v[len(v) // 2] if v else float("nan")
            print("\n── Arm A (ground truth exists) ──────────────────────────")
            print(f"  assemblies scored : {len(scored)}")
            print(f"  median RECALL     : {med('recall'):.3f}   "
                  f"(of TRUE CDS k-mer space recovered)")
            print(f"  median PRECISION  : {med('precision'):.3f}   "
                  f"(of PREDICTED k-mer space that is real CDS)")
            print(f"  median Jaccard    : {med('jaccard'):.3f}")
        novel = [r for r in rows if r["novelty"]]
        if novel:
            v = sorted(float(r["novelty"]) for r in novel)
            print("\n── novelty against db_v4 ────────────────────────────────")
            print(f"  median NOVELTY    : {v[len(v)//2]:.3f}   "
                  f"(of predicted k-mer space NOT already in db_v4)")
    finally:
        log.close()


if __name__ == "__main__":
    main()
