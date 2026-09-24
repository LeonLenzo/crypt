#!/usr/bin/env python3
"""
kraken_db_saturation.py — measure how many assemblies per species a Kraken2 build
actually needs, by k-mer accumulation.

Answers two questions that pull in opposite directions:

  1. DETECTION. How fast does a species' distinct 31-mer content saturate as
     assemblies are added? Past the plateau, another assembly adds no new k-mers
     for a read to match, so it buys nothing.
  2. DISCRIMINATION. As a species' k-mer set grows, how much of it is shared with
     its congeners in the same build? Shared k-mers are pushed to genus by
     Kraken2's LCA, so this fraction is the cost of the extra assemblies. The
     useful number is where the shared fraction starts growing faster than the
     unique fraction.

Kraken2 stores each distinct minimizer once, so assembly COUNT does not weight a
taxon. What extra assemblies buy is pan-genome breadth, which is what this measures.

Method: FracMinHash sketching (the sourmash estimator). Every canonical 31-mer is
hashed with splitmix64 and kept when hash < 2^64/scaled, so the sketch is a uniform
1/scaled sample of the k-mer set. Sketches of different assemblies sample the SAME
subspace, so union and intersection of sketches estimate union and intersection of
the full k-mer sets. Cardinality is |sketch| x scaled, with relative error ~1/sqrt(|sketch|).
At scaled=1000 a 35 Mbp CDS set gives ~35k hashes, so ~0.5%. Exact counting needs
KMC or jellyfish, neither of which is on Setonix.

Two stages, both resumable:
  sketch  one .npy per accession under data/sketches/ (the expensive stage)
  curve   accumulation curves + congener overlap from the cached sketches (seconds)

Run on Setonix: sbatch kraken/slurm/kraken_db_saturation.slurm
Run from crypt/:
    python kraken/db/kraken_db_saturation.py \\
        [--genomes-dir /scratch/pawsey1168/llenzo/crypt/kraken/output/db/search/data/cds/pathogen] \\
        [--stage {sketch,curve,all}] [--workers 16] [--scaled 1000] [--k 31]
        [--min-assemblies 4] [--permutations 25]

Outputs:
    kraken/output/db/saturation/data/sketches/{accession}.npy   (gitignored — cache)
    kraken/output/db/saturation/data/sketch_index.tsv           (tracked — what was sketched)
    kraken/output/db/saturation/data/saturation_curve.tsv       (tracked — the curves)
    kraken/output/db/saturation/data/saturation_summary.tsv     (tracked — one row per taxon)
"""

import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest

BUSCO_TSV = Path("kraken/output/db/busco/data/busco_scores.tsv")
DEFAULT_GENOMES_DIR = Path("kraken/output/db/search/data/cds/pathogen")

OUT_DIR      = Path("kraken/output/db/saturation")
DATA_DIR     = OUT_DIR / "data"
SKETCH_DIR   = DATA_DIR / "sketches"
INDEX_TSV    = DATA_DIR / "sketch_index.tsv"
CURVE_TSV    = DATA_DIR / "saturation_curve.tsv"
SUMMARY_TSV  = DATA_DIR / "saturation_summary.tsv"
BIAS_TSV     = DATA_DIR / "depth_bias.tsv"

# splitmix64 finaliser, as numpy uint64 so the arithmetic wraps rather than promoting
_GOLD = np.uint64(0x9E3779B97F4A7C15)
_C1   = np.uint64(0xBF58476D1CE4E5B9)
_C2   = np.uint64(0x94D049BB133111EB)
_S30, _S27, _S31, _TWO = np.uint64(30), np.uint64(27), np.uint64(31), np.uint64(2)
_THREE = np.uint64(3)

# byte -> 2-bit code, 255 for anything that is not an unambiguous base
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _c in enumerate(b"ACGT"):
    _LUT[_c] = _i
    _LUT[_c + 32] = _i


def _ts():
    return time.strftime("%H:%M:%S")


def _splitmix64(x):
    with np.errstate(over="ignore"):
        z = x + _GOLD
        z = (z ^ (z >> _S30)) * _C1
        z = (z ^ (z >> _S27)) * _C2
        return z ^ (z >> _S31)


def _fasta_codes(path):
    """Whole FASTA as 2-bit codes, headers and ambiguity codes becoming 255 breaks."""
    parts = []
    with open(path, "rb") as fh:
        for line in fh:
            parts.append(b"\n" if line[:1] == b">" else line.strip())
    return _LUT[np.frombuffer(b"".join(parts), dtype=np.uint8)]


def sketch_fasta(path, k=31, scaled=1000, chunk=4_000_000):
    """FracMinHash sketch: sorted unique uint64 hashes of canonical k-mers."""
    max_hash = np.uint64((1 << 64) // scaled)
    codes = _fasta_codes(path)
    valid = (codes != 255).astype(np.int8)
    # run boundaries of unambiguous sequence
    edges = np.flatnonzero(np.diff(np.concatenate(([0], valid, [0]))))
    kept = []
    for start, end in zip(edges[0::2], edges[1::2]):
        if end - start < k:
            continue
        run = codes[start:end].astype(np.uint64)
        n_windows = (end - start) - k + 1
        for off in range(0, n_windows, chunk):
            m = min(chunk, n_windows - off)
            sub = run[off:off + m + k - 1]
            fwd = np.zeros(m, dtype=np.uint64)
            for j in range(k):
                fwd = (fwd << _TWO) | sub[j:j + m]
            comp = _THREE - sub
            rev = np.zeros(m, dtype=np.uint64)
            for j in range(k):
                rev = (rev << _TWO) | comp[k - 1 - j:k - 1 - j + m]
            h = _splitmix64(np.minimum(fwd, rev))
            sel = h[h < max_hash]
            if sel.size:
                kept.append(sel)
    if not kept:
        return np.empty(0, dtype=np.uint64)
    return np.unique(np.concatenate(kept))


def _cds_path(genomes_dir, accession):
    p = genomes_dir / accession / "ncbi_dataset" / "data" / accession / "cds_from_genomic.fna"
    if p.exists():
        return p
    hits = list((genomes_dir / accession).rglob("cds_from_genomic.fna"))
    return hits[0] if hits else None


def _sketch_one(args):
    accession, cds, out_npy, k, scaled = args
    try:
        t0 = time.time()
        sk = sketch_fasta(Path(cds), k=k, scaled=scaled)
        np.save(out_npy, sk)
        return accession, len(sk), time.time() - t0, ""
    except Exception as exc:                                   # noqa: BLE001
        return accession, 0, 0.0, "{}: {}".format(type(exc).__name__, exc)


def load_selected():
    """Selected assemblies from busco_scores.tsv, grouped by taxon, best-first."""
    with open(BUSCO_TSV, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["selected"].strip().lower() in ("true", "1", "yes")]

    def quality(r):
        def num(key):
            try:
                return float(r.get(key) or 0)
            except ValueError:
                return 0.0
        # BUSCO completeness first, then contiguity; both descending
        return (-num("complete_pct"), -num("scaffold_n50_kb"), r["accession"])

    by_taxon = defaultdict(list)
    for r in rows:
        by_taxon[(r["taxid"], r["organism_name"])].append(r)
    for key in by_taxon:
        by_taxon[key].sort(key=quality)
    return by_taxon


def run_sketch(by_taxon, genomes_dir, k, scaled, workers):
    SKETCH_DIR.mkdir(parents=True, exist_ok=True)
    jobs, skipped, missing = [], 0, []
    for (taxid, name), rows in sorted(by_taxon.items()):
        for r in rows:
            acc = r["accession"]
            out = SKETCH_DIR / "{}.npy".format(acc)
            if out.exists() and out.stat().st_size > 0:
                skipped += 1
                continue
            cds = _cds_path(genomes_dir, acc)
            if cds is None:
                missing.append(acc)
                continue
            jobs.append((acc, str(cds), str(out), k, scaled))

    print("[{}] sketch: {} to do, {} cached, {} missing CDS".format(
        _ts(), len(jobs), skipped, len(missing)), flush=True)
    if missing:
        print("  missing: {}".format(", ".join(missing[:10])), flush=True)

    done = failed = 0
    if jobs:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(_sketch_one, j): j[0] for j in jobs}
            for fut in as_completed(futs):
                acc, n, secs, err = fut.result()
                done += 1
                if err:
                    failed += 1
                    print("[{}] {:>4}/{}  {}  FAILED  {}".format(
                        _ts(), done, len(jobs), acc, err), flush=True)
                elif done % 25 == 0 or done == len(jobs):
                    print("[{}] {:>4}/{}  {}  {:,} hashes  {:.1f}s".format(
                        _ts(), done, len(jobs), acc, n, secs), flush=True)
    print("[{}] sketch done: {} written, {} failed".format(_ts(), done - failed, failed),
          flush=True)

    with open(INDEX_TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["taxid", "organism_name", "accession", "complete_pct",
                    "scaffold_n50_kb", "n_hashes", "est_distinct_kmers"])
        for (taxid, name), rows in sorted(by_taxon.items()):
            for r in rows:
                p = SKETCH_DIR / "{}.npy".format(r["accession"])
                n = len(np.load(p)) if p.exists() else 0
                w.writerow([taxid, name, r["accession"], r.get("complete_pct", ""),
                            r.get("scaffold_n50_kb", ""), n, n * scaled])
    print("[{}] wrote {}".format(_ts(), INDEX_TSV), flush=True)


def _load_sketches(rows):
    out = []
    for r in rows:
        p = SKETCH_DIR / "{}.npy".format(r["accession"])
        if p.exists():
            out.append((r["accession"], np.load(p)))
    return out


def run_curve(by_taxon, scaled, min_assemblies, permutations, seed=1168):
    # congener pool: union of every sketch of every OTHER species in the same genus
    genus_members = defaultdict(list)
    for key in by_taxon:
        genus_members[key[1].split()[0]].append(key)

    targets = {k: v for k, v in by_taxon.items() if len(v) >= min_assemblies}
    print("[{}] curve: {} taxa with >= {} assemblies".format(
        _ts(), len(targets), min_assemblies), flush=True)

    rng = np.random.default_rng(seed)
    curve_rows, summary_rows = [], []

    for (taxid, name), rows in sorted(targets.items(), key=lambda kv: -len(kv[1])):
        sketches = _load_sketches(rows)
        if len(sketches) < min_assemblies:
            continue
        genus = name.split()[0]

        pool = [s for key in genus_members[genus] if key != (taxid, name)
                for _, s in _load_sketches(by_taxon[key])]
        congeners = np.unique(np.concatenate(pool)) if pool else np.empty(0, np.uint64)
        n_cong_taxa = len([k for k in genus_members[genus] if k != (taxid, name)])

        arrays = [s for _, s in sketches]
        n = len(arrays)

        def accumulate(order):
            cum = np.empty(0, dtype=np.uint64)
            sizes, shared = [], []
            for idx in order:
                cum = np.union1d(cum, arrays[idx])
                sizes.append(len(cum))
                shared.append(len(np.intersect1d(cum, congeners, assume_unique=True))
                              if congeners.size else 0)
            return sizes, shared

        q_sizes, q_shared = accumulate(range(n))

        perm_sizes = np.zeros(n)
        perm_shared = np.zeros(n)
        reps = min(permutations, 200)
        for _ in range(reps):
            s, h = accumulate(rng.permutation(n))
            perm_sizes += np.asarray(s, dtype=float)
            perm_shared += np.asarray(h, dtype=float)
        perm_sizes /= reps
        perm_shared /= reps

        for i in range(n):
            prev_q = q_sizes[i - 1] if i else 0
            prev_p = perm_sizes[i - 1] if i else 0.0
            curve_rows.append({
                "taxid": taxid, "organism_name": name, "genus": genus,
                "n_assemblies": i + 1,
                "quality_kmers": q_sizes[i] * scaled,
                "quality_new": (q_sizes[i] - prev_q) * scaled,
                "quality_frac_new": round((q_sizes[i] - prev_q) / q_sizes[i], 6) if q_sizes[i] else 0,
                "quality_shared_congener": q_shared[i] * scaled,
                "quality_frac_shared": round(q_shared[i] / q_sizes[i], 6) if q_sizes[i] else 0,
                "mean_kmers": int(round(perm_sizes[i] * scaled)),
                "mean_new": int(round((perm_sizes[i] - prev_p) * scaled)),
                "mean_frac_new": round((perm_sizes[i] - prev_p) / perm_sizes[i], 6) if perm_sizes[i] else 0,
                "mean_frac_shared": round(perm_shared[i] / perm_sizes[i], 6) if perm_sizes[i] else 0,
                "n_congener_taxa": n_cong_taxa,
                "permutations": reps,
            })

        # saturation point: first n where the mean curve gains under 1% / 2% / 5%
        def first_below(thresh):
            for i in range(1, n):
                if perm_sizes[i] and (perm_sizes[i] - perm_sizes[i - 1]) / perm_sizes[i] < thresh:
                    return i + 1
            return None

        summary_rows.append({
            "taxid": taxid, "organism_name": name, "genus": genus,
            "n_assemblies": n,
            "kmers_1_assembly": int(round(perm_sizes[0] * scaled)),
            "kmers_all_assemblies": int(round(perm_sizes[-1] * scaled)),
            "pangenome_ratio": round(perm_sizes[-1] / perm_sizes[0], 4) if perm_sizes[0] else 0,
            "n_for_95pct": next((i + 1 for i in range(n)
                                 if perm_sizes[i] >= 0.95 * perm_sizes[-1]), n),
            "n_for_99pct": next((i + 1 for i in range(n)
                                 if perm_sizes[i] >= 0.99 * perm_sizes[-1]), n),
            "saturate_5pct": first_below(0.05),
            "saturate_2pct": first_below(0.02),
            "saturate_1pct": first_below(0.01),
            "n_congener_taxa": n_cong_taxa,
            "frac_shared_1_assembly": round(perm_shared[0] / perm_sizes[0], 6) if perm_sizes[0] else 0,
            "frac_shared_all": round(perm_shared[-1] / perm_sizes[-1], 6) if perm_sizes[-1] else 0,
        })
        print("[{}] {:<44} n={:>3}  pan/single={:.2f}x  95% at n={}  shared {:.1%} -> {:.1%}".format(
            _ts(), name[:44], n, summary_rows[-1]["pangenome_ratio"],
            summary_rows[-1]["n_for_95pct"], summary_rows[-1]["frac_shared_1_assembly"],
            summary_rows[-1]["frac_shared_all"]), flush=True)

    for path, rows in ((CURVE_TSV, curve_rows), (SUMMARY_TSV, summary_rows)):
        if not rows:
            continue
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print("[{}] wrote {} ({} rows)".format(_ts(), path, len(rows)), flush=True)


def run_bias(scaled, min_deep=8, min_shallow=4):
    """Does uneven sampling depth between congeners misattribute reads?

    A k-mer truly present in both species A and B is labelled by Kraken2 from what
    the build OBSERVED, not from what is true. If A is sampled deeply and B by one
    assembly, a k-mer in both pangenomes but absent from B's single representative
    is seen only in A, so it is stored as A-specific and a read from B carrying it
    scores for A. That is the bias, and it is driven by the DIFFERENCE in depth,
    not by A's assembly count.

    For every congener pair (A deep, B shallow) this partitions B's pangenome
    against A at each depth of A:
      ambiguous    in A(d) and in B(1)         -> LCA sends it to genus, honest
      misassigned  in A(d), in B_full, not B(1) -> stored as A, pulls B's reads to A
    then repeats with B represented by 4 assemblies to show what closing the gap does.
    """
    by_taxon = load_selected()
    genus = defaultdict(list)
    for key in by_taxon:
        genus[key[1].split()[0]].append(key)

    def union(rows, upto=None):
        arrs = []
        for r in rows[:upto]:
            p = SKETCH_DIR / "{}.npy".format(r["accession"])
            if p.exists():
                arrs.append(np.load(p))
        return np.unique(np.concatenate(arrs)) if arrs else np.empty(0, np.uint64)

    out = []
    for g, keys in sorted(genus.items()):
        deep = [k for k in keys if len(by_taxon[k]) >= min_deep]
        shallow = [k for k in keys if len(by_taxon[k]) >= min_shallow]
        for a in deep:
            for b in shallow:
                if a == b:
                    continue
                arows, brows = by_taxon[a], by_taxon[b]
                b_full = union(brows)
                if b_full.size == 0:
                    continue
                b_at = {n: union(brows, n) for n in (1, 4) if len(brows) >= n}
                for d in range(1, len(arows) + 1):
                    a_d = union(arows, d)
                    hit = np.intersect1d(a_d, b_full, assume_unique=True)
                    row = {"genus": g, "deep_species": a[1], "shallow_species": b[1],
                           "deep_depth": d, "shallow_pangenome_kmers": b_full.size * scaled,
                           "overlap_kmers": hit.size * scaled}
                    for n, bn in b_at.items():
                        amb = np.intersect1d(hit, bn, assume_unique=True).size
                        row["ambiguous_at_b{}".format(n)] = amb * scaled
                        row["misassigned_at_b{}".format(n)] = (hit.size - amb) * scaled
                        row["frac_misassigned_at_b{}".format(n)] = (
                            round((hit.size - amb) / b_full.size, 6))
                    out.append(row)
                print("[{}] {:<26} vs {:<26} depths 1-{}".format(
                    _ts(), a[1][:26], b[1][:26], len(arows)), flush=True)

    if out:
        cols = sorted({k for r in out for k in r}, key=lambda c: (c[0].isdigit(), c))
        order = ["genus", "deep_species", "shallow_species", "deep_depth",
                 "shallow_pangenome_kmers", "overlap_kmers"]
        cols = order + [c for c in cols if c not in order]
        with open(BIAS_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(out)
        print("[{}] wrote {} ({} rows)".format(_ts(), BIAS_TSV, len(out)), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes-dir", default=str(DEFAULT_GENOMES_DIR))
    ap.add_argument("--stage", choices=("sketch", "curve", "bias", "all"), default="all")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--scaled", type=int, default=1000,
                    help="keep 1 hash in SCALED; lower is more accurate and slower to combine")
    ap.add_argument("--min-assemblies", type=int, default=4,
                    help="only curve taxa with at least this many selected assemblies")
    ap.add_argument("--permutations", type=int, default=25,
                    help="random assembly orders averaged for the unbiased curve")
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_db_saturation.log")
    link_latest(logs_base, log_dir / "kraken_db_saturation.log")
    sys.stdout = log

    try:
        print("[{}] k={} scaled={} stage={}".format(_ts(), args.k, args.scaled, args.stage),
              flush=True)
        by_taxon = load_selected()
        print("[{}] {} selected assemblies across {} taxa".format(
            _ts(), sum(len(v) for v in by_taxon.values()), len(by_taxon)), flush=True)

        if args.stage in ("sketch", "all"):
            run_sketch(by_taxon, Path(args.genomes_dir), args.k, args.scaled, args.workers)
        if args.stage in ("curve", "all"):
            run_curve(by_taxon, args.scaled, args.min_assemblies, args.permutations)
        if args.stage in ("bias", "all"):
            run_bias(args.scaled)
    finally:
        log.close()


if __name__ == "__main__":
    main()
