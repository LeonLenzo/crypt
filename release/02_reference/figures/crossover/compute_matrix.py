#!/usr/bin/env python3
# PROVENANCE for crossover_matrix.tsv (the crossover supplementary figure).
# Run on Setonix via compute_matrix.slurm (needs the CDS on scratch and py-numpy); the
# output was imported here as crossover_matrix.tsv. Not part of the build pipeline: it is a
# one-off measurement over the db1p CDS, kept beside the figure so the panel is reproducible.
"""Full pairwise shared-minimizer matrix for a taxon set (default: the 254 wheat-row taxa).

Same minimizer scheme as wp_shared_minimizers.py (canonical l=31 over k=35 windows, scrambled
order). Computes each species' minimizer set (union over its assemblies), then the full NxN
shared-minimizer matrix via parallel pairwise intersection over the upper triangle.

  python3 wp_matrix.py <meta.tsv> <out_matrix.tsv>
  meta.tsv: taxid<TAB>phylum<TAB>class<TAB>order<TAB>family<TAB>genus<TAB>name  (one per taxon)
"""
import sys, csv, os
from multiprocessing import Pool
import numpy as np

META, OUT = sys.argv[1], sys.argv[2]
ROOT = "/scratch/pawsey1168/llenzo/crypt"
CDS = os.path.join(ROOT, "release/02_reference/data/raw/cds")
PLAN = os.path.join(ROOT, "release/02_reference/data/curation/build_plan.tsv")
NPROC = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))

K, L = 35, 31
W = K - L + 1
POW = (np.uint64(4) ** np.arange(L - 1, -1, -1, dtype=np.uint64))
RPOW = POW[::-1].copy()
SCRAMBLE = np.uint64(0x9E3779B97F4A7C15)
U64MAX = np.uint64(0xFFFFFFFFFFFFFFFF)
CODE = np.full(256, 255, np.uint8)
for b, c in zip(b"ACGT", (0, 1, 2, 3)):
    CODE[b] = c


def seq_minimizers(seq):
    codes = CODE[np.frombuffer(seq.encode().upper(), np.uint8)]
    if codes.size < K:
        return np.empty(0, np.uint64)
    valid = codes != 255
    c = codes.astype(np.uint64); c[~valid] = 0
    comp = np.uint64(3) - c
    fwd = np.lib.stride_tricks.sliding_window_view(c, L) @ POW
    rc = np.lib.stride_tricks.sliding_window_view(comp, L) @ RPOW
    canon = np.minimum(fwd, rc)
    vwin = np.lib.stride_tricks.sliding_window_view(valid, L).all(axis=1)
    key = (canon * SCRAMBLE) & U64MAX
    key[~vwin] = U64MAX
    if key.size < W:
        return np.empty(0, np.uint64)
    kwin = np.lib.stride_tricks.sliding_window_view(key, W)
    cwin = np.lib.stride_tricks.sliding_window_view(canon, W)
    arg = kwin.argmin(axis=1)
    chosen = kwin[np.arange(kwin.shape[0]), arg]
    mins = cwin[np.arange(cwin.shape[0]), arg]
    return np.unique(mins[chosen != U64MAX])


def assembly_minimizers(acc):
    path = os.path.join(CDS, acc, "cds_from_genomic.fna")
    if not os.path.exists(path):
        return np.empty(0, np.uint64)
    parts, chunks = [], []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if chunks:
                    parts.append(seq_minimizers("".join(chunks)))
                chunks = []
            else:
                chunks.append(line.strip())
    if chunks:
        parts.append(seq_minimizers("".join(chunks)))
    return np.unique(np.concatenate(parts)) if parts else np.empty(0, np.uint64)


SPECIES = []        # global for phase-2 fork: list of uint64 arrays, indexed by position


def pair_block(args):
    i, js = args
    a = SPECIES[i]
    return i, [(j, int(np.intersect1d(a, SPECIES[j], assume_unique=True).size)) for j in js]


def main():
    meta = [r for r in csv.reader(open(META), delimiter="\t")]
    taxids = [r[0] for r in meta]
    # order by taxonomy so the heatmap is block-structured
    order = sorted(range(len(meta)), key=lambda k: tuple(meta[k][1:7]))
    meta = [meta[k] for k in order]
    taxids = [m[0] for m in meta]
    want = set(taxids)

    acc_by_tax = {}
    with open(PLAN) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["species_taxid"] in want:
                acc_by_tax.setdefault(row["species_taxid"], []).append(row["accession"])

    all_acc = [(t, a) for t in taxids for a in acc_by_tax.get(t, [])]
    print(f"{len(taxids)} taxa, {len(all_acc)} assemblies, {NPROC} procs", flush=True)
    with Pool(NPROC) as pool:
        res = pool.map(assembly_minimizers, [a for _, a in all_acc], chunksize=1)
    per = {}
    for (t, _), arr in zip(all_acc, res):
        per.setdefault(t, []).append(arr)

    global SPECIES
    SPECIES = [np.unique(np.concatenate(per[t])) if t in per else np.empty(0, np.uint64)
               for t in taxids]
    sizes = [s.size for s in SPECIES]
    print(f"minimizer sets built; total distinct (sum) {sum(sizes):,}", flush=True)

    n = len(taxids)
    mat = np.zeros((n, n), np.int64)
    for i in range(n):
        mat[i, i] = sizes[i]
    tasks = [(i, list(range(i + 1, n))) for i in range(n - 1)]
    with Pool(NPROC) as pool:
        for i, out in pool.map(pair_block, tasks, chunksize=1):
            for j, v in out:
                mat[i, j] = v
                mat[j, i] = v
    print("matrix complete", flush=True)

    with open(OUT, "w") as fh:
        fh.write("taxid\tphylum\tclass\torder\tfamily\tgenus\tname\tdistinct\t"
                 + "\t".join(taxids) + "\n")
        for i, m in enumerate(meta):
            fh.write("\t".join(m) + f"\t{sizes[i]}\t"
                     + "\t".join(str(int(x)) for x in mat[i]) + "\n")
    # crossover summary
    off = mat.copy(); np.fill_diagonal(off, 0)
    print(f"largest off-diagonal shared: {off.max():,}")
    print(f"matrix -> {OUT}")


if __name__ == "__main__":
    main()
