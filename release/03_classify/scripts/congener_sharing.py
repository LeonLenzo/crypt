#!/usr/bin/env python3
"""Pairwise minimizer sharing WITHIN each genus that the cohort calls.

    python3 release/03_classify/scripts/congener_sharing.py <meta.tsv> <out.tsv>
    (run on Setonix via scripts/slurm/congener_sharing.slurm; needs the CDS on scratch)

Writes one row per ordered congener pair: what fraction of B's own minimizers are also in A.
That fraction is the spillover coefficient. If A is present at coverage c and B is absent,
B is expected to score about share(B|A) * c purely from k-mers the two share, because a read
from A carrying a shared minimizer is indistinguishable from a read from B.

## Why this exists

The coverage floor is anchored on the off-host null, and that null cannot see this failure
mode. Measured on db2: four of the null's five genera hold one or two species and have blur
at or below 0.014, and they contribute ZERO rows where a congener is also present. The only
genus that does contribute is Puccinia, 1,289 rows, and Puccinia is the cleanest genus in the
database (blur 0.010, species separable across a 16-point ANI gap). So the null bounds
congener spillover for the best case in the reference and says nothing about Zymoseptoria
(blur 0.141), Fusarium (0.150), Alternaria (0.202) or Pyrenophora (0.254).

Zymoseptoria brevis is what exposed it: 24 calls, every one alongside Z. tritici, distinct
k-mers correlating with Z. tritici's at r = 0.868 and a near-constant ratio of 0.0205, and a
median coverage of 0.0245 against a floor of 0.02. It cleared the floor by a fifth.

## Only within-genus pairs, and that is not a shortcut

Spillover that matters comes from a close relative: 02_reference's crossover measurement
found sharing block-diagonal by genus, mean 0.015 within a genus against 0.0001 across. The
cross-genus term is two orders down and cannot lift a call off the floor. Restricting to
within-genus also makes this CHEAPER than the 254-taxon crossover run, 15,526 pairs against
32,131, while covering 487 species rather than 254.

## Same minimizer scheme as the crossover figure, deliberately

Canonical l=31 minimizers over k=35 windows with the same scramble, union over a species'
assemblies, read from the same CDS tree. The numbers are therefore comparable with
crossover_matrix.tsv rather than being a second, slightly different measurement of the same
thing. The code below is taken from 02_reference/figures/crossover/compute_matrix.py.
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

SPECIES = []        # global for the phase-2 fork: list of uint64 arrays, indexed by position


def pair_block(args):
    """Intersect one species against its remaining congeners."""
    i, js = args
    a = SPECIES[i]
    return i, [(j, int(np.intersect1d(a, SPECIES[j], assume_unique=True).size)) for j in js]


def main():
    meta = [r for r in csv.reader(open(META), delimiter="\t") if r and not r[0].startswith("#")]
    if meta and meta[0][0] == "species_taxid":
        meta = meta[1:]
    taxids = [r[0] for r in meta]
    genus = {r[0]: r[1] for r in meta}
    name = {r[0]: r[2] for r in meta}
    want = set(taxids)

    acc_by_tax = {}
    with open(PLAN) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["species_taxid"] in want:
                acc_by_tax.setdefault(row["species_taxid"], []).append(row["accession"])

    all_acc = [(t, a) for t in taxids for a in acc_by_tax.get(t, [])]
    print(f"{len(taxids)} taxa in {len(set(genus.values()))} genera, "
          f"{len(all_acc)} assemblies, {NPROC} procs", flush=True)
    with Pool(NPROC) as pool:
        res = pool.map(assembly_minimizers, [a for _, a in all_acc], chunksize=1)
    per = {}
    for (t, _), arr in zip(all_acc, res):
        per.setdefault(t, []).append(arr)

    global SPECIES
    SPECIES = [np.unique(np.concatenate(per[t])) if t in per else np.empty(0, np.uint64)
               for t in taxids]
    sizes = [s.size for s in SPECIES]
    missing = [taxids[i] for i, s in enumerate(sizes) if s == 0]
    print(f"minimizer sets built; {sum(sizes):,} distinct summed; "
          f"{len(missing)} taxa with no CDS on disk", flush=True)

    # ONLY within-genus pairs. Build the task list from the genus blocks, not from the full
    # upper triangle, or this is the 487x487 job it is explicitly not.
    idx_by_genus = {}
    for i, t in enumerate(taxids):
        idx_by_genus.setdefault(genus[t], []).append(i)
    tasks = []
    for g, idxs in idx_by_genus.items():
        for n, i in enumerate(idxs[:-1]):
            tasks.append((i, idxs[n + 1:]))
    npair = sum(len(js) for _, js in tasks)
    print(f"{npair:,} within-genus pairs", flush=True)

    shared = {}
    with Pool(NPROC) as pool:
        for i, out in pool.map(pair_block, tasks, chunksize=1):
            for j, v in out:
                shared[(i, j)] = v

    with open(OUT, "w") as fh:
        fh.write("genus\ttaxid_b\tname_b\ttaxid_a\tname_a\tdistinct_b\tdistinct_a\t"
                 "shared\tshare_b_given_a\n")
        for (i, j), v in sorted(shared.items()):
            # BOTH directions are written. Sharing is symmetric in COUNT but not as a
            # fraction: the same 50,000 shared minimizers are a fifth of a small genome and
            # a fiftieth of a large one, and the spillover a call suffers is governed by its
            # OWN denominator.
            for b, a in ((i, j), (j, i)):
                if not sizes[b]:
                    continue
                fh.write(f"{genus[taxids[b]]}\t{taxids[b]}\t{name[taxids[b]]}\t"
                         f"{taxids[a]}\t{name[taxids[a]]}\t{sizes[b]}\t{sizes[a]}\t"
                         f"{v}\t{v / sizes[b]:.6f}\n")
    top = sorted(((v / sizes[i], taxids[i], taxids[j]) for (i, j), v in shared.items()
                  if sizes[i]), reverse=True)[:5]
    print("highest share(B|A): " + ", ".join(f"{name[b]}|{name[a]} {s:.3f}"
                                             for s, b, a in top))
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
