#!/usr/bin/env python3
"""kallisto_neighbours.py — measure which species share k-mer space, to decide index membership.

The original plan expanded each flagged candidate to its congeners. That is the wrong unit.
Phylogenetic closeness is an amino-acid-level statement, while Kraken2 and kallisto both match
exact 31-mers, and at congeneric divergence synonymous substitutions break a 31-mer every few
bases. So containment predicts cross-assignment and taxonomic rank does not.

containment(A in B) = |A n B| / |A| : the fraction of A's k-mer space a read could match in B.

Reuses the saturation module's cached FracMinHash sketches (k=31, scaled=1000), so this is a
seconds-to-minutes calculation that downloads and sketches nothing. Requires numpy, so run it
through SLURM with cray-python loaded, not on a login node.

    sbatch 04_kallisto/slurm/kallisto_neighbours.slurm

What it found on 2026-10-02 (see 04_kallisto/README.md): cross-talk above 1% is almost
entirely within-genus, Ascochyta/P. nodorum is 0.13% so the wheat Ascochyta detections are not
P. nodorum bleed, and 11 of 45 flagged species have no neighbour above 1% at all.
"""
import argparse
import collections
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import NEIGHBOURS, PATHO_LIN, REF_V3, SKETCHES, flagged, sp2

SCALED = 1000   # must match what kraken_saturation.py sketched with
RANKS = ("genus", "family", "order", "class", "phylum")


def species_accessions():
    """species -> [accession] for everything with a sketch on disk."""
    acc = collections.defaultdict(list)
    for r in csv.DictReader(open(REF_V3), delimiter="\t"):
        s = sp2(r["organism_name"])
        if len(s.split()) >= 2 and (SKETCHES / f"{r['accession']}.npy").exists():
            acc[s].append(r["accession"])
    return acc


def lineage():
    return {r["species"]: r for r in csv.DictReader(open(PATHO_LIN), delimiter="\t")}


def relationship(a, b, lin):
    la, lb = lin.get(a), lin.get(b)
    if not la or not lb:
        return "unknown"
    for rank in RANKS:
        if la.get(rank) and la[rank] == lb.get(rank):
            return rank
    return "beyond"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-containment", type=float, default=0.01,
                    help="floor for recording a pair (default 0.01; below 1%% is not worth "
                         "an index slot)")
    ap.add_argument("--out", type=Path, default=NEIGHBOURS)
    args = ap.parse_args()

    acc = species_accessions()
    lin = lineage()
    targets = [s for s in sorted({r["species"] for r in flagged()}) if s in acc]
    print(f"{len(acc)} species with sketches; {len(targets)} flagged species to profile",
          flush=True)

    cache = {}

    def union(s):
        """Species-level sketch: the union over its accessions' sketches."""
        if s not in cache:
            arrs = [np.load(SKETCHES / f"{a}.npy") for a in acc[s]]
            cache[s] = arrs[0] if len(arrs) == 1 else np.unique(np.concatenate(arrs))
        return cache[s]

    others, rows = sorted(acc), []
    for i, a in enumerate(targets, 1):
        ua = union(a)
        na = len(ua)
        for b in others:
            if b == a:
                continue
            ub = union(b)
            inter = len(np.intersect1d(ua, ub, assume_unique=True))
            if not inter or inter / na < args.min_containment:
                continue
            rows.append({
                "species": a, "neighbour": b, "relationship": relationship(a, b, lin),
                "containment_a_in_b": round(inter / na, 5),
                "containment_b_in_a": round(inter / len(ub), 5),
                "shared_kmers": int(inter * SCALED),
                "kmers_a": int(na * SCALED), "kmers_b": int(len(ub) * SCALED),
                "n_acc_a": len(acc[a]), "n_acc_b": len(acc[b]),
            })
        print(f"[{i}/{len(targets)}] {a}: "
              f"{sum(1 for r in rows if r['species'] == a)} neighbours "
              f">={args.min_containment:.0%}", flush=True)
        cache.pop(a, None)   # a flagged species is only needed once as the subject

    if not rows:
        sys.exit("no pairs above the containment floor; nothing written")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["species"], -r["containment_a_in_b"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    by_rel = collections.Counter(r["relationship"] for r in rows)
    print(f"\n{args.out}: {len(rows)} pairs")
    print(f"  by relationship: {dict(by_rel)}")
    no_neighbour = sorted(set(targets) - {r["species"] for r in rows})
    if no_neighbour:
        print(f"  {len(no_neighbour)} flagged species have NO neighbour above the floor, so "
              f"the competitive-survival test cannot run for them:")
        for s in no_neighbour:
            print(f"    {s}")


if __name__ == "__main__":
    main()
