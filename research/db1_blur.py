#!/usr/bin/env python3
"""The per-genus BLUR INDEX: how much of a genus's sequence cannot resolve below genus.

Kraken2 stores every minimizer at the lowest taxonomic node consistent with all genomes
carrying it. A minimizer unique to one species sits at that SPECIES node; one shared by two
or more species in a genus is pushed up to the GENUS node. So the genus node is, literally,
the sequence the congeners hold in common and cannot be told apart by.

For each genus, from kraken2-inspect's per-node minimizer counts (column 3, minimizers
assigned to that node alone, NOT the cumulative clade count in column 2):

    blur = (minimizers at genus and sub-genus nodes)
           -------------------------------------------------------------
           (those) + (minimizers at species and sub-species nodes below)

blur -> 0  species separate cleanly; keeping congeners is free; db1 buys recall only.
blur -> 1  the genus has collapsed; a species call is really a genus call; db1 costs
           resolution and species-level numbers in that genus should be rolled up.

Input is a kraken2-inspect TSV (6 cols: pct, clade_minimizers, taxon_minimizers, rank,
taxid, indented_name). The tree is read from the leading-space indentation of the name, two
spaces per level, exactly as kraken2-inspect emits it.

    python3 db1_blur.py <db_inspect.tsv>                 all genera, most-blurred first
    python3 db1_blur.py <db_inspect.tsv> --min-mins 1e6  only genera with real sequence
    python3 db1_blur.py <db_inspect.tsv> --genus Puccinia Alternaria Fusarium
"""

import argparse
import sys
from pathlib import Path


def parse(path):
    """Yield (depth, rank, taxid, name, taxon_minimizers) for each non-comment row."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            pct, clade_m, taxon_m, rank, taxid, name = line.rstrip("\n").split("\t")
            depth = (len(name) - len(name.lstrip(" "))) // 2
            yield depth, rank, taxid, name.strip(), int(taxon_m)


def blur_by_genus(rows):
    """Walk the DFS tree. For each genus node, accumulate its subtree's minimizers,
    split into blurring (genus/sub-genus rank) and resolving (species/sub-species rank).

    A genus's subtree is every following row at greater depth, up to the next row at the
    genus's own depth or shallower. A stack of open genera handles nesting cleanly.
    """
    genera = {}                       # taxid -> dict
    stack = []                        # open genera: (depth, taxid)

    for depth, rank, taxid, name, mins in rows:
        # close any genera whose subtree we have left
        while stack and depth <= stack[-1][0]:
            stack.pop()

        if rank == "G":
            g = {"name": name, "taxid": taxid, "depth": depth,
                 "genus_node": mins, "blur_mins": mins, "species_mins": 0,
                 "n_species": 0}
            genera[taxid] = g
            stack.append((depth, taxid))
            continue

        if not stack:
            continue                  # a node above any genus (family and up): ignore
        g = genera[stack[-1][1]]
        if rank.startswith("G"):      # G1, G2: sub-genus, still cannot resolve a species
            g["blur_mins"] += mins
        elif rank.startswith("S"):    # S, S1, S2: species or finer, resolving
            g["species_mins"] += mins
            if rank == "S":
                g["n_species"] += 1
        # ranks between genus and species other than G*/S* are rare; ignore safely

    for g in genera.values():
        denom = g["blur_mins"] + g["species_mins"]
        g["total_mins"] = denom
        g["blur"] = g["blur_mins"] / denom if denom else 0.0
    return genera


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inspect", help="kraken2-inspect TSV")
    ap.add_argument("--min-mins", type=float, default=0,
                    help="skip genera with fewer than this many total minimizers")
    ap.add_argument("--genus", nargs="+", help="only these genera, in this order")
    ap.add_argument("--top", type=int, default=30, help="rows to print (default 30)")
    args = ap.parse_args()

    genera = blur_by_genus(parse(args.inspect))
    rows = [g for g in genera.values() if g["total_mins"] >= args.min_mins]

    if args.genus:
        want = {n.lower() for n in args.genus}
        rows = [g for g in rows if g["name"].lower() in want]
        rows.sort(key=lambda g: [n.lower() for n in args.genus].index(g["name"].lower()))
    else:
        rows.sort(key=lambda g: -g["blur"])
        rows = rows[:args.top]

    print(f"\n{len(genera)} genera in {Path(args.inspect).name}; "
          f"showing {len(rows)}\n")
    print(f"  {'genus':<22}{'blur':>7}{'spp':>5}{'genus-node':>13}"
          f"{'species':>13}{'total mins':>14}")
    print(f"  {'':<22}{'':>7}{'':>5}{'(blur)':>13}{'(resolving)':>13}{'':>14}")
    for g in rows:
        print(f"  {g['name']:<22}{g['blur']:>7.2%}{g['n_species']:>5}"
              f"{g['blur_mins']:>13,}{g['species_mins']:>13,}{g['total_mins']:>14,}")

    print("\n  blur -> 0: congeners separate, db1 buys recall only.")
    print("  blur -> 1: genus collapsed, species calls unreliable, db1 costs resolution.")


if __name__ == "__main__":
    main()
