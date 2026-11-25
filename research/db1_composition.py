#!/usr/bin/env python3
"""What is actually in db1_congener, and why did each genus qualify?

The scope rule reads USDA at GENUS, so db1 holds every assembly in a qualifying genus whose
species USDA never named on a plant. That lumps together two very different justifications,
and the question is which one carries the tier:

  A  COVER   the genus has USDA records that cannot be resolved to an NCBI species node, so
             the real organism behind those records may be any species in the genus. Adding
             congeners is covering a named-but-unresolvable observation.
  B  SPECULATIVE  every USDA record in the genus resolved cleanly and is already represented
             in db0. The db1 additions there are congeners never observed on a plant at all.

Reads  release/02_reference/data/curation/build_plan.tsv
       release/02_reference/data/curation/screen_species.tsv
Writes nothing. Prints the split.
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CUR = HERE.parent / "release/02_reference/data/curation"

UNRESOLVED_STATUS = {"no_ncbi_node"}


def genus_of(name):
    """First whitespace-delimited token of a USDA name. Crude, and that is the USDA rule."""
    return name.split()[0] if name.strip() else ""


def main():
    # 1. USDA evidence per genus
    resolved = Counter()        # genus -> USDA names that reached an NCBI species node
    unresolved = Counter()      # genus -> USDA names NCBI could not place
    bare_genus = Counter()      # genus -> USDA names that are only a genus
    with open(CUR / "screen_species.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = genus_of(r["usda_name"])
            if not g:
                continue
            if r["rank"] == "genus":
                bare_genus[g] += 1
            elif r["status"] in UNRESOLVED_STATUS or not r["ncbi_taxid"]:
                unresolved[g] += 1
            else:
                resolved[g] += 1

    # 2. db1 assemblies, grouped by genus
    db1 = defaultdict(list)
    db0_genera = set()
    with open(CUR / "build_plan.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["db_tier"] == "db1_congener":
                db1[r["genus"]].append(r["organism_name"])
            elif r["db_tier"] == "db0_baseline":
                db0_genera.add(r["genus"])

    n_asm = sum(len(v) for v in db1.values())
    print(f"db1_congener: {n_asm:,} assemblies across {len(db1):,} genera\n")

    # 3. the split
    cover_g, spec_g, nogenus_g = [], [], []
    for g in db1:
        if unresolved[g] or bare_genus[g]:
            cover_g.append(g)
        elif resolved[g]:
            spec_g.append(g)
        else:
            nogenus_g.append(g)       # genus not found in the USDA screen at all

    def block(label, genera, note):
        asm = sum(len(db1[g]) for g in genera)
        spp = len({s for g in genera for s in db1[g]})
        print(f"  {label:<34} {asm:>6,} assemblies  {spp:>5,} species  "
              f"{len(genera):>4} genera  {asm/n_asm:>6.1%}")
        print(f"  {'':<34} {note}")
        return genera, asm

    print("  justification for each db1 assembly, by its genus:\n")
    cover, n_cover = block("A  covers an unresolvable record", cover_g,
                           "USDA named something in this genus that NCBI cannot place")
    spec, n_spec = block("B  speculative congener", spec_g,
                         "every USDA record here resolved and is already in db0")
    _, n_none = block("?  genus absent from USDA screen", nogenus_g,
                      "admitted by an exception or a name the screen spells differently")

    # 4. the genera doing the most work, each way
    for label, genera in (("A  cover", cover), ("B  speculative", spec)):
        top = sorted(genera, key=lambda g: -len(db1[g]))[:12]
        print(f"\n  largest {label} genera:\n")
        for g in top:
            print(f"    {g:<24} {len(db1[g]):>5} assemblies   "
                  f"USDA resolved {resolved[g]:>4}, unresolvable {unresolved[g]:>4}")

    print(f"\n  {n_cover:,} of {n_asm:,} db1 assemblies ({n_cover/n_asm:.1%}) stand in for a "
          f"USDA record that exists but cannot be resolved.")
    print(f"  {n_spec:,} ({n_spec/n_asm:.1%}) are congeners with no unresolved USDA record "
          f"behind them.")


if __name__ == "__main__":
    main()
