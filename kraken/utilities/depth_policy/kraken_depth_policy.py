#!/usr/bin/env python3
"""
kraken_db_depth_policy.py — decide the per-species assembly depth policy for a
Kraken2 build: raise the shallow species (a floor), or trim the deep ones (a cap)?

Both change the same quantity, the DIFFERENCE in sampling depth between congeners,
which is what actually misassigns reads: a k-mer present in two species' pangenomes
but absent from the shallow one's included assemblies is stored as the deep species,
so the shallow species' reads score for the deep one.

This script does no sketching. It reads the cached outputs of
kraken_db_saturation.py and answers three questions the build design turns on:

  1. What does a cap BUY?   For every congener pair, misassignment with the deep
     species truncated to --cap, against its full depth. Sweeping deep_depth is
     already what depth_bias.tsv records, so the cap needs no new measurement.

  2. What does a cap COST?  The fraction of a deep species' observed pangenome
     retained at --cap, from saturation_curve.tsv. Pangenomes are open (Heaps
     gamma ~0.16), so this is a real loss, not a redundancy trim.

  3. Is the cap even NECESSARY?  A pair only needs the cap if the floor cannot
     fix it, i.e. if the shallow partner has fewer than --floor assemblies in
     existence. holobase supplies that count.

Output: data/depth_policy.tsv, one row per implicated congener pair.

Usage (reads from the repo layout by default; --bias/--curve accept local copies):
    python kraken/utilities/depth_policy.py
    python kraken/utilities/depth_policy.py --cap 5 --floor 4
"""

import argparse
import csv
import sqlite3
import statistics as st
from collections import defaultdict
from pathlib import Path

OUT_DIR   = Path("kraken/utilities/saturation")
DATA_DIR  = OUT_DIR / "data"
BIAS_TSV  = DATA_DIR / "depth_bias.tsv"
CURVE_TSV = DATA_DIR / "saturation_curve.tsv"
POLICY_TSV = DATA_DIR / "depth_policy.tsv"
HOLOBASE  = Path.home() / "data_analysis" / "_refdata" / "holobase.db"

# A pair below this is noise: the median pair sits at 0.08%, so sub-1% differences
# do not change any call. Only pairs above it can justify a depth policy at all.
IMPLICATED = 0.01


def load_pairs(bias_tsv):
    """depth_bias.tsv -> {(deep, shallow): {deep_depth: row}}."""
    pairs = defaultdict(dict)
    with open(bias_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pairs[(r["deep_species"], r["shallow_species"])][int(r["deep_depth"])] = r
    return pairs


def load_retained(curve_tsv, cap):
    """saturation_curve.tsv -> {species: fraction of observed pangenome kept at cap}.

    quality_kmers is cumulative over assemblies in quality order, so the value at
    n=cap over the value at full depth is what a cap would retain.
    """
    cum = defaultdict(dict)
    with open(curve_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cum[r["organism_name"]][int(r["n_assemblies"])] = float(r["quality_kmers"])
    out = {}
    for name, d in cum.items():
        full = max(d)
        if full > cap and cap in d and d[full] > 0:
            out[name] = (d[cap] / d[full], full)
    return out


def load_available(db_path):
    """holobase -> {species name: assemblies in existence}, max over homonym rows."""
    if not db_path.exists():
        return {}
    con = sqlite3.connect(db_path)
    q = """SELECT t.scientific_name, COUNT(a.accession)
             FROM taxon t LEFT JOIN assembly a ON a.taxon_id = t.taxon_id
            WHERE t.rank = 'species' GROUP BY t.taxon_id"""
    out = {}
    for name, n in con.execute(q):
        if n > (out.get(name) or 0):
            out[name] = n
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bias", default=str(BIAS_TSV))
    ap.add_argument("--curve", default=str(CURVE_TSV))
    ap.add_argument("--holobase", default=str(HOLOBASE))
    ap.add_argument("--out", default=str(POLICY_TSV))
    ap.add_argument("--cap", type=int, default=5,
                    help="depth a cap would truncate deep species to (default 5)")
    ap.add_argument("--floor", type=int, default=4,
                    help="depth a floor would raise shallow species to (default 4)")
    args = ap.parse_args()

    pairs = load_pairs(args.bias)
    retained = load_retained(args.curve, args.cap)
    available = load_available(Path(args.holobase))

    rows, cap_helps, cap_hurts = [], 0, 0
    for (deep, shallow), depths in pairs.items():
        full = max(depths)
        if full <= args.cap:
            continue                      # the cap cannot touch this pair
        at_cap = max(d for d in depths if d <= args.cap)

        # b1 = shallow species left at one assembly, the state a broad sweep puts
        # 70% of new species in. b4 = shallow species raised to the floor.
        f_cap  = float(depths[at_cap]["frac_misassigned_at_b1"])
        f_full = float(depths[full]["frac_misassigned_at_b1"])
        if f_full <= IMPLICATED:
            continue

        if f_cap < f_full:
            cap_helps += 1
        elif f_cap > f_full:
            cap_hurts += 1

        keep, n_asm = retained.get(deep, (None, full))
        n_avail = available.get(shallow)
        rows.append({
            "deep_species": deep,
            "shallow_species": shallow,
            "deep_assemblies": n_asm,
            "misassigned_deep_full": round(f_full, 6),
            "misassigned_deep_capped": round(f_cap, 6),
            "cap_benefit_pp": round((f_full - f_cap) * 100, 3),
            "misassigned_shallow_at_floor": round(
                float(depths[full]["frac_misassigned_at_b4"]), 6),
            "floor_benefit_pp": round(
                (f_full - float(depths[full]["frac_misassigned_at_b4"])) * 100, 3),
            "deep_pangenome_kept_at_cap": round(keep, 4) if keep is not None else "",
            "deep_pangenome_lost_at_cap_pp": round((1 - keep) * 100, 1)
                                             if keep is not None else "",
            "shallow_assemblies_available": n_avail if n_avail is not None else "",
            "floor_reachable": "" if n_avail is None else str(n_avail >= args.floor),
        })

    rows.sort(key=lambda r: -r["misassigned_deep_full"])
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    cap_b   = [r["cap_benefit_pp"] for r in rows]
    floor_b = [r["floor_benefit_pp"] for r in rows]
    lost    = [r["deep_pangenome_lost_at_cap_pp"] for r in rows
               if r["deep_pangenome_lost_at_cap_pp"] != ""]
    reach   = [r for r in rows if r["floor_reachable"] == "True"]

    print(f"implicated pairs (>{IMPLICATED*100:g}% misassigned at full depth): {len(rows)}")
    print(f"  deep species involved: "
          f"{len({r['deep_species'] for r in rows})}")
    print(f"\ncap to {args.cap}:")
    print(f"  reduces misassignment in {cap_helps} pairs, increases it in {cap_hurts}")
    print(f"  median benefit  {st.median(cap_b):.2f} pp   max {max(cap_b):.2f} pp")
    print(f"  median cost     {st.median(lost):.1f}% of the deep species' pangenome"
          f"   max {max(lost):.1f}%")
    print(f"\nfloor at {args.floor}:")
    print(f"  median benefit  {st.median(floor_b):.2f} pp   max {max(floor_b):.2f} pp")
    print(f"  cost            none, it only adds assemblies")
    print(f"  reachable for   {len(reach)}/{len(rows)} pairs")
    print(f"\nwrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
