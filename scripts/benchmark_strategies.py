#!/usr/bin/env python3
"""
benchmark_strategies.py - measure PMID discovery yield on unlabelled BioProjects.

Tests each strategy in 03b_fetch_literature.py independently on BioProjects that currently
have NO primary_pmid in the cache (i.e., the actual hard cases the pipeline
needs to solve). Reports what fraction of unlabelled entries each strategy
can find a PMID for, and how long it takes.

This is a discovery benchmark, not a recall benchmark: there is no ground truth
for unlabelled entries, so "yield" = any PMID returned.

Usage:
  python benchmark_strategies.py           # default 20 random unlabelled entries
  python benchmark_strategies.py --n 50    # sample 50

Run from crypt/:
  python benchmark_strategies.py
"""

import argparse
import importlib.util
import json
import random
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path


# ── Load 03b_fetch_literature.py as a module (main() is guarded by __name__=="__main__") ──

def _load_mod():
    root = Path(__file__).parent.parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    spec = importlib.util.spec_from_file_location(
        "fetch_lit03b", root / "03b_fetch_literature.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── BioProject XML fetch (shared across strategies that need bp_xml_pmids) ───

def _fetch_bp_xml_pmids(mod, uid: str) -> set[str]:
    bp_xml_pmids: set[str] = set()
    try:
        raw  = mod._get("efetch.fcgi", db="bioproject", id=uid)
        root = ET.fromstring(raw)
        for pub in root.findall(".//Publication"):
            pid = (pub.get("id") or "").strip()
            if pid.isdigit():
                bp_xml_pmids.add(pid)
    except Exception as e:
        print(f"    WARNING: efetch failed for uid={uid}: {e}", flush=True)
    return bp_xml_pmids


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--n", type=int, default=20,
        help="number of unlabelled BioProjects to sample (default: 20)",
    )
    ap.add_argument(
        "--seed", type=int, default=42,
        help="random seed for sampling (default: 42)",
    )
    args = ap.parse_args()

    # ── Load cache and runs.tsv ──────────────────────────────────────────────
    cache_path = Path("output/03b_fetch_literature/data/lit_cache.json")
    crypt_path = Path("output/02_filter_runs/data/runs.tsv")

    cache: dict = {}
    if cache_path.exists():
        with open(cache_path) as f:
            cache = json.load(f)

    # All BioProjects in runs.tsv
    all_bps: set[str] = set()
    if not crypt_path.exists():
        sys.exit(f"ERROR: {crypt_path} not found — run 02_filter_runs.py first")
    with open(crypt_path) as f:
        import csv
        for row in csv.DictReader(f, delimiter="\t"):
            bp = row.get("BioProject", "").strip()
            if bp:
                all_bps.add(bp)

    # Unlabelled = no primary_pmid in cache
    unlabelled = [
        bp for bp in sorted(all_bps)
        if not cache.get(bp, {}).get("primary_pmid")
    ]

    if not unlabelled:
        sys.exit("No unlabelled BioProjects found — all entries already have PMIDs!")

    random.seed(args.seed)
    sample = random.sample(unlabelled, min(args.n, len(unlabelled)))
    sample.sort()

    print(f"Benchmark: {len(sample)} unlabelled BioProjects "
          f"(sampled from {len(unlabelled)} with no PMID)", flush=True)
    print(flush=True)

    # ── Load 03b_fetch_literature.py ─────────────────────────────────────────
    try:
        mod = _load_mod()
    except Exception as e:
        sys.exit(f"ERROR: could not import 03b_fetch_literature.py: {e}")

    strategies = mod._STRATEGIES
    strat_results: dict[str, list[tuple[bool, float, list[str]]]] = {
        name: [] for name, _ in strategies
    }

    # ── Run benchmark ────────────────────────────────────────────────────────
    for idx, acc in enumerate(sample, 1):
        print(f"Testing [{idx:>3}/{len(sample)}] {acc} ...", flush=True)

        try:
            uid = mod._bp_uid(acc)
        except Exception as e:
            print(f"  WARNING: _bp_uid failed: {e} — skipping", flush=True)
            for name, _ in strategies:
                strat_results[name].append((False, 0.0, []))
            continue

        if uid is None:
            print(f"  WARNING: uid not found for {acc} — skipping", flush=True)
            for name, _ in strategies:
                strat_results[name].append((False, 0.0, []))
            continue

        bp_xml_pmids = _fetch_bp_xml_pmids(mod, uid)

        for name, fn in strategies:
            t0 = time.monotonic()
            try:
                pmids = fn(acc, uid, bp_xml_pmids)
            except Exception as e:
                print(f"  WARNING: strategy '{name}' raised: {e}", flush=True)
                pmids = []
            elapsed = time.monotonic() - t0
            hit = bool(pmids)
            strat_results[name].append((hit, elapsed, pmids[:3]))
            status = "HIT " if hit else "miss"
            print(f"  [{status}] {name:<22}  {elapsed:.2f}s"
                  + (f"  → {pmids[:3]}" if hit else ""),
                  flush=True)

        print(flush=True)

    # ── Results table ─────────────────────────────────────────────────────────
    n_total = len(sample)
    print("Discovery yield (sorted by yield):")
    sep = "─" * 57
    print(sep)
    print(f"{'Strategy':<24}  {'Found':>5}  {'Yield':>6}  {'Avg time':>9}")
    print(sep)

    stats = []
    for name, _ in strategies:
        results = strat_results[name]
        n = len(results)
        if n == 0:
            found, yield_pct, avg_time = 0, 0.0, 0.0
        else:
            found    = sum(1 for hit, _, _ in results if hit)
            yield_pct = 100.0 * found / n
            avg_time  = sum(t for _, t, _ in results) / n
        stats.append((name, found, yield_pct, avg_time))

    stats_sorted = sorted(stats, key=lambda x: (-x[2], x[3]))

    for name, found, yield_pct, avg_time in stats_sorted:
        print(f"  {name:<22}  {found:>5}/{n_total}  {yield_pct:>5.1f}%  {avg_time:>8.2f}s")

    print(sep)
    print(flush=True)

    bp_xml_entry = next((s for s in stats_sorted if s[0] == "BioProject XML"), None)
    others       = [s for s in stats_sorted if s[0] != "BioProject XML"]

    print("Suggested strategy order (BioProject XML always first — it's free):")
    order = []
    if bp_xml_entry:
        order.append(bp_xml_entry)
    order.extend(others)

    for rank, (name, found, yield_pct, _) in enumerate(order, 1):
        print(f"  {rank}. {name}  ({yield_pct:.1f}% yield on unlabelled entries)")

    print(flush=True)


if __name__ == "__main__":
    main()
