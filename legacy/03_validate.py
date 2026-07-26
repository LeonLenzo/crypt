#!/usr/bin/env python3
"""
03_validate.py — empirical threshold calibration from STAT kingdom distributions.

Reads confirmed runs from output/02_stat/{mode}_confirmed.json, streams
stat_cache.jsonl to extract per-run kingdom percentages, species-level
detection percentages, and absolute k-mer counts (pct × analyzed / 100).

For each kingdom this script answers two questions:
  1. At what percentage does kingdom-level signal separate from background noise?
     → breakpoint table + ASCII histogram of kingdom pct across all confirmed runs
  2. Does read depth create false confidence at low percentages?
     → scatter summary of pct vs analyzed, showing absolute count at threshold

Output (output/03_validate/):
  {mode}_kingdom_dist.tsv    one row per run: run, analyzed, fungi_pct, virus_pct, ...
  {mode}_species_dist.tsv    one row per species detection: run, kingdom, name, pct, abs_count
  {mode}_summary.txt         breakpoint tables + threshold recommendations

Usage:
  python 03_validate.py --mode mal
  python 03_validate.py --mode hal
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

from _util import _Tee, load_json

# ── Paths ──────────────────────────────────────────────────────────────────────

DB_PATH         = Path("output/00_build/data/phibase_db.json")
STAT_JSONL_PATH = Path("output/02_stat/data/stat_cache.jsonl")
OUT_DIR         = Path("output/03_validate")

KINGDOMS = ["Fungi", "Viruses", "Bacteria", "Oomycota", "Nematoda"]

# Current hardcoded thresholds (what we're evaluating)
CURRENT_THRESHOLDS: dict[str, float] = {
    "Fungi":    1.0,
    "Viruses":  5.0,
    "Bacteria": 5.0,
    "Oomycota": 0.5,
    "Nematoda": 1.0,
}

ABS_MIN_PCT = 0.1
LEAF_FRAC   = 0.75

MYCOVIRUS_KEYWORDS = {
    "mycovirus", "mitovirus", "hypovirus", "chrysovirus", "partitivirus",
    "totivirus", "endornavirus", "victorivirus", "botybirnavirus",
    "mitoviridae", "narnaviridae", "hypoviridae", "mycovirales",
}

# Breakpoints for the summary table
BREAKPOINTS = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

# Minimum absolute k-mer count floor candidates to evaluate
ABS_COUNT_FLOORS = [100, 250, 500, 1000, 2000, 5000]


# ── STAT helpers (mirrors 04_crypt.py) ────────────────────────────────────────

def _analyzed(stat_data: list) -> int:
    t = stat_data[0].get("tax_totals", {})
    return t.get("analysed", 0) or t.get("analyzed", 0)


def _table(stat_data: list) -> list:
    return stat_data[0].get("tax_table", [])


def _node_count(table: list, node: str) -> int:
    for e in table:
        if e.get("org") == node:
            return e.get("total_count", 0)
    return 0


def _is_mycovirus(name: str) -> bool:
    n = name.lower()
    return any(kw in n for kw in MYCOVIRUS_KEYWORDS)


def specific_hits(table: list, node: str, analyzed: int,
                  min_pct: float = ABS_MIN_PCT) -> list[tuple[str, float]]:
    top = _node_count(table, node)
    if not top:
        return []
    min_count = max(1, int(analyzed * min_pct / 100))
    under = [
        (e["org"], e["total_count"]) for e in table
        if 0 < e.get("total_count", 0) <= top
        and e.get("org") != node
        and e.get("total_count", 0) >= min_count
    ]
    if not under:
        return [(node, top / analyzed * 100)]
    count_vals  = sorted({c for _, c in under}, reverse=True)
    leaf_counts = {
        c for c in count_vals
        if not any(c2 < c and c2 >= LEAF_FRAC * c for c2 in count_vals)
    }
    best: dict[int, tuple[str, float]] = {}
    for name, count in under:
        if count not in leaf_counts:
            continue
        pct = count / analyzed * 100
        nw  = len(name.split())
        if count not in best:
            best[count] = (name, pct)
        elif nw == 2 and len(best[count][0].split()) != 2:
            best[count] = (name, pct)
    return sorted(best.values(), key=lambda x: -x[1])


# ── Distribution helpers ───────────────────────────────────────────────────────

def _percentile(sorted_vals: list[float], p: float) -> float:
    if not sorted_vals:
        return 0.0
    idx = (len(sorted_vals) - 1) * p / 100
    lo, hi = int(idx), min(int(idx) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (idx - lo)


def _ascii_histogram(values: list[float], title: str, width: int = 50) -> str:
    """Log-scale histogram of values > 0, binned by powers of 2."""
    pos = [v for v in values if v > 0]
    if not pos:
        return f"{title}: no values > 0\n"

    # Bins: 0.1–0.2, 0.2–0.4, 0.4–0.8, 0.8–1.6, ... up to max
    lo = math.floor(math.log2(min(pos))) if min(pos) >= 1 else -4
    hi = math.ceil(math.log2(max(pos))) + 1 if max(pos) >= 1 else 1
    edges = [2 ** i for i in range(lo, hi + 1)]
    edges = [e for e in edges if e <= max(pos) * 2]

    counts = [0] * (len(edges) - 1)
    for v in pos:
        for i in range(len(edges) - 1):
            if edges[i] <= v < edges[i + 1]:
                counts[i] += 1
                break
        else:
            counts[-1] += 1

    max_count = max(counts) if counts else 1
    lines = [f"\n{title} (n={len(pos):,} runs with signal > 0, total={len(values):,})"]
    for i, c in enumerate(counts):
        label = f"{edges[i]:>7.2f}–{edges[i+1]:<7.2f}%"
        bar   = "█" * int(c / max_count * width)
        lines.append(f"  {label} | {bar} {c:,}")
    return "\n".join(lines) + "\n"


def _gap_threshold(sorted_vals: list[float]) -> float | None:
    """
    Find the largest relative gap in sorted positive values.
    Returns the value at the lower end of the gap, or None if distribution
    is too sparse to call a threshold.
    """
    pos = [v for v in sorted_vals if v > 0]
    if len(pos) < 10:
        return None
    # Look for gaps in log space between consecutive distinct values
    log_vals = sorted(set(round(math.log10(v), 3) for v in pos))
    if len(log_vals) < 3:
        return None
    max_gap = 0.0
    gap_lo  = None
    for i in range(len(log_vals) - 1):
        gap = log_vals[i + 1] - log_vals[i]
        if gap > max_gap:
            max_gap = gap
            gap_lo  = 10 ** log_vals[i]
    return gap_lo


def _breakpoint_table(kingdom: str, kpcts: list[float],
                      analyzed_vals: list[int]) -> str:
    """
    Summary table: at each threshold breakpoint, how many runs pass,
    and what is the median absolute k-mer count for passing runs.
    """
    n_total  = len(kpcts)
    cur_thresh = CURRENT_THRESHOLDS[kingdom]
    lines = [
        f"\n{kingdom} (n={n_total:,} confirmed runs, current threshold={cur_thresh}%)"
    ]
    lines.append(
        f"  {'Thresh%':>8}  {'Runs≥thresh':>12}  {'% of total':>10}  "
        f"{'Median abs count':>17}  {'P10 abs count':>14}"
    )
    lines.append("  " + "-" * 68)

    for bp in BREAKPOINTS:
        passing   = [(p, a) for p, a in zip(kpcts, analyzed_vals) if p >= bp]
        n_passing = len(passing)
        if n_passing == 0:
            lines.append(
                f"  {bp:>7.2f}%  {0:>12,}  {0:>9.1f}%  {'—':>17}  {'—':>14}"
            )
            continue
        abs_counts = sorted(p * a / 100 for p, a in passing)
        med_abs    = _percentile(abs_counts, 50)
        p10_abs    = _percentile(abs_counts, 10)
        marker     = " ← current" if bp == cur_thresh else ""
        lines.append(
            f"  {bp:>7.2f}%  {n_passing:>12,}  {n_passing/n_total*100:>9.1f}%  "
            f"{med_abs:>17,.0f}  {p10_abs:>14,.0f}{marker}"
        )
    return "\n".join(lines)


def _abs_count_table(kingdom: str, species_rows: list[dict]) -> str:
    """
    For detections at the current threshold, show what fraction survive
    various absolute count floors.
    """
    cur_thresh = CURRENT_THRESHOLDS[kingdom]
    at_thresh  = [r for r in species_rows
                  if r["kingdom"] == kingdom and r["pct"] >= cur_thresh]
    if not at_thresh:
        return f"\n{kingdom}: no species detections at current threshold"

    abs_counts = sorted(r["abs_count"] for r in at_thresh)
    lines = [
        f"\n{kingdom} species-level detections at ≥{cur_thresh}% "
        f"(n={len(at_thresh):,}): abs count floor impact"
    ]
    lines.append(f"  {'Min abs count':>14}  {'Survives':>10}  {'% kept':>8}")
    lines.append("  " + "-" * 38)
    for floor in ABS_COUNT_FLOORS:
        n_survive = sum(1 for c in abs_counts if c >= floor)
        lines.append(
            f"  {floor:>14,}  {n_survive:>10,}  {n_survive/len(at_thresh)*100:>7.1f}%"
        )
    return "\n".join(lines)


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True)
    args = ap.parse_args()
    mode = args.mode

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "logs").mkdir(parents=True, exist_ok=True)
    log = _Tee(OUT_DIR / "logs" / f"{mode}.log")
    sys.stdout = log

    try:
        confirmed_path = Path(f"output/02_stat/data/{mode}_confirmed.json")
        if not confirmed_path.exists():
            raise FileNotFoundError(
                f"{confirmed_path} not found — run 02_stat.py --mode {mode} first")

        confirmed   = json.loads(confirmed_path.read_text())
        needed      = set(confirmed.keys())
        n_confirmed = len(needed)
        print(f"Mode: {mode.upper()}")
        print(f"Confirmed runs: {n_confirmed:,}", flush=True)

        # ── Stream stat_cache.jsonl ──────────────────────────────────────────
        print(f"Streaming {STAT_JSONL_PATH} ...", flush=True)

        # Per-run kingdom pcts
        kingdom_rows: list[dict] = []
        # Per-species detections (at any pct > ABS_MIN_PCT)
        species_rows: list[dict] = []

        n_found  = 0
        n_nodata = 0

        with open(STAT_JSONL_PATH) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    acc, data_str = line.split("\t", 1)
                except ValueError:
                    continue
                if acc not in needed:
                    continue
                try:
                    stat_data = json.loads(data_str)
                except json.JSONDecodeError:
                    continue

                analyzed = _analyzed(stat_data)
                if not analyzed:
                    n_nodata += 1
                    continue

                table = _table(stat_data)
                n_found += 1

                row = {"run": acc, "analyzed": analyzed}
                for k in KINGDOMS:
                    kc     = _node_count(table, k)
                    kpct   = kc / analyzed * 100 if analyzed else 0.0
                    row[f"{k.lower()}_pct"] = round(kpct, 4)

                    # Species-level hits under this kingdom (at ABS_MIN_PCT floor)
                    for name, pct in specific_hits(table, k, analyzed):
                        if k == "Viruses" and _is_mycovirus(name):
                            continue
                        abs_count = pct * analyzed / 100
                        species_rows.append({
                            "run":       acc,
                            "kingdom":   k,
                            "name":      name,
                            "pct":       round(pct, 4),
                            "abs_count": round(abs_count, 1),
                            "analyzed":  analyzed,
                        })

                kingdom_rows.append(row)

                if n_found % 500 == 0:
                    print(f"  [{n_found:,}/{n_confirmed:,}] processed", flush=True)

        print(f"Processed: {n_found:,}  |  no STAT data: {n_nodata:,}", flush=True)

        # ── Write TSVs ───────────────────────────────────────────────────────
        kd_path  = OUT_DIR / "data" / f"{mode}_kingdom_dist.tsv"
        sp_path  = OUT_DIR / "data" / f"{mode}_species_dist.tsv"

        kd_fields = ["run", "analyzed"] + [f"{k.lower()}_pct" for k in KINGDOMS]
        with open(kd_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=kd_fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(kingdom_rows)

        sp_fields = ["run", "kingdom", "name", "pct", "abs_count", "analyzed"]
        with open(sp_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=sp_fields, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(species_rows)

        print(f"Written: {kd_path}  ({len(kingdom_rows):,} rows)")
        print(f"Written: {sp_path}  ({len(species_rows):,} rows)", flush=True)

        # ── Build summary ────────────────────────────────────────────────────
        summary_lines = [
            f"── 03_validate {mode.upper()} summary ──────────────────────",
            f"Confirmed runs processed: {n_found:,}",
            f"No STAT data:             {n_nodata:,}",
            "",
            "=" * 72,
            "KINGDOM-LEVEL BREAKPOINT TABLES",
            "=" * 72,
        ]

        analyzed_vals = [r["analyzed"] for r in kingdom_rows]

        for k in KINGDOMS:
            kpcts = [r[f"{k.lower()}_pct"] for r in kingdom_rows]
            summary_lines.append(_breakpoint_table(k, kpcts, analyzed_vals))

        summary_lines += [
            "",
            "=" * 72,
            "KINGDOM PERCENTAGE DISTRIBUTIONS (log2 bins)",
            "=" * 72,
        ]

        for k in KINGDOMS:
            kpcts = [r[f"{k.lower()}_pct"] for r in kingdom_rows]
            summary_lines.append(_ascii_histogram(kpcts, k))

            # Gap-based threshold suggestion
            sorted_pos = sorted(v for v in kpcts if v > 0)
            gap_val    = _gap_threshold(sorted_pos)
            if gap_val is not None:
                summary_lines.append(
                    f"  → Largest log-space gap below: {gap_val:.3f}%  "
                    f"(suggested empirical threshold)"
                )
            else:
                summary_lines.append(
                    f"  → Distribution too sparse / unimodal for gap detection"
                )

        summary_lines += [
            "",
            "=" * 72,
            "ABSOLUTE COUNT FLOOR ANALYSIS (species-level detections)",
            "=" * 72,
        ]

        for k in KINGDOMS:
            summary_lines.append(_abs_count_table(k, species_rows))

        summary_lines += [
            "",
            "=" * 72,
            "READ DEPTH SUMMARY (all confirmed runs)",
            "=" * 72,
        ]
        sorted_analyzed = sorted(analyzed_vals)
        summary_lines.append(
            f"  analyzed reads — "
            f"p10={_percentile(sorted_analyzed, 10):,.0f}  "
            f"median={_percentile(sorted_analyzed, 50):,.0f}  "
            f"p90={_percentile(sorted_analyzed, 90):,.0f}  "
            f"max={max(sorted_analyzed):,.0f}"
        )

        summary = "\n".join(summary_lines) + "\n"
        (OUT_DIR / "logs" / f"{mode}_summary.txt").write_text(summary)
        print(f"\n{summary}")

    finally:
        log.close()


if __name__ == "__main__":
    main()
