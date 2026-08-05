#!/usr/bin/env python3
"""
figure/scatter/prep_scatter.py
Prepare scatter_data.tsv for scatter.R.

ALL points (background and foreground) use kingdom-level total_count from
stat_cache so that gate lines (x=1, y=1) consistently reflect the actual
gate criteria used in 02_filter_runs.py:
  MAL gate: Viridiplantae total_count / analyzed >= 1%  → x-axis
  HAL gate: (Fungi|Oomycota|Nematoda) total_count / analyzed >= 1%  → y-axis

Foreground (coloured): biosample_representative gate-pass runs from runs.tsv.
Background (grey):     all other stat_cache runs.

Run from crypt/: python figure/scatter/prep_scatter.py
"""

import csv
import json
from pathlib import Path

STAT_CACHE = Path("output/01_fetch_runs/data/stat_cache.jsonl")
RUNS_TSV   = Path("output/02_filter_runs/data/runs.tsv")
MAL_RUNS   = Path("output/01_fetch_runs/data/mal_runs.json")
HAL_RUNS   = Path("output/01_fetch_runs/data/hal_runs.json")
OUT_TSV    = Path("output/figures/scatter/scatter_data.tsv")

# ── Load gate-pass metadata from runs.tsv ─────────────────────────────────────

print("Loading runs.tsv ...", end=" ", flush=True)
gate_pass_all   = set()        # all gate-pass Run accessions (any biosample_rep value)
biosample_rep   = {}           # Run → {mode, coinf}  for biosample_representative rows only

with open(RUNS_TSV) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        run = r["Run"]
        gate_pass_all.add(run)
        if r.get("biosample_representative") == "True":
            biosample_rep[run] = {
                "mode":  r["mode"].upper(),
                "coinf": "True" if r["co_infection_flag"] != "single" else "False",
            }

print(f"{len(gate_pass_all):,} gate-pass runs, {len(biosample_rep):,} biosample-rep")

# ── Stream stat_cache — compute kingdom-level pcts for every run ───────────────

VIRIDIPLANTAE = "viridiplantae"
EUK_KINGDOMS  = {"fungi", "oomycota", "nematoda"}

def parse_stat(data: list) -> tuple[float, float]:
    """Kingdom-level host_pct and euk_pct — consistent with actual gate logic."""
    if not data or data[0] is None:
        return 0.0, 0.0
    totals   = data[0].get("tax_totals", {})
    analyzed = totals.get("analysed", 0) or totals.get("analyzed", 0)
    if not analyzed:
        return 0.0, 0.0
    host_cnt = euk_cnt = 0
    for entry in data[0].get("tax_table", []):
        org = entry.get("org", "").lower()
        tc  = entry.get("total_count", 0)
        if org == VIRIDIPLANTAE:
            host_cnt += tc
        elif org in EUK_KINGDOMS:
            euk_cnt  += tc
    return host_cnt / analyzed * 100, euk_cnt / analyzed * 100

# ── Load mode keys for background run assignment ───────────────────────────────

print("Loading mal_runs.json keys ...", end=" ", flush=True)
with open(MAL_RUNS) as f:
    mal_run_keys = set(json.load(f).keys())
print(f"{len(mal_run_keys):,}")

print("Loading hal_runs.json keys ...", end=" ", flush=True)
with open(HAL_RUNS) as f:
    hal_run_keys = set(json.load(f).keys())
print(f"{len(hal_run_keys):,}")

# ── Stream stat_cache ──────────────────────────────────────────────────────────

print("Streaming stat_cache.jsonl:", end=" ", flush=True)
fg_rows = []   # foreground: biosample_rep gate-pass
bg_rows = []   # background: everything else
seen_fg = set()
n_lines = 0

with open(STAT_CACHE) as f:
    for line in f:
        n_lines += 1
        if n_lines % 100_000 == 0:
            print(f"{n_lines // 1_000}k", end=" ", flush=True)
        try:
            acc, js = line.split("\t", 1)
        except ValueError:
            continue
        try:
            data = json.loads(js)
        except json.JSONDecodeError:
            continue
        host_pct, euk_pct = parse_stat(data)
        if host_pct == 0 and euk_pct == 0:
            continue
        if acc in biosample_rep:
            meta = biosample_rep[acc]
            fg_rows.append((meta["mode"], meta["coinf"], host_pct, euk_pct))
            seen_fg.add(acc)
        elif acc not in gate_pass_all:
            in_mal = acc in mal_run_keys
            in_hal = acc in hal_run_keys
            # only show in grey if it would genuinely fail the gate for that mode
            if in_mal and host_pct < 1.0:
                bg_rows.append(("MAL", host_pct, euk_pct))
            if in_hal and euk_pct < 1.0:
                bg_rows.append(("HAL", host_pct, euk_pct))

print(f"\n{n_lines:,} stat_cache lines")

missing = set(biosample_rep) - seen_fg
if missing:
    print(f"WARNING: {len(missing):,} biosample_rep runs not found in stat_cache")

print(f"Foreground: {len(fg_rows):,}  Background: {len(bg_rows):,}")

# ── Write output ──────────────────────────────────────────────────────────────

print(f"Writing {OUT_TSV} ...", end=" ", flush=True)
with open(OUT_TSV, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["layer", "mode", "coinf", "host_pct", "euk_pct"])
    for mode, coinf, host_pct, euk_pct in fg_rows:
        w.writerow(["gate", mode, coinf, host_pct, euk_pct])
    for mode, host_pct, euk_pct in bg_rows:
        w.writerow(["background", mode, "", host_pct, euk_pct])
print("done")
print(f"Written: {OUT_TSV}")
