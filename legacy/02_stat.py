#!/usr/bin/env python3
"""
02_stat.py — fetch NCBI STAT taxonomy for all runs and apply retention gate.

Reads run accessions from output/01_sra/{mode}_runs.json (written by 01_sra.py),
fetches STAT for any not already cached, then applies the mode-specific
retention gate:

  MAL  keep runs where Viridiplantae reads >= MAL_MIN_HOST_PCT (1%).
       Confirms the pathogen was sequenced in host tissue, not pure culture.

  HAL  keep runs where >= 1 PHI-base plant pathogen is detected in STAT
       at >= HAL_MIN_PATHOGEN_PCT (1%). Confirms a pathogen is present in
       the host transcriptome.

Cache design (append-only, crash-safe, constant memory):
  stat_cache.jsonl      one line per fetched run: accession<TAB>json_data
                        shared between MAL and HAL — a run only fetched once
  stat_cache_index.txt  accessions only (~7 MB for 559k runs), read on startup
                        to resume without loading the full data file

Retention gate is applied by streaming stat_cache.jsonl line by line —
memory usage is O(1) regardless of dataset size.

Confirmed runs are written to output/02_stat/{mode}_confirmed.json.

Usage:
  python 02_stat.py --mode mal
  python 02_stat.py --mode hal
"""

import argparse
import json
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from _util import _Tee, load_json, save_json, http_get

# ── Settings ──────────────────────────────────────────────────────────────────

DB_PATH          = Path("output/00_build/data/phibase_db.json")
OUT_DIR          = Path("output/02_stat")
STAT_JSONL_PATH  = OUT_DIR / "data/stat_cache.jsonl"
STAT_INDEX_PATH  = OUT_DIR / "data/stat_cache_index.txt"

STAT_URL             = ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/"
                        "run_taxonomy?acc={run}&cluster_name=public")

MAL_MIN_HOST_PCT     = 1.0   # Viridiplantae % floor for in-planta confirmation
HAL_MIN_PATHOGEN_PCT = 1.0   # minimum % to count a PHI-base pathogen as detected
CACHE_MIN_PCT        = 0.1   # filter tax_table at write time; matches ABS_MIN_PCT in 03_crypt.py
LOG_EVERY            = 200   # print progress after this many new fetches
MAX_WORKERS          = 32

API_KEY = os.environ.get("NCBI_API_KEY", "")
RATE    = 15.0 if API_KEY else 5.0
HEADERS = {"User-Agent": "crypt/02_stat (leon.lenzo@curtin.edu.au)"}

HOST_NODE = "Viridiplantae"

_KINGDOM_KEYS = ("fungal_to_seed", "bacterial_to_seed", "oomycete_to_seed",
                 "nematode_to_seed", "virus_to_seed")


# ── HTTP ──────────────────────────────────────────────────────────────────────

def _get(url: str) -> bytes:
    return http_get(url, HEADERS)


# ── STAT fetch ────────────────────────────────────────────────────────────────

def _fetch_stat(run: str) -> list | None:
    url = STAT_URL.format(run=run)
    try:
        body = _get(url)
        return json.loads(body) if body else None
    except Exception as e:
        print(f"  STAT error {run}: {e}", flush=True)
        return None


_rate_lock = threading.Lock()
_rate_last: list[float] = [0.0]


def _throttled_fetch(run: str) -> tuple[str, list | None]:
    gap = 1.0 / RATE
    with _rate_lock:
        wait = _rate_last[0] + gap - time.time()
        if wait > 0:
            time.sleep(wait)
        _rate_last[0] = time.time()
    return run, _fetch_stat(run)


def _load_index(path: Path) -> set[str]:
    """Return set of already-fetched run accessions from the index file."""
    if not path.exists():
        return set()
    return set(path.read_text().split())


def _filter_stat(result: list | None) -> list | None:
    """Drop tax_table entries below CACHE_MIN_PCT before writing to cache."""
    if not result or not isinstance(result, list):
        return result
    block    = result[0]
    totals   = block.get("tax_totals", {})
    analyzed = totals.get("analysed", 0) or totals.get("analyzed", 0)
    if not analyzed:
        return result
    table    = block.get("tax_table", [])
    filtered = [e for e in table if e.get("total_count", 0) / analyzed * 100 >= CACHE_MIN_PCT]
    return [{"tax_totals": totals, "tax_table": filtered}]


def fetch_stat_append(to_fetch: list[str],
                      jsonl_path: Path, index_path: Path) -> None:
    """
    Fetch STAT for each accession in to_fetch; append immediately to jsonl_path
    and index_path. Each result is flushed to disk before moving to the next,
    so a crash loses at most whatever is in-flight across the worker threads.
    """
    if not to_fetch:
        print("  All runs already in cache.", flush=True)
        return

    print(f"  Fetching STAT for {len(to_fetch):,} runs …", flush=True)

    done = 0
    with open(jsonl_path, "a") as jf, open(index_path, "a") as ix:
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = {pool.submit(_throttled_fetch, acc): acc for acc in to_fetch}
            for future in as_completed(futures):
                run, result = future.result()
                jf.write(run + "\t" + json.dumps(_filter_stat(result)) + "\n")
                jf.flush()
                ix.write(run + "\n")
                ix.flush()
                done += 1
                if done % LOG_EVERY == 0 or done == len(to_fetch):
                    print(f"  [{done:,}/{len(to_fetch):,}] fetched  {time.strftime('%H:%M:%S')}", flush=True)


# ── STAT parsing ──────────────────────────────────────────────────────────────

def _parse_stat(data) -> dict | None:
    """
    Extract kingdom percentages from a STAT response.
    Returns None if data is absent, malformed, or has zero analysed reads.
    """
    if not data or not isinstance(data, list):
        return None
    block    = data[0]
    totals   = block.get("tax_totals", {})
    analyzed = totals.get("analysed", 0) or totals.get("analyzed", 0)
    if not analyzed:
        return None

    pct = {e["org"]: e.get("total_count", 0) / analyzed * 100
           for e in block.get("tax_table", []) if e.get("org")}

    return {
        "analyzed":      analyzed,
        "host_pct":      pct.get(HOST_NODE,    0.0),
        "fungi_pct":     pct.get("Fungi",      0.0),
        "virus_pct":     pct.get("Viruses",    0.0),
        "bacteria_pct":  pct.get("Bacteria",   0.0),
        "oomycete_pct":  pct.get("Oomycota",   0.0),
        "nematode_pct":  pct.get("Nematoda",   0.0),
        "_table":        block.get("tax_table", []),
    }


# ── Retention gates ───────────────────────────────────────────────────────────

def passes_mal_gate(stat: dict) -> bool:
    return stat["host_pct"] >= MAL_MIN_HOST_PCT


def passes_hal_gate(stat: dict, name_to_taxid: dict,
                    pathogen_taxids: set[int]) -> bool:
    analyzed = stat["analyzed"]
    for entry in stat["_table"]:
        name = entry.get("org", "")
        pct  = entry.get("total_count", 0) / analyzed * 100
        if pct < HAL_MIN_PATHOGEN_PCT:
            continue
        taxid = name_to_taxid.get(name.lower())
        if taxid and taxid in pathogen_taxids:
            return True
    return False


# ── PHI-base helpers (HAL gate only) ─────────────────────────────────────────

def _load_phibase_for_hal(db_path: Path) -> tuple[dict, set[int]]:
    if not db_path.exists():
        raise FileNotFoundError(
            f"PHI-base DB not found: {db_path}\nRun:  python3 00_build.py")
    with open(db_path) as f:
        raw = json.load(f)
    name_to_taxid   = raw["name_to_taxid"]
    pathogen_taxids = {int(k) for key in _KINGDOM_KEYS for k in raw[key]}
    print(f"PHI-base: {len(name_to_taxid):,} names, "
          f"{len(pathogen_taxids):,} pathogen taxids", flush=True)
    return name_to_taxid, pathogen_taxids


# ── Gate pass (streaming) ─────────────────────────────────────────────────────

def apply_gate(jsonl_path: Path, runs: dict, mode: str,
               name_to_taxid: dict | None,
               pathogen_taxids: set | None) -> tuple[dict, int, int]:
    """
    Stream stat_cache.jsonl line by line; apply retention gate.
    Returns (confirmed, n_no_stat, n_fail).
    Only processes lines whose accession appears in `runs` (this mode's runs).
    """
    confirmed  = {}
    n_no_stat  = n_fail = 0

    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                acc, data_str = line.split("\t", 1)
                data = json.loads(data_str)
            except (ValueError, json.JSONDecodeError):
                continue

            run_row = runs.get(acc)
            if run_row is None:
                continue  # fetched for other mode; skip

            stat = _parse_stat(data)
            if stat is None:
                n_no_stat += 1
                continue

            passes = (passes_mal_gate(stat) if mode == "mal"
                      else passes_hal_gate(stat, name_to_taxid, pathogen_taxids))
            if passes:
                row = dict(run_row)
                row.update({
                    "host_pct":     round(stat["host_pct"],     2),
                    "fungi_pct":    round(stat["fungi_pct"],    2),
                    "virus_pct":    round(stat["virus_pct"],    2),
                    "bacteria_pct": round(stat["bacteria_pct"], 2),
                    "oomycete_pct": round(stat["oomycete_pct"], 2),
                    "nematode_pct": round(stat["nematode_pct"], 2),
                    "analyzed":     stat["analyzed"],
                })
                confirmed[acc] = row
            else:
                n_fail += 1

    return confirmed, n_no_stat, n_fail


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True,
                    help="mal = microbe-as-library  |  hal = host-as-library")
    args = ap.parse_args()

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "logs").mkdir(parents=True, exist_ok=True)
    log = _Tee(OUT_DIR / "logs" / f"{args.mode}.log")
    sys.stdout = log

    try:
        runs_path      = Path(f"output/01_sra/data/{args.mode}_runs.json")
        confirmed_path = OUT_DIR / "data" / f"{args.mode}_confirmed.json"

        print(f"Mode: {args.mode.upper()}", flush=True)
        print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)

        if not runs_path.exists():
            raise FileNotFoundError(
                f"{runs_path} not found — run 01_sra.py --mode {args.mode} first")

        runs       = load_json(runs_path)
        accessions = list(runs.keys())
        print(f"Runs from 01_sra:    {len(accessions):,}", flush=True)

        # Resume: read index (accessions only) — fast even for 559k entries
        fetched    = _load_index(STAT_INDEX_PATH)
        to_fetch   = [a for a in accessions if a not in fetched]
        n_cached   = len(accessions) - len(to_fetch)
        print(f"STAT cache (shared): {len(fetched):,} entries total, "
              f"{n_cached:,} for this mode already cached", flush=True)

        fetch_stat_append(to_fetch, STAT_JSONL_PATH, STAT_INDEX_PATH)

        # Gate pass — PHI-base only needed here for HAL
        name_to_taxid = pathogen_taxids = None
        if args.mode == "hal":
            name_to_taxid, pathogen_taxids = _load_phibase_for_hal(DB_PATH)

        print(f"\nApplying {args.mode.upper()} gate (streaming stat_cache.jsonl) …",
              flush=True)
        confirmed, n_no_stat, n_fail = apply_gate(
            STAT_JSONL_PATH, runs, args.mode, name_to_taxid, pathogen_taxids)

        save_json(confirmed, confirmed_path)

        n_total     = len(accessions)
        n_stat      = n_total - n_no_stat
        n_confirmed = len(confirmed)
        gate_desc   = (f"≥{MAL_MIN_HOST_PCT}% Viridiplantae reads" if args.mode == "mal"
                       else f"PHI-base pathogen ≥{HAL_MIN_PATHOGEN_PCT}% detected")

        summary = (
            f"── 02_stat {args.mode.upper()} summary ─────────────────────────\n"
            f"Mode:              {args.mode.upper()}\n"
            f"\n"
            f"Total runs:        {n_total:>8,}\n"
            f"With STAT data:    {n_stat:>8,}  ({n_stat/max(n_total,1)*100:.1f}%)\n"
            f"No STAT data:      {n_no_stat:>8,}\n"
            f"Failed gate:       {n_fail:>8,}\n"
            f"Passed gate:       {n_confirmed:>8,}  "
            f"({n_confirmed/max(n_stat,1)*100:.1f}% of screened)\n"
            f"Gate:              {gate_desc}\n"
            f"\n"
            f"Confirmed: {confirmed_path}\n"
            f"Log:       {OUT_DIR / 'logs' / f'{args.mode}.log'}\n"
        )
        (OUT_DIR / "logs" / f"{args.mode}_summary.txt").write_text(summary)
        print(f"\n{summary}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
