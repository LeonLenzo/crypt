#!/usr/bin/env python3
"""
02_stat.py — fetch NCBI STAT taxonomy for all runs and apply retention gate.

Reads run accessions from cache/{mode}_runs.json (written by 01_sra.py),
fetches STAT for any not already cached, then applies the mode-specific
retention gate:

  MAL  keep runs where Viridiplantae reads >= MAL_MIN_HOST_PCT (1%).
       Confirms the pathogen was sequenced in host tissue, not pure culture.

  HAL  keep runs where >= 1 PHI-base plant pathogen is detected in STAT
       at >= HAL_MIN_PATHOGEN_PCT (1%). Confirms a pathogen is present in
       the host transcriptome.

STAT responses are cached in cache/stat_cache.json (shared across both
modes — a run only needs to be fetched once).

Confirmed runs are written to cache/{mode}_confirmed.json.

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

from _util import _Tee, http_get, load_json, save_json

# ── Settings ──────────────────────────────────────────────────────────────────

DB_PATH         = Path("output/00_build/phibase_db.json")
OUT_DIR         = Path("output/02_stat")
STAT_CACHE_PATH = OUT_DIR / "stat_cache.json"

STAT_URL             = ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/"
                        "run_taxonomy?acc={run}&cluster_name=public")

MAL_MIN_HOST_PCT     = 1.0   # Viridiplantae % floor for in-planta confirmation
HAL_MIN_PATHOGEN_PCT = 1.0   # minimum % to count a PHI-base pathogen as detected
SAVE_EVERY           = 200   # save stat_cache after this many new fetches
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


def fetch_stat_parallel(accessions: list[str], stat_cache: dict,
                        stat_path: Path) -> None:
    """Fetch STAT for accessions not already in stat_cache, saving periodically."""
    to_fetch = [a for a in accessions if a not in stat_cache]
    if not to_fetch:
        print(f"  All {len(accessions):,} runs already in STAT cache.", flush=True)
        return

    print(f"  Fetching STAT for {len(to_fetch):,} runs "
          f"({len(accessions) - len(to_fetch):,} cached) …", flush=True)

    done = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(_throttled_fetch, acc): acc for acc in to_fetch}
        for future in as_completed(futures):
            run, result = future.result()
            stat_cache[run] = result
            done += 1
            if done % SAVE_EVERY == 0 or done == len(to_fetch):
                save_json(stat_cache, stat_path)
                print(f"  [{done:,}/{len(to_fetch):,}] fetched", flush=True)


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
    """MAL: plant reads must be present to confirm the run is in planta."""
    return stat["host_pct"] >= MAL_MIN_HOST_PCT


def passes_hal_gate(stat: dict, name_to_taxid: dict,
                    pathogen_taxids: set[int]) -> bool:
    """HAL: at least one PHI-base plant pathogen must be detected in STAT."""
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


# ── PHI-base helpers (HAL only) ───────────────────────────────────────────────

def _load_phibase_for_hal(db_path: Path) -> tuple[dict, set[int]]:
    """Return (name_to_taxid, pathogen_taxids) from phibase_db.json."""
    if not db_path.exists():
        raise FileNotFoundError(
            f"PHI-base DB not found: {db_path}\nRun:  python3 00_build.py")
    with open(db_path) as f:
        raw = json.load(f)
    name_to_taxid   = raw["name_to_taxid"]
    pathogen_taxids = {
        int(k)
        for key in _KINGDOM_KEYS
        for k in raw[key]
    }
    print(f"PHI-base: {len(name_to_taxid):,} names, "
          f"{len(pathogen_taxids):,} pathogen taxids", flush=True)
    return name_to_taxid, pathogen_taxids


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True,
                    help="mal = microbe-as-library  |  hal = host-as-library")
    args = ap.parse_args()

    log = _Tee(OUT_DIR / f"{args.mode}.log")
    sys.stdout = log

    try:
        runs_path      = Path(f"output/01_sra/{args.mode}_runs.json")
        confirmed_path = OUT_DIR / f"{args.mode}_confirmed.json"

        print(f"Mode: {args.mode.upper()}", flush=True)

        if not runs_path.exists():
            raise FileNotFoundError(
                f"{runs_path} not found — run 01_sra.py --mode {args.mode} first")

        runs       = load_json(runs_path)
        accessions = list(runs.keys())
        print(f"Runs from 01_sra:    {len(accessions):,}", flush=True)

        name_to_taxid = pathogen_taxids = None
        if args.mode == "hal":
            name_to_taxid, pathogen_taxids = _load_phibase_for_hal(DB_PATH)

        stat_cache = load_json(STAT_CACHE_PATH)
        print(f"STAT cache (shared): {len(stat_cache):,} entries", flush=True)
        fetch_stat_parallel(accessions, stat_cache, STAT_CACHE_PATH)

        print(f"\nApplying {args.mode.upper()} retention gate …", flush=True)
        confirmed = {}
        n_no_stat = n_fail = 0

        for acc, run_row in runs.items():
            stat = _parse_stat(stat_cache.get(acc))
            if stat is None:
                n_no_stat += 1
                continue
            passes = (passes_mal_gate(stat) if args.mode == "mal"
                      else passes_hal_gate(stat, name_to_taxid, pathogen_taxids))
            if passes:
                row = dict(run_row)
                row["host_pct"]      = round(stat["host_pct"],      2)
                row["fungi_pct"]     = round(stat["fungi_pct"],     2)
                row["virus_pct"]     = round(stat["virus_pct"],     2)
                row["bacteria_pct"]  = round(stat["bacteria_pct"],  2)
                row["oomycete_pct"]  = round(stat["oomycete_pct"],  2)
                row["nematode_pct"]  = round(stat["nematode_pct"],  2)
                row["analyzed"]      = stat["analyzed"]
                confirmed[acc] = row
            else:
                n_fail += 1

        save_json(confirmed, confirmed_path)

        n_total     = len(accessions)
        n_stat      = n_total - n_no_stat
        n_confirmed = len(confirmed)

        gate_desc = (f"≥{MAL_MIN_HOST_PCT}% Viridiplantae reads" if args.mode == "mal"
                     else f"PHI-base pathogen ≥{HAL_MIN_PATHOGEN_PCT}% detected")
        summary = (
            f"── 02_stat {args.mode.upper()} summary ─────────────────────────\n"
            f"Mode:              {args.mode.upper()}\n"
            f"\n"
            f"Total runs:        {n_total:>8,}\n"
            f"With STAT data:    {n_stat:>8,}  ({n_stat/max(n_total,1)*100:.1f}%)\n"
            f"No STAT data:      {n_no_stat:>8,}\n"
            f"Passed gate:       {n_confirmed:>8,}  "
            f"({n_confirmed/max(n_stat,1)*100:.1f}% of screened)\n"
            f"Gate:              {gate_desc}\n"
            f"\n"
            f"Confirmed: {confirmed_path}\n"
            f"Log:       {OUT_DIR / f'{args.mode}.log'}\n"
        )
        (OUT_DIR / f"{args.mode}_summary.txt").write_text(summary)
        print(f"\n{summary}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
