#!/usr/bin/env python3
"""
01_fetch_runs.py — fetch SRA run IDs and NCBI STAT taxonomy for co-infection mining.

Merges the SRA ID fetch (01_sra.py) and STAT data fetch (02_stat.py) into a
single data-gathering step.  No retention gate is applied here — that belongs
in 02_filter_runs.py.

Two modes:
  mal  microbe-as-library: PHI-base plant pathogen as library organism.
  hal  host-as-library:    PHI-base plant host as library organism.

Output (output/01_fetch_runs/):
  data/{mode}_runs.json      RunInfo rows keyed by Run accession
  data/{mode}_uids.json      fetched UID set (resumability)
  data/stat_cache.jsonl      STAT responses, append-only, shared MAL+HAL
  data/stat_cache_index.txt  accessions only — read on startup to resume
  logs/{mode}.log

Do NOT run MAL and HAL simultaneously — they share stat_cache.jsonl.
Chain HAL after MAL:
  until grep -q 'MAL summary' output/01_fetch_runs/logs/latest/mal_summary.txt; do sleep 30; done
  python 01_fetch_runs.py --mode hal

Usage:
  python 01_fetch_runs.py --mode mal
  python 01_fetch_runs.py --mode hal
  python 01_fetch_runs.py --mode mal --count   (hit counts only, no fetch)
"""

import argparse
import csv
import io
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from _util import _Tee, http_get, link_latest, load_json, make_log_dir, save_json

# ── Settings ──────────────────────────────────────────────────────────────────

DB_PATH  = Path("output/00_build/data/phibase_db.json")
OUT_DIR  = Path("output/01_fetch_runs")

LIBRARY_STRAT    = "RNA-Seq"
MAL_BATCH        = 50
RUNINFO_BATCH    = 500
SRA_MAX_WORKERS  = 8
SAVE_EVERY       = 20

STAT_URL         = ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/"
                    "run_taxonomy?acc={run}&cluster_name=public")
STAT_MAX_WORKERS = 32
CACHE_MIN_PCT    = 0.1   # filter tax_table at write time; matches ABS_MIN_PCT in 02_filter_runs.py
LOG_EVERY        = 200

API_KEY  = os.environ.get("NCBI_API_KEY", "")
RATE     = 9.0  if API_KEY else 2.5
STAT_RATE = 15.0 if API_KEY else 5.0
HEADERS  = {"User-Agent": "crypt/01_fetch_runs (leon.lenzo@curtin.edu.au)"}
ENTREZ   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_MAL_KINGDOM_KEYS = [
    ("Fungi",    "fungal_to_seed"),
    ("Bacteria", "bacterial_to_seed"),
    ("Oomycota", "oomycete_to_seed"),
    ("Nematoda", "nematode_to_seed"),
    ("Viruses",  "virus_to_seed"),
]


# ── HTTP helpers ──────────────────────────────────────────────────────────────

def _get(url: str) -> bytes:
    return http_get(url, HEADERS)


def _post(endpoint: str, retries: int = 5, **params) -> bytes:
    if API_KEY:
        params["api_key"] = API_KEY
    url  = f"{ENTREZ}/{endpoint}"
    body = urllib.parse.urlencode(params).encode()
    hdrs = {**HEADERS, "Content-Type": "application/x-www-form-urlencoded"}
    delay = 1.0
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, data=body, headers=hdrs)
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(delay); delay *= 2
            else:
                raise
        except Exception:
            time.sleep(delay); delay *= 2
    raise RuntimeError(f"POST failed after {retries} retries: {endpoint}")


def _url(endpoint: str, **params) -> str:
    if API_KEY:
        params["api_key"] = API_KEY
    return f"{ENTREZ}/{endpoint}?{urllib.parse.urlencode(params)}"


def _esearch_count(db: str, term: str) -> int:
    url  = _url("esearch.fcgi", db=db, term=term, retmax=0, retmode="json")
    data = json.loads(_get(url))["esearchresult"]
    time.sleep(1.0 / RATE)
    return int(data["count"])


def _esearch(db: str, term: str, retmax: int = 10_000) -> list[str]:
    url   = _url("esearch.fcgi", db=db, term=term,
                 usehistory="y", retmax=0, retmode="json")
    data  = json.loads(_get(url))["esearchresult"]
    total = int(data["count"])
    wk, qk = data["webenv"], data["querykey"]
    uids: list[str] = []
    for start in range(0, total, retmax):
        url = _url("efetch.fcgi", db=db, query_key=qk, WebEnv=wk,
                   retstart=start, retmax=retmax, rettype="uilist", retmode="text")
        uids.extend(_get(url).decode().strip().splitlines())
        time.sleep(1.0 / RATE)
    return uids


# ── DB helpers ────────────────────────────────────────────────────────────────

def _load_db() -> dict:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"PHI-base DB not found: {DB_PATH}\nRun:  python3 00_build.py")
    with open(DB_PATH) as f:
        return json.load(f)


def _load_seeds(mode: str) -> list[int]:
    raw = _load_db()
    if mode == "mal":
        seeds = sorted(int(k) for k in raw["pathogen_to_hosts"])
        label = "seed pathogen species"
    else:
        seeds = sorted(set(raw["host_to_seed"].values()))
        label = "seed host species"
    print(f"PHI-base: {len(seeds)} {label}", flush=True)
    return seeds


def _load_seeds_by_kingdom() -> dict[str, list[int]]:
    raw = _load_db()
    return {
        kingdom: sorted(set(raw[db_key].values()))
        for kingdom, db_key in _MAL_KINGDOM_KEYS
    }


# ── SRA query ─────────────────────────────────────────────────────────────────

def _build_queries(seeds: list[int]) -> list[str]:
    queries: list[str] = []
    for i in range(0, len(seeds), MAL_BATCH):
        batch    = seeds[i:i + MAL_BATCH]
        or_terms = " OR ".join(f"txid{t}[Organism:exp]" for t in batch)
        queries.append(f"({or_terms}) AND {LIBRARY_STRAT}[Strategy]")
    return queries


def count_hits(seeds: list[int]) -> int:
    queries = _build_queries(seeds)
    total   = 0
    for i, query in enumerate(queries, 1):
        n      = _esearch_count("sra", query)
        total += n
        print(f"  [{i}/{len(queries)}] {n:,} hits  (running total {total:,})", flush=True)
    return total


def fetch_uids(seeds: list[int]) -> list[str]:
    queries  = _build_queries(seeds)
    all_uids: list[str] = []
    for i, query in enumerate(queries, 1):
        uids = _esearch("sra", query)
        all_uids.extend(uids)
        print(f"  [{i}/{len(queries)}] +{len(uids):,} UIDs  "
              f"(total {len(all_uids):,})", flush=True)
    return all_uids


# ── RunInfo fetch ─────────────────────────────────────────────────────────────

_ri_lock = threading.Lock()
_ri_last: list[float] = [0.0]


def _fetch_runinfo_batch(batch: list[str]) -> tuple[list[str], list[dict]]:
    gap = 1.0 / RATE
    with _ri_lock:
        wait = _ri_last[0] + gap - time.time()
        if wait > 0:
            time.sleep(wait)
        _ri_last[0] = time.time()
    text = _post("efetch.fcgi", db="sra", id=",".join(batch),
                 rettype="runinfo", retmode="text").decode(errors="replace")
    rows = [r for r in csv.DictReader(io.StringIO(text)) if r.get("Run")]
    return batch, rows


def fetch_runinfo(uids: list[str], run_cache: dict, uid_cache: dict,
                  run_path: Path, uid_path: Path) -> int:
    to_fetch = [u for u in uids if u not in uid_cache]
    if not to_fetch:
        print(f"  All {len(uids):,} UIDs already fetched.", flush=True)
        return 0
    batches = [to_fetch[i:i + RUNINFO_BATCH]
               for i in range(0, len(to_fetch), RUNINFO_BATCH)]
    print(f"  Fetching RunInfo for {len(to_fetch):,} UIDs in "
          f"{len(batches):,} batches ({len(uids) - len(to_fetch):,} cached) …",
          flush=True)
    added = done = 0
    with ThreadPoolExecutor(max_workers=SRA_MAX_WORKERS) as pool:
        futures = {pool.submit(_fetch_runinfo_batch, b): b for b in batches}
        for future in as_completed(futures):
            batch_uids, rows = future.result()
            for row in rows:
                acc = row["Run"]
                if acc not in run_cache:
                    run_cache[acc] = row
                    added += 1
            for u in batch_uids:
                uid_cache[u] = True
            done += 1
            if done % SAVE_EVERY == 0 or done == len(batches):
                save_json(run_cache, run_path)
                save_json(uid_cache, uid_path)
                print(f"  [{done:,}/{len(batches):,} batches]  "
                      f"{added:,} new runs", flush=True)
    return added


# ── STAT fetch ────────────────────────────────────────────────────────────────

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
    filtered = [e for e in table
                if e.get("total_count", 0) / analyzed * 100 >= CACHE_MIN_PCT]
    return [{"tax_totals": totals, "tax_table": filtered}]


def _fetch_stat(run: str) -> list | None:
    url = STAT_URL.format(run=run)
    try:
        body = http_get(url, HEADERS)
        return json.loads(body) if body else None
    except Exception as e:
        print(f"  STAT error {run}: {e}", flush=True)
        return None


_stat_lock = threading.Lock()
_stat_last: list[float] = [0.0]


def _throttled_stat(run: str) -> tuple[str, list | None]:
    gap = 1.0 / STAT_RATE
    with _stat_lock:
        wait = _stat_last[0] + gap - time.time()
        if wait > 0:
            time.sleep(wait)
        _stat_last[0] = time.time()
    return run, _fetch_stat(run)


def _load_index(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return set(path.read_text().split())


def fetch_stat_append(to_fetch: list[str],
                      jsonl_path: Path, index_path: Path) -> None:
    if not to_fetch:
        print("  All runs already in STAT cache.", flush=True)
        return
    print(f"  Fetching STAT for {len(to_fetch):,} runs …", flush=True)
    done = 0
    with open(jsonl_path, "a") as jf, open(index_path, "a") as ix:
        with ThreadPoolExecutor(max_workers=STAT_MAX_WORKERS) as pool:
            futures = {pool.submit(_throttled_stat, acc): acc for acc in to_fetch}
            for future in as_completed(futures):
                run, result = future.result()
                jf.write(run + "\t" + json.dumps(_filter_stat(result)) + "\n")
                jf.flush()
                ix.write(run + "\n")
                ix.flush()
                done += 1
                if done % LOG_EVERY == 0 or done == len(to_fetch):
                    print(f"  [{done:,}/{len(to_fetch):,}] fetched  "
                          f"{time.strftime('%H:%M:%S')}", flush=True)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True,
                    help="mal = microbe-as-library  |  hal = host-as-library")
    ap.add_argument("--count", action="store_true",
                    help="print hit counts only; do not fetch")
    args = ap.parse_args()

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)

    stat_jsonl = OUT_DIR / "data" / "stat_cache.jsonl"
    stat_index = OUT_DIR / "data" / "stat_cache_index.txt"
    run_path   = OUT_DIR / "data" / f"{args.mode}_runs.json"
    uid_path   = OUT_DIR / "data" / f"{args.mode}_uids.json"

    log = _Tee(log_dir / f"{args.mode}.log")
    link_latest(logs_base, log_dir / f"{args.mode}.log")
    sys.stdout = log

    try:
        print(f"Mode: {args.mode.upper()}", flush=True)
        print(f"NCBI_API_KEY set: {'yes' if API_KEY else 'no'}", flush=True)
        seeds = _load_seeds(args.mode)

        if args.count:
            if args.mode == "mal":
                by_kingdom = _load_seeds_by_kingdom()
                grand = 0
                for kingdom, k_seeds in by_kingdom.items():
                    print(f"\n── {kingdom} ({len(k_seeds)} seeds) ──", flush=True)
                    n = count_hits(k_seeds)
                    print(f"  {kingdom} subtotal: {n:,}", flush=True)
                    grand += n
                print(f"\nMAL grand total (batched, may overcount): {grand:,}",
                      flush=True)
            else:
                n = count_hits(seeds)
                print(f"\nHAL total (batched, may overcount): {n:,}", flush=True)
            return

        # ── Step 1: SRA fetch ─────────────────────────────────────────────────
        print("\n── Step 1: SRA run fetch ────────────────────────────────────",
              flush=True)
        run_cache = load_json(run_path)
        uid_cache = load_json(uid_path)
        print(f"Cached: {len(run_cache):,} runs  |  {len(uid_cache):,} UIDs fetched",
              flush=True)

        uids  = list(set(fetch_uids(seeds)))
        added = fetch_runinfo(uids, run_cache, uid_cache, run_path, uid_path)
        print(f"SRA: {len(run_cache):,} total runs  ({added:,} new)", flush=True)

        # ── Step 2: STAT fetch ────────────────────────────────────────────────
        print("\n── Step 2: STAT taxonomy fetch ──────────────────────────────",
              flush=True)
        accessions = list(run_cache.keys())
        fetched    = _load_index(stat_index)
        to_fetch   = [a for a in accessions if a not in fetched]
        n_cached   = len(accessions) - len(to_fetch)
        print(f"STAT cache: {len(fetched):,} total entries, "
              f"{n_cached:,} for this mode already cached", flush=True)

        fetch_stat_append(to_fetch, stat_jsonl, stat_index)

        summary = (
            f"── 01_fetch_runs {args.mode.upper()} summary ──────────────────────────\n"
            f"Mode:              {args.mode.upper()}\n"
            f"\n"
            f"Unique UIDs found: {len(uids):>8,}\n"
            f"Runs fetched:      {len(run_cache):>8,}  ({added:,} new)\n"
            f"STAT newly fetched:{len(to_fetch):>8,}\n"
            f"\n"
            f"Runs: {run_path}\n"
            f"STAT: {stat_jsonl}\n"
            f"Log:  {log_dir / f'{args.mode}.log'}\n"
        )
        summary_path = log_dir / f"{args.mode}_summary.txt"
        summary_path.write_text(summary)
        link_latest(logs_base, summary_path)
        print(f"\n{summary}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
