#!/usr/bin/env python3
"""
01_sra.py — fetch SRA run IDs for cryptic co-infection mining.

Two modes:
  mal  microbe-as-library: PHI-base plant pathogen as library organism.
  hal  host-as-library: PHI-base plant host as library organism.

Output (output/01_sra/):
  {mode}_runs.json   RunInfo rows keyed by Run accession
  {mode}_uids.json   fetched UID set (for resumability)
  {mode}.log         full run log

Usage:
  python 01_sra.py --mode mal
  python 01_sra.py --mode hal
  python 01_sra.py --mode mal --count   (hit counts only, no fetch)
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

# ── Settings ──────────────────────────────────────────────────────────────────

DB_PATH       = Path("output/00_build/phibase_db.json")
OUT_DIR       = Path("output/01_sra")
LIBRARY_STRAT = "RNA-Seq"
MAL_BATCH     = 50    # seed taxids per esearch query
RUNINFO_BATCH = 500   # UIDs per RunInfo efetch
MAX_WORKERS   = 8     # parallel RunInfo fetch workers
SAVE_EVERY    = 20    # save caches after this many batches

API_KEY = os.environ.get("NCBI_API_KEY", "")
RATE    = 9.0 if API_KEY else 2.5   # requests/second
HEADERS = {"User-Agent": "crypt/01_sra (leon.lenzo@curtin.edu.au)"}
ENTREZ  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


# ── Logging ───────────────────────────────────────────────────────────────────

class _Tee:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._log    = open(path, "w", buffering=1)
        self._stdout = sys.__stdout__

    def write(self, text: str) -> None:
        self._stdout.write(text)
        self._log.write(text)

    def flush(self) -> None:
        self._stdout.flush()
        self._log.flush()

    def close(self) -> None:
        sys.stdout = self._stdout
        self._log.close()


# ── HTTP ──────────────────────────────────────────────────────────────────────

def _get(url: str, retries: int = 5) -> bytes:
    delay = 1.0
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(delay); delay *= 2
            else:
                raise
        except Exception:
            time.sleep(delay); delay *= 2
    raise RuntimeError(f"Failed after {retries} retries: {url}")


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
    """Return hit count for term without fetching any UIDs."""
    url  = _url("esearch.fcgi", db=db, term=term, retmax=0, retmode="json")
    data = json.loads(_get(url))["esearchresult"]
    time.sleep(1.0 / RATE)
    return int(data["count"])


def _esearch(db: str, term: str, retmax: int = 10_000) -> list[str]:
    """Return all UIDs matching term, paginating via WebEnv."""
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


# ── Cache helpers ─────────────────────────────────────────────────────────────

def _load(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def _save(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


# ── Query functions ───────────────────────────────────────────────────────────

def _load_db() -> dict:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"PHI-base DB not found: {DB_PATH}\nRun:  python3 00_build.py")
    with open(DB_PATH) as f:
        return json.load(f)


def _load_seed_pathogens() -> list[int]:
    raw   = _load_db()
    seeds = sorted(int(k) for k in raw["pathogen_to_hosts"].keys())
    print(f"PHI-base: {len(seeds)} seed pathogen species", flush=True)
    return seeds


def _load_seed_hosts() -> list[int]:
    raw   = _load_db()
    seeds = sorted(int(k) for k in raw["host_to_pathogens"].keys())
    print(f"PHI-base: {len(seeds)} seed host species", flush=True)
    return seeds


def mal_count() -> int:
    """MAL: sum esearch hit counts across batches (may overcount cross-batch duplicates)."""
    seeds     = _load_seed_pathogens()
    n_batches = (len(seeds) + MAL_BATCH - 1) // MAL_BATCH
    total     = 0

    for i in range(0, len(seeds), MAL_BATCH):
        batch    = seeds[i:i + MAL_BATCH]
        or_terms = " OR ".join(f"txid{t}[Organism:exp]" for t in batch)
        query    = f"({or_terms}) AND {LIBRARY_STRAT}[Strategy]"
        n     = _esearch_count("sra", query)
        total += n
        bn    = i // MAL_BATCH + 1
        print(f"  [{bn}/{n_batches}] {n:,} hits  (running total {total:,})",
              flush=True)

    print(f"\nMAL total (batched, may overcount): {total:,}", flush=True)
    return total


def mal_uids() -> list[str]:
    """MAL: query each batch of PHI-base seed pathogen taxids × plant [Host]."""
    seeds     = _load_seed_pathogens()
    n_batches = (len(seeds) + MAL_BATCH - 1) // MAL_BATCH
    all_uids: list[str] = []

    for i in range(0, len(seeds), MAL_BATCH):
        batch    = seeds[i:i + MAL_BATCH]
        or_terms = " OR ".join(f"txid{t}[Organism:exp]" for t in batch)
        query    = f"({or_terms}) AND {LIBRARY_STRAT}[Strategy]"
        uids     = _esearch("sra", query)
        all_uids.extend(uids)
        bn       = i // MAL_BATCH + 1
        print(f"  [{bn}/{n_batches}] +{len(uids):,} UIDs  "
              f"(total {len(all_uids):,})", flush=True)

    return all_uids


def hal_count() -> int:
    """HAL: sum esearch hit counts for PHI-base plant host taxids."""
    seeds     = _load_seed_hosts()
    n_batches = (len(seeds) + MAL_BATCH - 1) // MAL_BATCH
    total     = 0

    for i in range(0, len(seeds), MAL_BATCH):
        batch    = seeds[i:i + MAL_BATCH]
        or_terms = " OR ".join(f"txid{t}[Organism:exp]" for t in batch)
        query    = f"({or_terms}) AND {LIBRARY_STRAT}[Strategy]"
        n        = _esearch_count("sra", query)
        total   += n
        bn       = i // MAL_BATCH + 1
        print(f"  [{bn}/{n_batches}] {n:,} hits  (running total {total:,})",
              flush=True)

    print(f"\nHAL total (batched, may overcount): {total:,}", flush=True)
    return total


def hal_uids() -> list[str]:
    """HAL: PHI-base plant host taxids as library organism."""
    seeds     = _load_seed_hosts()
    n_batches = (len(seeds) + MAL_BATCH - 1) // MAL_BATCH
    all_uids: list[str] = []

    for i in range(0, len(seeds), MAL_BATCH):
        batch    = seeds[i:i + MAL_BATCH]
        or_terms = " OR ".join(f"txid{t}[Organism:exp]" for t in batch)
        query    = f"({or_terms}) AND {LIBRARY_STRAT}[Strategy]"
        uids     = _esearch("sra", query)
        all_uids.extend(uids)
        bn       = i // MAL_BATCH + 1
        print(f"  [{bn}/{n_batches}] +{len(uids):,} UIDs  "
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
    """Fetch RunInfo in parallel for UIDs not already in uid_cache."""
    to_fetch = [u for u in uids if u not in uid_cache]
    if not to_fetch:
        print(f"  All {len(uids):,} UIDs already fetched.", flush=True)
        return 0

    batches = [to_fetch[i:i + RUNINFO_BATCH]
               for i in range(0, len(to_fetch), RUNINFO_BATCH)]
    print(f"  Fetching RunInfo for {len(to_fetch):,} UIDs in "
          f"{len(batches):,} batches ({len(uids) - len(to_fetch):,} cached) …",
          flush=True)

    added = 0
    done  = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
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
                _save(run_cache, run_path)
                _save(uid_cache, uid_path)
                print(f"  [{done:,}/{len(batches):,} batches]  "
                      f"{added:,} new runs", flush=True)

    return added


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True,
                    help="mal = microbe-as-library  |  hal = host-as-library")
    ap.add_argument("--count", action="store_true",
                    help="print hit counts only; do not fetch run IDs or RunInfo")
    args = ap.parse_args()

    log = _Tee(OUT_DIR / f"{args.mode}.log")
    sys.stdout = log

    try:
        print(f"Mode: {args.mode.upper()}", flush=True)

        if args.count:
            mal_count() if args.mode == "mal" else hal_count()
            print(f"\nLog → {OUT_DIR / f'{args.mode}.log'}")
            return

        run_path  = OUT_DIR / f"{args.mode}_runs.json"
        uid_path  = OUT_DIR / f"{args.mode}_uids.json"
        run_cache = _load(run_path)
        uid_cache = _load(uid_path)
        print(f"Cached: {len(run_cache):,} runs  |  {len(uid_cache):,} UIDs fetched",
              flush=True)

        uids = mal_uids() if args.mode == "mal" else hal_uids()
        uids = list(set(uids))
        print(f"\nTotal unique UIDs: {len(uids):,}", flush=True)

        added = fetch_runinfo(uids, run_cache, uid_cache, run_path, uid_path)

        summary = (
            f"── 01_sra {args.mode.upper()} summary ──────────────────────────\n"
            f"Mode:              {args.mode.upper()}\n"
            f"\n"
            f"Unique UIDs found: {len(uids):>8,}\n"
            f"New runs fetched:  {added:>8,}\n"
            f"Total runs cached: {len(run_cache):>8,}\n"
            f"\n"
            f"Runs: {run_path}\n"
            f"Log:  {OUT_DIR / f'{args.mode}.log'}\n"
        )
        (OUT_DIR / f"{args.mode}_summary.txt").write_text(summary)
        print(f"\n{summary}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
