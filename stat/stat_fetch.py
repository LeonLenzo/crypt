#!/usr/bin/env python3
"""
stat/stat_fetch.py — fetch SRA run IDs + NCBI STAT taxonomy.

Two modes:
  mal   Microbe-as-library: PHI-base plant pathogen as library organism.
  hal   Host-as-library:    PHI-base plant host as library organism.

No retention gate is applied here — that belongs in stat/stat_filter.py.

Output (stat/output/stat_fetch/data/):
  stat_cache.jsonl           Unified cache: RunInfo + STAT, one JSON line per run.
                             Line format: {acc}\t{json}
                             JSON keys: mode, BioSample, BioProject, SRAStudy,
                             Platform, LibrarySource, ScientificName, _stat
  stat_cache_index.txt       Accessions with STAT fetched — STAT resume checkpoint.
  {mode}_uid_checkpoint.txt  UIDs with RunInfo fetched — RunInfo resume checkpoint.
  {mode}_accessions.txt      All accessions in this mode (used by stat_filter.py).

Do NOT run MAL and HAL simultaneously — both write to stat_cache.jsonl.
Chain: finish MAL, then run HAL.

Run from crypt/:
    python stat/stat_fetch.py --mode mal
    python stat/stat_fetch.py --mode hal
    python stat/stat_fetch.py --mode mal --count
"""

import argparse
import csv
import io
import json
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, http_get, link_latest, make_log_dir

# ── Constants ─────────────────────────────────────────────────────────────────

DB_PATH       = Path("stat/output/stat_build/data/phibase_db.json")
OUT_DIR       = Path("stat/output/stat_fetch")

LIBRARY_STRAT = "RNA-Seq"
MAL_BATCH     = 50         # taxids per esearch OR-query
RUNINFO_BATCH = 500        # UIDs per efetch RunInfo request
SRA_WORKERS   = 8          # parallel RunInfo batch threads
STAT_WORKERS  = 32         # parallel STAT fetch threads
SAVE_EVERY    = 20         # RunInfo batches between checkpoint flushes
LOG_EVERY     = 200        # STAT fetches between progress lines
CACHE_MIN_PCT = 0.1        # prune tax_table entries below this % at write time

STAT_URL = ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/"
            "run_taxonomy?acc={run}&cluster_name=public")
ENTREZ   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

API_KEY   = os.environ.get("NCBI_API_KEY", "")
RATE      = 9.0  if API_KEY else 2.5
STAT_RATE = 15.0 if API_KEY else 5.0
HEADERS   = {"User-Agent": "crypt/stat_fetch (leon.lenzo@curtin.edu.au)"}

_RUNINFO_FIELDS = (
    "BioSample", "BioProject", "SRAStudy", "Platform",
    "LibrarySource", "ScientificName",
)


# ── HTTP / Entrez helpers ─────────────────────────────────────────────────────

def _get(url: str) -> bytes:
    return http_get(url, HEADERS)


def _post(endpoint: str, retries: int = 5, **params) -> bytes:
    import urllib.error, urllib.parse, urllib.request
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
    import urllib.parse
    if API_KEY:
        params["api_key"] = API_KEY
    return f"{ENTREZ}/{endpoint}?{urllib.parse.urlencode(params)}"


# ── esearch ───────────────────────────────────────────────────────────────────

def _esearch_count(term: str) -> int:
    data = json.loads(_get(_url("esearch.fcgi", db="sra", term=term,
                                retmax=0, retmode="json")))["esearchresult"]
    time.sleep(1.0 / RATE)
    return int(data["count"])


def _esearch_uids(term: str, retmax: int = 10_000) -> list[str]:
    data  = json.loads(_get(_url("esearch.fcgi", db="sra", term=term,
                                 usehistory="y", retmax=0,
                                 retmode="json")))["esearchresult"]
    total = int(data["count"])
    wk, qk = data["webenv"], data["querykey"]
    uids: list[str] = []
    for start in range(0, total, retmax):
        uids.extend(
            _get(_url("efetch.fcgi", db="sra", query_key=qk, WebEnv=wk,
                      retstart=start, retmax=retmax,
                      rettype="uilist", retmode="text")).decode().strip().splitlines()
        )
        time.sleep(1.0 / RATE)
    return uids


def fetch_uids(queries: list[str]) -> list[str]:
    all_uids: list[str] = []
    for i, query in enumerate(queries, 1):
        uids = _esearch_uids(query)
        all_uids.extend(uids)
        print(f"  [{i}/{len(queries)}] +{len(uids):,} UIDs  "
              f"(total {len(all_uids):,})", flush=True)
    return list(set(all_uids))


def count_hits(queries: list[str]) -> int:
    total = 0
    for i, query in enumerate(queries, 1):
        n = _esearch_count(query)
        total += n
        print(f"  [{i}/{len(queries)}] {n:,}  (running {total:,})", flush=True)
    return total


# ── RunInfo fetch ─────────────────────────────────────────────────────────────

_ri_lock = threading.Lock()
_ri_last: list[float] = [0.0]


def _fetch_runinfo_batch(uids: list[str]) -> tuple[list[str], list[dict]]:
    gap = 1.0 / RATE
    with _ri_lock:
        wait = _ri_last[0] + gap - time.time()
        if wait > 0:
            time.sleep(wait)
        _ri_last[0] = time.time()
    text = _post("efetch.fcgi", db="sra", id=",".join(uids),
                 rettype="runinfo", retmode="text").decode(errors="replace")
    rows = [r for r in csv.DictReader(io.StringIO(text)) if r.get("Run")]
    return uids, rows


def fetch_runinfo(
    all_uids: list[str],
    uid_done: set[str],
    uid_checkpoint: Path,
    accessions_file: Path,
) -> dict[str, dict]:
    """Fetch RunInfo for UIDs not in uid_done.
    Appends processed UIDs to uid_checkpoint and new accessions to accessions_file.
    Returns {acc: runinfo_dict} for newly fetched runs only."""
    to_fetch = [u for u in all_uids if u not in uid_done]
    if not to_fetch:
        print(f"  All {len(all_uids):,} UIDs already in checkpoint.", flush=True)
        return {}

    batches = [to_fetch[i:i + RUNINFO_BATCH]
               for i in range(0, len(to_fetch), RUNINFO_BATCH)]
    print(f"  Fetching RunInfo: {len(to_fetch):,} UIDs  "
          f"({len(all_uids) - len(to_fetch):,} cached)  "
          f"{len(batches):,} batches …", flush=True)

    new_meta: dict[str, dict] = {}
    done = 0

    with (open(uid_checkpoint, "a") as uid_f,
          open(accessions_file, "a") as acc_f,
          ThreadPoolExecutor(max_workers=SRA_WORKERS) as pool):
        futures = {pool.submit(_fetch_runinfo_batch, b): b for b in batches}
        for future in as_completed(futures):
            uids_batch, rows = future.result()
            for row in rows:
                acc = row["Run"]
                if acc not in new_meta:
                    new_meta[acc] = {k: row.get(k, "") for k in _RUNINFO_FIELDS}
                    acc_f.write(acc + "\n")
            for uid in uids_batch:
                uid_f.write(uid + "\n")
            done += 1
            if done % SAVE_EVERY == 0 or done == len(batches):
                uid_f.flush()
                acc_f.flush()
                print(f"  [{done:,}/{len(batches):,} batches]  "
                      f"{len(new_meta):,} new runs", flush=True)

    return new_meta


# ── STAT fetch ────────────────────────────────────────────────────────────────

def _filter_stat(result: list | None) -> list | None:
    """Prune tax_table entries below CACHE_MIN_PCT at write time."""
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


_stat_lock = threading.Lock()
_stat_last: list[float] = [0.0]


def _throttled_stat(run: str) -> tuple[str, list | None]:
    gap = 1.0 / STAT_RATE
    with _stat_lock:
        wait = _stat_last[0] + gap - time.time()
        if wait > 0:
            time.sleep(wait)
        _stat_last[0] = time.time()
    url = STAT_URL.format(run=run)
    try:
        body = http_get(url, HEADERS)
        return run, json.loads(body) if body else None
    except Exception as e:
        print(f"  STAT error {run}: {e}", flush=True)
        return run, None


def fetch_stat(
    run_meta: dict[str, dict],
    mode: str,
    stat_cache: Path,
    stat_index: Path,
    stat_done: set[str],
) -> int:
    """Fetch STAT for runs not in stat_done.
    Writes unified cache entries (RunInfo + STAT) to stat_cache.
    Returns number of newly fetched runs."""
    to_fetch = [acc for acc in run_meta if acc not in stat_done]
    if not to_fetch:
        print("  All runs already in STAT cache.", flush=True)
        return 0

    print(f"  Fetching STAT for {len(to_fetch):,} runs …", flush=True)
    done = 0

    with (open(stat_cache, "a") as jf,
          open(stat_index, "a") as ix,
          ThreadPoolExecutor(max_workers=STAT_WORKERS) as pool):
        futures = {pool.submit(_throttled_stat, acc): acc for acc in to_fetch}
        for future in as_completed(futures):
            run, result = future.result()
            entry = {"mode": mode, **run_meta[run], "_stat": _filter_stat(result)}
            jf.write(run + "\t" + json.dumps(entry) + "\n")
            jf.flush()
            ix.write(run + "\n")
            ix.flush()
            done += 1
            if done % LOG_EVERY == 0 or done == len(to_fetch):
                print(f"  [{done:,}/{len(to_fetch):,}]  "
                      f"{time.strftime('%H:%M:%S')}", flush=True)

    return done


# ── Seed + query helpers ──────────────────────────────────────────────────────

def _load_seeds(mode: str) -> list[int]:
    if not DB_PATH.exists():
        sys.exit(f"PHI-base DB not found: {DB_PATH}\nRun: python stat/stat_build.py")
    with open(DB_PATH) as f:
        raw = json.load(f)
    if mode == "mal":
        seeds = sorted(int(k) for k in raw["pathogen_to_hosts"])
        print(f"PHI-base: {len(seeds)} seed pathogen species", flush=True)
    else:
        seeds = sorted(set(raw["host_to_seed"].values()))
        print(f"PHI-base: {len(seeds)} seed host species", flush=True)
    return seeds


def _build_queries(seeds: list[int]) -> list[str]:
    queries = []
    for i in range(0, len(seeds), MAL_BATCH):
        batch    = seeds[i:i + MAL_BATCH]
        or_terms = " OR ".join(f"txid{t}[Organism:exp]" for t in batch)
        queries.append(f"({or_terms}) AND {LIBRARY_STRAT}[Strategy]")
    return queries


# ── Resume helpers ────────────────────────────────────────────────────────────

def _load_text_set(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return set(path.read_text().split())


def _load_runmeta_from_cache(cache_path: Path, mode_accs: set[str]) -> dict[str, dict]:
    """Extract RunInfo from new-format cache entries (have 'BioSample' key).
    Only loads entries whose accession is in mode_accs to limit memory use."""
    result: dict[str, dict] = {}
    if not cache_path.exists() or not mode_accs:
        return result
    with open(cache_path) as f:
        for line in f:
            if "\t" not in line:
                continue
            acc, _, rest = line.partition("\t")
            if acc not in mode_accs:
                continue
            try:
                data = json.loads(rest)
                if "BioSample" in data:   # new-format entry
                    result[acc] = {k: data.get(k, "") for k in _RUNINFO_FIELDS}
            except json.JSONDecodeError:
                pass
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["mal", "hal"], required=True)
    ap.add_argument("--count", action="store_true",
                    help="Print hit counts only; do not fetch.")
    args = ap.parse_args()
    mode = args.mode

    (OUT_DIR / "data").mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / f"{mode}.log")
    link_latest(logs_base, log_dir / f"{mode}.log")
    sys.stdout = log

    data_dir       = OUT_DIR / "data"
    stat_cache     = data_dir / "stat_cache.jsonl"
    stat_index     = data_dir / "stat_cache_index.txt"
    uid_checkpoint = data_dir / f"{mode}_uid_checkpoint.txt"
    accs_file      = data_dir / f"{mode}_accessions.txt"

    try:
        print(f"Mode: {mode.upper()}", flush=True)
        print(f"NCBI_API_KEY: {'set' if API_KEY else 'not set'}", flush=True)

        seeds   = _load_seeds(mode)
        queries = _build_queries(seeds)

        if args.count:
            n = count_hits(queries)
            print(f"\n{mode.upper()} total: {n:,}", flush=True)
            return

        # ── Resume state ──────────────────────────────────────────────────────
        print("\n── Resume state ─────────────────────────────────────────────",
              flush=True)
        stat_done    = _load_text_set(stat_index)
        uid_done     = _load_text_set(uid_checkpoint)
        mode_accs    = _load_text_set(accs_file)
        print(f"  STAT cache:     {len(stat_done):,} accessions", flush=True)
        print(f"  UID checkpoint: {len(uid_done):,} UIDs", flush=True)
        print(f"  Mode accessions:{len(mode_accs):,}", flush=True)

        run_meta = _load_runmeta_from_cache(stat_cache, mode_accs)
        print(f"  RunInfo loaded: {len(run_meta):,} runs", flush=True)

        # ── Step 1: Fetch UIDs ────────────────────────────────────────────────
        print("\n── Step 1: Fetch UIDs ───────────────────────────────────────",
              flush=True)
        all_uids = fetch_uids(queries)
        print(f"Unique UIDs: {len(all_uids):,}", flush=True)

        # ── Step 2: Fetch RunInfo ─────────────────────────────────────────────
        print("\n── Step 2: Fetch RunInfo ────────────────────────────────────",
              flush=True)
        new_meta = fetch_runinfo(all_uids, uid_done, uid_checkpoint, accs_file)
        run_meta.update(new_meta)
        print(f"Total RunInfo: {len(run_meta):,} runs", flush=True)

        # ── Step 3: Fetch STAT ────────────────────────────────────────────────
        print("\n── Step 3: Fetch STAT taxonomy ──────────────────────────────",
              flush=True)
        n_new = fetch_stat(run_meta, mode, stat_cache, stat_index, stat_done)

        # ── Summary ───────────────────────────────────────────────────────────
        summary = (
            f"── stat_fetch {mode.upper()} summary ──────────────────────────\n"
            f"Mode:               {mode.upper()}\n"
            f"\n"
            f"Unique UIDs:        {len(all_uids):>8,}\n"
            f"RunInfo in memory:  {len(run_meta):>8,}\n"
            f"STAT newly fetched: {n_new:>8,}\n"
            f"\n"
            f"Cache:    {stat_cache}\n"
            f"Mode accs:{accs_file}\n"
            f"Log:      {log_dir / f'{mode}.log'}\n"
        )
        summary_path = log_dir / f"{mode}_summary.txt"
        summary_path.write_text(summary)
        link_latest(logs_base, summary_path)
        print(f"\n{summary}")

    finally:
        log.close()


if __name__ == "__main__":
    main()
