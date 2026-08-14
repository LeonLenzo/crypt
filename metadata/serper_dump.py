#!/usr/bin/env python3
"""
serper_dump.py — query Serper (Google) for no-PMID BioProjects and dump raw results.

Saves raw search results to metadata/output/serper/serper_results.tsv for manual
inspection before committing to automated PMID extraction.

Usage:
    python metadata/serper_dump.py              # all no-PMID BPs
    python metadata/serper_dump.py --limit 50   # first N for a quick look
"""

import argparse
import csv
import json
import os
import time
import urllib.request
from pathlib import Path

LIT_JSON    = Path("metadata/output/fetch_lit/data/literature.json")
CACHE_PATH  = Path("metadata/output/fetch_lit/data/lit_cache.json")
OUT_DIR     = Path("metadata/output/serper")
OUT_TSV     = OUT_DIR / "serper_results.tsv"

SERPER_API_KEY = os.environ.get("SERPER_API_KEY", "")

OUT_COLS = ["bioproject", "rank", "title", "link", "snippet"]

_last: float = 0.0

def _wait() -> None:
    global _last
    gap  = 1.0
    wait = _last + gap - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last = time.monotonic()


def serper_query(bp: str, num: int = 10) -> list[dict]:
    _wait()
    body = json.dumps({"q": f'"{bp}"', "num": num}).encode()
    req  = urllib.request.Request(
        "https://google.serper.dev/search", data=body,
        headers={"X-API-KEY": SERPER_API_KEY, "Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            data = json.loads(r.read())
        return data.get("organic", [])
    except Exception as e:
        print(f"  WARNING: {bp} query failed: {e}", flush=True)
        return []


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=0, help="Max BPs to query (0 = all)")
    ap.add_argument("--num",   type=int, default=10, help="Results per query (max 10)")
    args = ap.parse_args()

    if not SERPER_API_KEY:
        raise SystemExit("SERPER_API_KEY not set")

    lit = json.loads(LIT_JSON.read_text())
    no_pmid = [bp for bp, v in lit.items() if not v.get("primary_pmid")]
    print(f"No-PMID BioProjects: {len(no_pmid)}", flush=True)

    if args.limit:
        no_pmid = no_pmid[:args.limit]
        print(f"Limiting to first {args.limit}", flush=True)

    # Skip any already in the dump
    existing: set[str] = set()
    if OUT_TSV.exists():
        with open(OUT_TSV) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                existing.add(row["bioproject"])
        print(f"Already dumped: {len(existing)} BPs — skipping", flush=True)

    to_query = [bp for bp in no_pmid if bp not in existing]
    print(f"To query: {len(to_query)}", flush=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_header = not OUT_TSV.exists()

    with open(OUT_TSV, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        if write_header:
            w.writeheader()

        for i, bp in enumerate(to_query, 1):
            print(f"[{i:>4}/{len(to_query)}] {bp}", end=" ", flush=True)
            results = serper_query(bp, num=args.num)
            print(f"→ {len(results)} results", flush=True)
            for rank, it in enumerate(results, 1):
                w.writerow({
                    "bioproject": bp,
                    "rank":       rank,
                    "title":      it.get("title", ""),
                    "link":       it.get("link", ""),
                    "snippet":    (it.get("snippet", "") or "").replace("\n", " "),
                })
            fh.flush()

            if not results:
                w.writerow({"bioproject": bp, "rank": 0,
                            "title": "", "link": "", "snippet": ""})
                fh.flush()

    n_with_results = sum(1 for bp in no_pmid if bp in existing or bp in set(to_query))
    print(f"\nDone. Written to {OUT_TSV}")
    print(f"Rows: {sum(1 for _ in open(OUT_TSV)) - 1}")


if __name__ == "__main__":
    main()
