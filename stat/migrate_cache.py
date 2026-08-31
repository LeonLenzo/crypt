#!/usr/bin/env python3
"""
stat/migrate_cache.py — migrate stat_cache.jsonl from old array format to new object format.

Old format (per line):  {acc}\t[{stat_array}]
New format (per line):  {acc}\t{"mode":..., "BioSample":..., ..., "_stat":[...]}

Sources for RunInfo:
  stat/output/stat_fetch/data/{mal,hal,aerial}_runs.json  → acc → RunInfo fields
  stat/output/stat_fetch/data/{mal,hal}_uids.json         → uid → bool (done set)

Outputs (overwrites in place):
  stat/output/stat_fetch/data/stat_cache.jsonl    — rewritten to new format
  stat/output/stat_fetch/data/mal_accessions.txt  — one accession per line
  stat/output/stat_fetch/data/hal_accessions.txt  — one accession per line
  stat/output/stat_fetch/data/mal_uid_checkpoint.txt  — one UID per line
  stat/output/stat_fetch/data/hal_uid_checkpoint.txt  — one UID per line

Aerial entries are converted to mode="aerial" but no aerial_accessions.txt is written
(aerial mode is dropped; stat_filter.py will not process these runs).

Run from crypt/:
    python stat/migrate_cache.py
"""

import json
import os
import sys
from pathlib import Path

DATA = Path("stat/output/stat_fetch/data")
CACHE = DATA / "stat_cache.jsonl"
TMP   = DATA / "stat_cache.jsonl.tmp"

_RUNINFO_FIELDS = (
    "BioSample", "BioProject", "SRAStudy",
    "Platform", "LibrarySource", "ScientificName",
)


def load_runs(mode: str) -> dict[str, dict]:
    p = DATA / f"{mode}_runs.json"
    if not p.exists():
        print(f"  [warn] {p} not found — {mode} RunInfo unavailable", file=sys.stderr)
        return {}
    with open(p) as f:
        return json.load(f)


def load_uids(mode: str) -> set[str]:
    p = DATA / f"{mode}_uids.json"
    if not p.exists():
        return set()
    with open(p) as f:
        return set(json.load(f).keys())


def write_accessions(mode: str, runs: dict[str, dict]) -> None:
    p = DATA / f"{mode}_accessions.txt"
    with open(p, "w") as f:
        for acc in sorted(runs):
            f.write(acc + "\n")
    print(f"  wrote {len(runs):,} → {p}")


def write_uid_checkpoint(mode: str, uids: set[str]) -> None:
    p = DATA / f"{mode}_uid_checkpoint.txt"
    with open(p, "w") as f:
        for uid in sorted(uids, key=lambda u: int(u)):
            f.write(uid + "\n")
    print(f"  wrote {len(uids):,} UIDs → {p}")


def main() -> None:
    print("Loading RunInfo …")
    mal_runs  = load_runs("mal")
    hal_runs  = load_runs("hal")
    aer_runs  = load_runs("aerial")

    # acc → (mode, runinfo_dict)
    acc_map: dict[str, tuple[str, dict]] = {}
    for mode, runs in (("mal", mal_runs), ("hal", hal_runs), ("aerial", aer_runs)):
        for acc, ri in runs.items():
            acc_map[acc] = (mode, ri)

    print(f"  MAL {len(mal_runs):,}  HAL {len(hal_runs):,}  aerial {len(aer_runs):,}")

    print("Loading UIDs …")
    mal_uids = load_uids("mal")
    hal_uids = load_uids("hal")
    print(f"  MAL {len(mal_uids):,}  HAL {len(hal_uids):,}")

    print(f"Migrating {CACHE} …")
    n_old = n_new = n_skip = n_missing = 0

    with open(CACHE) as src, open(TMP, "w") as dst:
        for lineno, raw in enumerate(src, 1):
            raw = raw.rstrip("\n")
            if not raw:
                continue
            tab = raw.index("\t")
            acc  = raw[:tab]
            rest = raw[tab + 1:]

            if rest.startswith("{"):
                # Already new format — pass through unchanged.
                dst.write(raw + "\n")
                n_new += 1
                continue

            if rest == "null":
                # Old format: STAT returned no data for this run.
                stat_array = None
            elif rest.startswith("["):
                # Old format: rest is a JSON array (the STAT tax_totals list).
                stat_array = json.loads(rest)
            else:
                print(f"  [warn] line {lineno}: unrecognised format, skipping")
                n_skip += 1
                continue

            if acc not in acc_map:
                # Not in any runs.json — keep old format, can't enrich.
                dst.write(raw + "\n")
                n_missing += 1
                continue

            mode, ri = acc_map[acc]
            entry = {"mode": mode}
            for field in _RUNINFO_FIELDS:
                entry[field] = ri.get(field, "")
            entry["_stat"] = stat_array

            dst.write(acc + "\t" + json.dumps(entry, separators=(",", ":")) + "\n")
            n_old += 1

            if lineno % 50_000 == 0:
                print(f"  … {lineno:,} lines processed", flush=True)

    print(f"  converted {n_old:,} old  |  kept {n_new:,} new  |  "
          f"no-RunInfo {n_missing:,}  |  skipped {n_skip:,}")

    # Atomically replace original.
    os.replace(TMP, CACHE)
    print(f"  stat_cache.jsonl updated in place")

    print("Writing accessions files …")
    write_accessions("mal", mal_runs)
    write_accessions("hal", hal_runs)

    print("Writing UID checkpoint files …")
    write_uid_checkpoint("mal", mal_uids)
    write_uid_checkpoint("hal", hal_uids)

    print("Done.")


if __name__ == "__main__":
    main()
