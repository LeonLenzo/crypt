#!/usr/bin/env python3
"""
kraken_run_select.py — select target BioSamples and download reads.
Submodule 2, step 1 of 3 (select -> split -> assign).

Host reference genome resolution/download moved to kraken_run_split.py
(2026-09-07 refactor — see kraken/README.md's "select vs split ownership"
note). Originally built with genome download here too, matching
kraken_db_search.py's "the ONE place that fetches" pattern — but that made
`select` a heavy, slow script doing multi-GB network I/O as a side effect of
what's supposed to be a cheap filtering step, and meant re-running it under
a different --setting/--limit risked re-triggering large downloads. `split`
already needs the host genome on disk to build its BBMap index, so it's the
more natural owner of "make sure the genome exists" — `select` now only
resolves and records which taxids are needed, it doesn't fetch anything for
them itself.

  1. SELECT — filter metadata/output/meta_classify/data/samples.tsv down to the
     target BioSamples: default is every field, aerial-tissue BioSample
     (2,719 as of 2026-08-28) — NOT narrowed to the subset STAT/LLM already
     flagged as cryptic-co-infected, since that would defeat the point of a
     validation pass meant to catch what STAT missed (the documented PST/rust
     blind spot especially shows up as "clean" per STAT, not flagged).
     --cryptic-only narrows to the already-flagged subset for a quick smoke
     test; --setting broadens/narrows beyond field. Then join against
     stat/output/stat_filter/data/runs.tsv on BioSample to
     get every Run accession belonging to each target sample. Host taxid comes
     straight from samples.tsv's llm_host_resolved_taxid column — resolved by
     meta_classify.py itself (deterministic name->taxid lookup right where the
     name was extracted; see _util.resolve_taxon_name) — this script just reads
     it, no resolution logic here. Host, not pathogen, and from the manuscript
     (LLM), not STAT — STAT's inferred host can be wrong/generic; the
     author-stated host is the ground truth for what to remove as background.
     Every named candidate taxid (not just the resolved one) is recorded per
     row in candidate_host_taxids — kraken_run_split.py needs the full set,
     since an unresolved multi-host sample's true host is still covered as
     long as it's one of the named candidates.

  2. DOWNLOAD — prefetch + fasterq-dump each Run's full reads (no subsampling
     — scratch space isn't the constraint, download time is, and fasterq-dump
     is the fastest path to full files).

Run on Setonix (requires sra-tools: prefetch, fasterq-dump):
    python kraken/run/kraken_run_select.py
    python kraken/run/kraken_run_select.py --setting field,greenhouse --limit 50   # broader/smaller test

Output:
    kraken/output/run/select/data/run_list.tsv   (tracked — Run/BioSample/host/
        candidate_host_taxids/status; kraken_run_split.py consumes this to
        resolve+download host genomes itself)
    kraken/output/run/select/data/reads/{run}_1.fastq.gz [+ _2]  (gitignored)
"""

import argparse
import csv
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest

SAMPLES_TSV = Path("metadata/output/meta_classify/data/samples.tsv")
RUNS_TSV    = Path("stat/output/stat_filter/data/runs.tsv")

OUT_DIR   = Path("kraken/output/run/select")
DATA_DIR  = OUT_DIR / "data"
RUN_LIST  = DATA_DIR / "run_list.tsv"
READS_DIR = DATA_DIR / "reads"

RUN_LIST_COLS = [
    "Run", "BioSample", "BioProject", "llm_host_resolved", "host_taxid",
    "candidate_host_taxids", "download_status",
]

_CRYPTIC_PMS = {"partial_match_plus_undeclared", "no_match_stat_found_different",
                "stat_only_no_named"}


# ── filter (mirrors metadata/figures/sample_funnel_v3.py) ────────────────────

def ti(r):
    t = r.get("llm_tissue", "")
    return t if t in ("aerial", "non-aerial") else "unclear"


def se(r):
    s = r.get("llm_study_setting", "")
    return s if s in ("field", "greenhouse", "growth_chamber",
                      "detached_leaf_assay", "in_vitro") else "unclear"


def is_cryptic(r) -> bool:
    stress = r.get("llm_stress", "")
    if stress == "biotic":
        return r.get("pathogen_match_status", "") in _CRYPTIC_PMS
    if stress in ("abiotic", "none"):
        try:
            return int(r.get("n_pathogens") or 0) > 0
        except ValueError:
            return False
    return False


def select_targets(rows: list, settings: set, require_cryptic: bool,
                   aerial_only: bool) -> list:
    out = []
    for r in rows:
        if aerial_only and ti(r) != "aerial":
            continue
        if se(r) not in settings:
            continue
        if require_cryptic and not is_cryptic(r):
            continue
        out.append(r)
    return out


# ── read download (prefetch + fasterq-dump, full files) ──────────────────────

def download_run(run: str, dest_dir: Path) -> str:
    """prefetch + fasterq-dump a run's full reads to dest_dir, gzipped.
    Returns 'ok', 'cached', or 'failed'. Resumable: skips if gz output exists."""
    r1 = dest_dir / f"{run}_1.fastq.gz"
    r2 = dest_dir / f"{run}_2.fastq.gz"
    se_ = dest_dir / f"{run}.fastq.gz"
    if (r1.exists() and r2.exists()) or se_.exists():
        return "cached"

    dest_dir.mkdir(parents=True, exist_ok=True)
    sra_dir = dest_dir / f"_{run}_sra"
    try:
        r = subprocess.run(
            ["prefetch", "--max-size", "50g", "-O", str(sra_dir), run],
            capture_output=True, text=True, timeout=3600,
        )
        if r.returncode != 0:
            return "failed"
        sra_file = sra_dir / run / f"{run}.sra"
        if not sra_file.exists():
            candidates = list(sra_dir.glob(f"**/{run}.sra"))
            if not candidates:
                return "failed"
            sra_file = candidates[0]

        r = subprocess.run(
            ["fasterq-dump", "--split-files", "-O", str(dest_dir), str(sra_file)],
            capture_output=True, text=True, timeout=7200,
        )
        if r.returncode != 0:
            return "failed"

        raw_r1 = dest_dir / f"{run}_1.fastq"
        raw_r2 = dest_dir / f"{run}_2.fastq"
        raw_se = dest_dir / f"{run}.fastq"
        ok = False
        if raw_r1.exists() and raw_r2.exists():
            for raw, gz in [(raw_r1, r1), (raw_r2, r2)]:
                subprocess.run(["gzip", "-f", str(raw)], timeout=1800)
            ok = r1.exists() and r2.exists()
        elif raw_se.exists():
            subprocess.run(["gzip", "-f", str(raw_se)], timeout=1800)
            ok = se_.exists()
        return "ok" if ok else "failed"
    finally:
        if sra_dir.exists():
            subprocess.run(["rm", "-rf", str(sra_dir)])


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--setting", default="field",
                    help="Comma-separated llm_study_setting values to include "
                         "(default: field). Use 'all' for every setting.")
    ap.add_argument("--aerial-only", dest="aerial_only", action="store_true", default=True)
    ap.add_argument("--no-aerial-only", dest="aerial_only", action="store_false",
                    help="Include non-aerial tissue too (default: aerial only)")
    ap.add_argument("--cryptic-only", dest="require_cryptic", action="store_true", default=False,
                    help="Narrow to only the samples STAT/LLM already flagged as "
                         "cryptic-co-infected, instead of every field/aerial sample "
                         "(default: off — restricting to already-flagged samples defeats "
                         "the point of a validation pass meant to catch what STAT missed, "
                         "e.g. the documented PST/rust blind spot; useful for a quick "
                         "smoke test only)")
    ap.add_argument("--limit", type=int, default=None,
                    help="Process at most N target BioSamples (for testing). "
                         "NOTE: samples.tsv is not randomly ordered — rows from "
                         "the same BioProject sit together (confirmed 2026-09-08: "
                         "the first 100 field/aerial rows are ALL one wheat "
                         "BioProject). Pair with --shuffle for a real "
                         "cross-cohort sample instead of --limit N's plain "
                         "file-order prefix.")
    ap.add_argument("--shuffle", action="store_true",
                    help="Shuffle target BioSamples (fixed --seed, reproducible) "
                         "before applying --limit, so a small --limit batch spans "
                         "many BioProjects/hosts instead of just the first ones "
                         "in samples.tsv's file order.")
    ap.add_argument("--seed", type=int, default=42,
                    help="Random seed for --shuffle (default 42, for reproducibility)")
    ap.add_argument("--download", action="store_true",
                    help="After selecting + resolving hosts, download reads")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_run_select.log")
    link_latest(logs_base, log_dir / "kraken_run_select.log")
    sys.stdout = log

    try:
        settings = ({"field", "greenhouse", "growth_chamber", "detached_leaf_assay", "in_vitro"}
                    if args.setting == "all" else set(args.setting.split(",")))

        samples = list(csv.DictReader(open(SAMPLES_TSV), delimiter="\t"))
        targets = select_targets(samples, settings, args.require_cryptic, args.aerial_only)
        if args.shuffle:
            random.Random(args.seed).shuffle(targets)
        if args.limit:
            targets = targets[:args.limit]
        print(f"Target BioSamples: {len(targets)} (settings={sorted(settings)}, "
              f"cryptic_only={args.require_cryptic}, aerial_only={args.aerial_only})")

        # BioSample -> Run accessions
        bs_to_runs = {}
        with open(RUNS_TSV, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                bs_to_runs.setdefault(row["BioSample"], []).append(row["Run"])

        # kraken_run_split.py builds ONE bbmap.sh index PER distinct host taxid
        # named here (not one combined multi-reference index — some host
        # genomes are huge, e.g. wheat 14.5Gb, pine 22.4Gb; combining ~90 of them
        # into a single BBSplit index would be several hundred Gb and unbuildable
        # on any Setonix node), and it resolves+downloads each one itself (see
        # module docstring — moved out of this script 2026-09-07). A sample
        # doesn't need its own host confidently resolved to still get correctly
        # host-filtered: kraken_run_split.py runs each candidate's individual
        # index separately and picks the one with the highest mapping rate as
        # the confirmed host — an independent, read-level confirmation (often
        # more reliable than metadata guessing). So: record every CANDIDATE in
        # llm_named_hosts_taxids, not just the single llm_host_resolved_taxid —
        # an unresolved multi-host sample's true host is still covered as long
        # as it's one of the named candidates. host_taxid (singular, resolved)
        # is kept per-row for provenance/QC only.
        run_rows = []
        host_taxids_needed = set()
        n_resolved, n_unresolved = 0, 0
        for r in targets:
            bs = r["BioSample"]
            host_resolved = r.get("llm_host_resolved", "")
            resolved_taxid = (r.get("llm_host_resolved_taxid", "") or "").strip()
            candidate_taxids = sorted({t for t in
                (x.strip() for x in (r.get("llm_named_hosts_taxids", "") or "").split(";"))
                if t})
            if resolved_taxid:
                n_resolved += 1
            else:
                n_unresolved += 1
            host_taxids_needed.update(int(t) for t in candidate_taxids)
            for run in bs_to_runs.get(bs, []):
                run_rows.append({
                    "Run": run, "BioSample": bs, "BioProject": r.get("BioProject", ""),
                    "llm_host_resolved": host_resolved,
                    "host_taxid": resolved_taxid,
                    "candidate_host_taxids": "; ".join(candidate_taxids),
                    "download_status": "",
                })

        print(f"Host resolution: {n_resolved} confidently resolved, {n_unresolved} "
              f"not per-sample-resolved (still covered by per-candidate indices "
              f"below as long as their true host is a named candidate)")
        print(f"Per-taxid host indices needed: {len(host_taxids_needed)} distinct taxids "
              f"(resolved+downloaded by kraken_run_split.py, not here)")
        print(f"Runs to process: {len(run_rows)}")

        # ── download reads ─────────────────────────────────────────────────────
        if args.download:
            READS_DIR.mkdir(parents=True, exist_ok=True)
            print(f"\nDownloading reads for {len(run_rows)} runs …")
            t0 = time.time()
            n_ok = n_cached = n_fail = 0
            # Runs are unique per row here (one row per Run), so this index is 1:1 —
            # avoids an O(n) scan of run_rows per completion (O(n^2) overall).
            rows_by_run = {row["Run"]: row for row in run_rows}

            def work(row):
                return row["Run"], download_run(row["Run"], READS_DIR)

            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                futs = {pool.submit(work, row): row["Run"] for row in run_rows}
                for done, fut in enumerate(as_completed(futs), 1):
                    run, status = fut.result()
                    rows_by_run[run]["download_status"] = status
                    if status == "ok":
                        n_ok += 1
                    elif status == "cached":
                        n_cached += 1
                    else:
                        n_fail += 1
                    if done % 10 == 0 or done == len(run_rows):
                        elapsed = time.time() - t0
                        rate = done / elapsed if elapsed > 0 else 0
                        print(f"  [{done}/{len(run_rows)}] ok={n_ok} cached={n_cached} "
                              f"fail={n_fail}  rate={rate:.2f}/s", flush=True)

        # ── write run list ─────────────────────────────────────────────────────
        with open(RUN_LIST, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=RUN_LIST_COLS, delimiter="\t")
            w.writeheader()
            w.writerows(run_rows)
        print(f"\nOutput: {RUN_LIST}")
        if not args.download:
            print("Re-run with --download to fetch reads.")
        else:
            print("Next: python kraken/run/kraken_run_split.py "
                  "(resolves + downloads host genomes, builds indices, splits)")
    finally:
        log.close()


if __name__ == "__main__":
    main()
