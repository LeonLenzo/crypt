#!/usr/bin/env python3
"""
kraken_run_split.py — host-read removal via per-taxid BBMap indices.
Submodule 2, step 2 of 3 (select -> split -> assign).

Design (agreed 2026-09-01, see kraken/README.md's "kraken_run_split.py design"
section for the full rationale): ONE bbmap.sh index PER distinct host taxid, not
one combined multi-reference BBSplit index. The 2,719-sample field/aerial cohort
names 116 candidate host taxids; 94 have an NCBI genomic assembly totalling
~275Gb, dominated by outliers (Pinus radiata 22.4Gb, Allium sativum 16.5Gb,
Triticum aestivum 14.5Gb, Avena insularis 14.2Gb, ...) that can't just be
excluded — they're the dominant crops in the cohort, not obscure edge cases.
A single combined index at that scale isn't buildable on any Setonix node;
indexing each genome separately is trivial even for the largest one alone.

Two stages, each independently resumable:

  1. BUILD INDEX — one bbmap.sh index per distinct host taxid named in
     kraken_run_select.py's host_taxid_to_accession.json. The genomic FASTA is
     already extracted (kraken_db_search.download_cds() unzips on download) —
     no separate "extract" stage needed. Skips any taxid whose index already
     exists. Embarrassingly parallel across taxids, but kept to a modest
     default worker count (--workers) since a handful of these genomes are
     tens of Gb and each bbmap.sh index build gets its own -Xmx allocation —
     too much concurrency risks overcommitting node memory.

  2. SPLIT — for each run in run_list.tsv, align its reads against EVERY named
     candidate's index separately (not one combined pass). The candidate with
     the highest mapped fraction is the confirmed host — an independent,
     read-level confirmation, often more reliable than metadata alone. Its
     unmapped reads are the pathogen-enriched output for kraken_run_assign.py.
     Other candidates' mapped fractions are kept as QC/confirmation metadata
     only — no cross-candidate intersection of unmapped reads. Mapped fraction
     is computed directly from read counts (input reads vs. unmapped-output
     reads), not parsed from BBMap's statsfile text, so it doesn't depend on
     guessing an exact key format.

Run on Setonix (requires the bbmap module):
    module load bbmap/38.96--h5c4e2a8_0
    python kraken/run/kraken_run_split.py --build-index      # stage 1 only
    python kraken/run/kraken_run_split.py                    # both stages
    python kraken/run/kraken_run_split.py --limit 5           # smoke test

Output:
    kraken/output/run/split/data/index/{taxid}/               (gitignored —
        one BBMap index dir per host taxid)
    kraken/output/run/split/data/reads/{run}_{1,2}.fastq.gz    (gitignored —
        confirmed-host-removed, pathogen-enriched reads)
    kraken/output/run/split/data/split_results.tsv             (tracked — one
        row per run: confirmed host taxid + mapped %, every candidate's mapped
        % for QC, agreement with meta_classify.py's llm_host_resolved)
"""

import argparse
import csv
import gzip
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _util import _Tee, make_log_dir, link_latest, load_json, save_json

# ── paths ──────────────────────────────────────────────────────────────────────

SELECT_DATA     = Path("kraken/output/run/select/data")
RUN_LIST        = SELECT_DATA / "run_list.tsv"
SELECT_READS_DIR = SELECT_DATA / "reads"
HOST_TAXID_MAP  = SELECT_DATA / "host_taxid_to_accession.json"
HOST_CDS_DIR    = Path("kraken/output/db/search/data/cds/host")

OUT_DIR        = Path("kraken/output/run/split")
DATA_DIR       = OUT_DIR / "data"
INDEX_DIR      = DATA_DIR / "index"
SPLIT_READS_DIR = DATA_DIR / "reads"
SPLIT_RESULTS  = DATA_DIR / "split_results.tsv"

SPLIT_RESULTS_COLS = [
    "Run", "BioSample", "BioProject", "llm_host_resolved",
    "confirmed_host_taxid", "confirmed_host_mapped_pct", "agrees_with_llm",
    "candidate_results", "status",
]

# BBMap needs its whole index (roughly proportional to reference size) resident
# in memory. 48g comfortably covers even the largest single host genome in this
# cohort (~22Gb Pinus radiata) with headroom — sized for ONE genome at a time,
# since indices are per-taxid now, not combined.
BBMAP_XMX = "48g"

_FASTA_EXCLUDE_SUFFIXES = ("_clean.fna", "_combined.fna", ".tagged.fna")


# ── reference lookup ──────────────────────────────────────────────────────────

def find_host_fasta(accession: str) -> Path | None:
    """Locate the already-extracted genomic .fna for an accession.
    kraken_db_search.download_cds() (called by kraken_run_select.py) unzips on
    download, so this is just a lookup, never a fetch."""
    d = HOST_CDS_DIR / accession
    if not d.is_dir():
        return None
    fnas = [f for f in d.glob("**/*.fna") if not f.name.endswith(_FASTA_EXCLUDE_SUFFIXES)]
    return fnas[0] if fnas else None


# ── stage 1: build one index per taxid ────────────────────────────────────────

def index_ready(taxid: str) -> bool:
    """BBMap writes ref/genome/1/summary.txt on a successful index build."""
    return (INDEX_DIR / taxid / "ref" / "genome" / "1" / "summary.txt").exists()


def build_index(taxid: str, accession: str) -> str:
    """Build a BBMap index for one host taxid. Resumable: skips if already
    built. Returns 'ok', 'cached', or 'failed'."""
    if index_ready(taxid):
        return "cached"
    fasta = find_host_fasta(accession)
    if not fasta:
        return "failed"
    idx_dir = INDEX_DIR / taxid
    idx_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["bbmap.sh", f"-Xmx{BBMAP_XMX}",
           f"ref={fasta}", f"path={idx_dir}", "build=1", "overwrite=t"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not index_ready(taxid):
        print(f"  taxid {taxid}: index build FAILED\n{r.stderr[-1500:]}", flush=True)
        return "failed"
    return "ok"


# ── stage 2: split ─────────────────────────────────────────────────────────────

def _count_fastq_reads(path: Path) -> int:
    """Count reads in a (possibly gzipped) FASTQ file — 4 lines per read."""
    if not path.exists():
        return 0
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rt") as f:
        for _ in f:
            n += 1
    return n // 4


def align_against_taxid(taxid: str, run: str, r1: Path, r2: Path | None,
                        work_dir: Path) -> dict:
    """Align one run's reads against one taxid's pre-built index. Returns
    {'mapped_pct': float, 'unmapped_r1': Path, 'unmapped_r2': Path|None,
    'n_input': int, 'n_unmapped': int} or {} on failure."""
    if not index_ready(taxid):
        return {}
    idx_dir = INDEX_DIR / taxid
    out_u1 = work_dir / f"{run}__{taxid}__unmapped_1.fastq.gz"
    cmd = ["bbmap.sh", f"-Xmx{BBMAP_XMX}", f"path={idx_dir}", "build=1",
           "overwrite=t", "statsfile=stderr"]
    if r2 is not None:
        out_u2 = work_dir / f"{run}__{taxid}__unmapped_2.fastq.gz"
        cmd += [f"in={r1}", f"in2={r2}", f"outu={out_u1}", f"outu2={out_u2}"]
    else:
        out_u2 = None
        cmd += [f"in={r1}", f"outu={out_u1}"]

    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        return {}

    n_input = _count_fastq_reads(r1) + (_count_fastq_reads(r2) if r2 else 0)
    n_unmapped = _count_fastq_reads(out_u1) + (_count_fastq_reads(out_u2) if out_u2 else 0)
    mapped_pct = (n_input - n_unmapped) / n_input * 100 if n_input else 0.0
    return {"mapped_pct": mapped_pct, "unmapped_r1": out_u1, "unmapped_r2": out_u2,
            "n_input": n_input, "n_unmapped": n_unmapped}


def split_run(row: dict, work_dir: Path) -> dict:
    """Run stage 2 for one Run: align against every candidate host's index,
    pick the highest-mapping candidate as confirmed host, keep its unmapped
    reads as the final pathogen-enriched output. Returns a SPLIT_RESULTS row."""
    run = row["Run"]
    r1 = SELECT_READS_DIR / f"{run}_1.fastq.gz"
    r2 = SELECT_READS_DIR / f"{run}_2.fastq.gz"
    se = SELECT_READS_DIR / f"{run}.fastq.gz"
    paired = r1.exists() and r2.exists()
    if not paired and not se.exists():
        return {**{k: "" for k in SPLIT_RESULTS_COLS}, "Run": run,
                "BioSample": row.get("BioSample", ""), "status": "no_reads"}

    candidates = sorted({t.strip() for t in (row.get("candidate_host_taxids", "") or "").split(";") if t.strip()})
    results = {}
    for taxid in candidates:
        if paired:
            res = align_against_taxid(taxid, run, r1, r2, work_dir)
        else:
            res = align_against_taxid(taxid, run, se, None, work_dir)
        if res:
            results[taxid] = res

    if not results:
        return {**{k: "" for k in SPLIT_RESULTS_COLS}, "Run": run,
                "BioSample": row.get("BioSample", ""),
                "BioProject": row.get("BioProject", ""),
                "llm_host_resolved": row.get("llm_host_resolved", ""),
                "status": "no_index_available"}

    winner = max(results, key=lambda t: results[t]["mapped_pct"])
    win = results[winner]

    # Promote the winner's unmapped reads to the final output location.
    SPLIT_READS_DIR.mkdir(parents=True, exist_ok=True)
    final_r1 = SPLIT_READS_DIR / f"{run}_1.fastq.gz"
    win["unmapped_r1"].rename(final_r1)
    final_r2 = None
    if win.get("unmapped_r2"):
        final_r2 = SPLIT_READS_DIR / f"{run}_2.fastq.gz"
        win["unmapped_r2"].rename(final_r2)

    # Clean up the non-winning candidates' unmapped FASTQ (kept only stats).
    for t, res in results.items():
        if t == winner:
            continue
        for key in ("unmapped_r1", "unmapped_r2"):
            p = res.get(key)
            if p and p.exists():
                p.unlink()

    candidate_summary = "; ".join(
        f"{t}:{results[t]['mapped_pct']:.1f}%" for t in sorted(results))
    resolved_taxid = (row.get("host_taxid", "") or "").strip()
    agrees = "" if not resolved_taxid else str(resolved_taxid == winner)

    return {
        "Run": run, "BioSample": row.get("BioSample", ""),
        "BioProject": row.get("BioProject", ""),
        "llm_host_resolved": row.get("llm_host_resolved", ""),
        "confirmed_host_taxid": winner,
        "confirmed_host_mapped_pct": f"{win['mapped_pct']:.2f}",
        "agrees_with_llm": agrees,
        "candidate_results": candidate_summary,
        "status": "ok",
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-index", action="store_true",
                    help="Stop after building indices; don't run the split stage")
    ap.add_argument("--limit", type=int, default=None,
                    help="Process at most N runs in the split stage (for testing)")
    ap.add_argument("--workers", type=int, default=4,
                    help="Parallel bbmap.sh jobs (default 4 — kept modest since "
                         "some host genomes are tens of Gb and each job gets its "
                         "own -Xmx allocation; too much concurrency risks "
                         "overcommitting node memory)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_run_split.log")
    link_latest(logs_base, log_dir / "kraken_run_split.log")
    sys.stdout = log

    try:
        if not RUN_LIST.exists():
            sys.exit(f"Error: {RUN_LIST} not found — run kraken_run_select.py --download first.")
        taxid_to_accession = {str(k): v for k, v in load_json(HOST_TAXID_MAP).items()}
        if not taxid_to_accession:
            sys.exit(f"Error: {HOST_TAXID_MAP} not found/empty — run "
                     f"kraken_run_select.py --download first.")

        # ── stage 1: build indices ─────────────────────────────────────────
        print(f"Building indices for {len(taxid_to_accession)} distinct host taxids "
              f"({args.workers} workers) …", flush=True)
        n_ok = n_cached = n_fail = 0
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(build_index, taxid, acc): taxid
                    for taxid, acc in taxid_to_accession.items() if acc}
            for done, fut in enumerate(as_completed(futs), 1):
                taxid = futs[fut]
                status = fut.result()
                if status == "ok":
                    n_ok += 1
                elif status == "cached":
                    n_cached += 1
                else:
                    n_fail += 1
                if done % 10 == 0 or done == len(futs):
                    elapsed = time.time() - t0
                    print(f"  [{done}/{len(futs)}] built={n_ok} cached={n_cached} "
                          f"failed={n_fail}  ({elapsed:.0f}s)", flush=True)
        print(f"Index build complete: {n_ok} built, {n_cached} already cached, "
              f"{n_fail} failed", flush=True)

        if args.build_index:
            print("\n--build-index: stopping before the split stage.")
            return

        # ── stage 2: split ──────────────────────────────────────────────────
        with open(RUN_LIST, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if args.limit:
            rows = rows[:args.limit]
        print(f"\nSplitting {len(rows)} runs against their candidate host "
              f"indices ({args.workers} workers) …", flush=True)

        work_dir = DATA_DIR / "_tmp"
        work_dir.mkdir(parents=True, exist_ok=True)
        results = []
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(split_run, row, work_dir): row["Run"] for row in rows}
            for done, fut in enumerate(as_completed(futs), 1):
                results.append(fut.result())
                if done % 10 == 0 or done == len(rows):
                    elapsed = time.time() - t0
                    n_status = {}
                    for r in results:
                        n_status[r["status"]] = n_status.get(r["status"], 0) + 1
                    print(f"  [{done}/{len(rows)}] {n_status}  ({elapsed:.0f}s)", flush=True)

        with open(SPLIT_RESULTS, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=SPLIT_RESULTS_COLS, delimiter="\t")
            w.writeheader()
            w.writerows(results)
        print(f"\nOutput: {SPLIT_RESULTS}")
        print(f"Next: python kraken/run/kraken_run_assign.py")
    finally:
        log.close()


if __name__ == "__main__":
    main()
