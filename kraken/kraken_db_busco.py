#!/usr/bin/env python3
"""
kraken_db_busco.py — BUSCO-score every candidate assembly from
kraken_db_search.py and apply completeness thresholds. Submodule 1, step 2 of 3
(search → busco → build).

Reads CDS FASTA that kraken_db_search.py already downloaded to
kraken_db_search/data/cds_v2/{accession}/ — this script never downloads anything
itself. Runs BUSCO in transcriptome mode per accession, then applies per-kingdom
completeness thresholds with a per-taxid fallback (if nothing passes for a taxid,
keep the best-scoring option anyway, flagged, so kraken_db_build.py always has
something to build with).

Output is ONE merged table (candidate metadata + BUSCO score + pass/fail +
final selection decision) rather than the old two-file busco_scores.tsv /
ref_screen.tsv split — `selected=True` rows are exactly what
kraken_db_build.py should add to the Kraken2 library.

Run on Setonix: sbatch kraken/slurm/kraken_db_busco.slurm
Run from crypt/:
    python kraken/kraken_db_busco.py \\
        --busco-db-path /scratch/pawsey1168/llenzo/crypt/kraken/output/kraken_db_busco/data/busco_downloads \\
        [--genomes-dir /scratch/pawsey1168/llenzo/crypt/kraken/output/kraken_db_search/data/cds_v2] \\
        [--busco-out   /scratch/pawsey1168/llenzo/kraken/busco] \\
        [--workers 16] [--cpus-per-busco 8]
        [--fungi-threshold 50.0] [--oomycete-threshold 65.0]

Outputs:
    kraken/output/kraken_db_busco/data/busco_scores.tsv        (tracked — final merged table)
    kraken/output/kraken_db_busco/data/busco_scan_cache.tsv    (tracked — resumable raw score cache)
    {busco-out}/{accession}/                                    (BUSCO run dirs, deleted after scoring)
"""

import argparse
import csv
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _util import _Tee, make_log_dir, link_latest

CANDIDATES_TSV = Path("kraken/output/kraken_db_search/data/ref_candidates.tsv")
DEFAULT_GENOMES_DIR = Path("kraken/output/kraken_db_search/data/cds_v2")
DEFAULT_BUSCO_DB    = Path("kraken/output/kraken_db_busco/data/busco_downloads")
DEFAULT_BUSCO_OUT   = Path("/scratch/pawsey1168/llenzo/kraken/busco")  # ephemeral, cleared per-run

OUT_DIR    = Path("kraken/output/kraken_db_busco")
DATA_DIR   = OUT_DIR / "data"
SCAN_CACHE = DATA_DIR / "busco_scan_cache.tsv"
FINAL_TSV  = DATA_DIR / "busco_scores.tsv"

SCAN_COLS = [
    "accession", "organism_name", "kingdom", "busco_lineage",
    "complete_pct", "single_pct", "duplicated_pct",
    "fragmented_pct", "missing_pct", "n_markers",
    "status",   # pass | fail | no_cds | busco_error
]

FINAL_COLS = [
    "taxid", "organism_name", "kingdom", "source", "accession",
    "assembly_level", "release_date", "has_annotation", "country",
    "protein_coding_genes", "scaffold_n50_kb", "total_length_mb",
    "fasta_type", "busco_lineage", "diversity_rank", "diversity_reason",
    "complete_pct", "single_pct", "duplicated_pct", "fragmented_pct",
    "missing_pct", "n_markers", "busco_status",
    "selected", "selected_reason",
]

# Regex for the BUSCO summary line: C:91.2%[S:88.4%,D:2.8%],F:3.1%,M:5.7%,n:758
_BUSCO_RE = re.compile(
    r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],"
    r"F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)"
)


def _ts() -> str:
    return time.strftime("%H:%M:%S")


def _mb(path: Path) -> str:
    try:
        return f"{path.stat().st_size / 1e6:.1f} MB"
    except OSError:
        return "? MB"


# ── FASTA prep (concat + clean; CDS itself already downloaded by kraken_db_search.py) ─

def concat_fnas(fnas: list, dest: Path, accession: str) -> Path:
    if len(fnas) == 1:
        return fnas[0]
    print(f"  [{_ts()}] {accession}  concat: merging {len(fnas)} files → {dest.name}", flush=True)
    with open(dest, "wb") as fout:
        for fna in sorted(fnas):
            with open(fna, "rb") as fin:
                shutil.copyfileobj(fin, fout)
    return dest


def clean_fasta_for_busco(fna: Path, dest: Path, accession: str) -> Path:
    """Strip kraken:taxid|N| prefix (added by kraken_db_build.py on a prior run,
    if this DB is being rebuilt in place); return original if already clean."""
    if dest.exists():
        return dest
    with open(fna, "rb") as f:
        head = f.read(30)
    if b">kraken:taxid|" not in head:
        return fna
    print(f"  [{_ts()}] {accession}  clean: stripping kraken:taxid headers...", flush=True)
    content = fna.read_bytes()
    cleaned = re.sub(rb">kraken:taxid\|\d+\|", b">", content)
    dest.write_bytes(cleaned)
    return dest


# ── BUSCO run + parse ─────────────────────────────────────────────────────────

def run_busco(fna: Path, lineage: str, out_name: str,
              out_path: Path, db_path: Path, cpus: int):
    accession = out_name
    run_dir   = out_path / out_name

    existing = list(run_dir.glob(f"short_summary.specific.{lineage}.*.txt"))
    if not existing:
        existing = list(run_dir.glob(f"run_{lineage}/short_summary.txt"))

    if existing:
        print(f"  [{_ts()}] {accession}  busco: summary already exists"
              f" ({existing[0].name}), skipping run", flush=True)
    else:
        run_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "busco", "-m", "transcriptome", "-i", str(fna), "-l", lineage,
            "--out", out_name, "--out_path", str(out_path), "--cpu", str(cpus),
            "--offline", "--download_path", str(db_path), "--force",
        ]
        print(f"  [{_ts()}] {accession}  busco: starting {lineage} on {fna.name} ({_mb(fna)})...",
              flush=True)
        t0 = time.time()
        # stdout/stderr → DEVNULL + proc.wait() (no pipe threads) avoids the deadlock
        # where MetaEuk/Diamond children keep captured pipes open after the parent is
        # killed. start_new_session + killpg on timeout kills the whole process tree.
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                start_new_session=True)
        try:
            proc.wait(timeout=3600)
        except subprocess.TimeoutExpired:
            elapsed = time.time() - t0
            print(f"  [{_ts()}] {accession}  busco: TIMEOUT after {elapsed:.0f}s — killing process group",
                  flush=True)
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
            return None

        elapsed = time.time() - t0
        if proc.returncode != 0:
            print(f"  [{_ts()}] {accession}  busco: FAILED rc={proc.returncode} after {elapsed:.0f}s",
                  flush=True)
            return None
        print(f"  [{_ts()}] {accession}  busco: completed in {elapsed:.0f}s (rc=0)", flush=True)

        existing = list(run_dir.glob(f"short_summary.specific.{lineage}.*.txt"))
        if not existing:
            existing = list(run_dir.glob(f"run_{lineage}/short_summary.txt"))

    if not existing:
        print(f"  [{_ts()}] {accession}  busco: no summary found after run", flush=True)
        return None
    return parse_busco_summary(existing[0])


def parse_busco_summary(summary_path: Path):
    m = _BUSCO_RE.search(summary_path.read_text())
    if not m:
        return None
    return {
        "complete_pct": float(m.group(1)), "single_pct": float(m.group(2)),
        "duplicated_pct": float(m.group(3)), "fragmented_pct": float(m.group(4)),
        "missing_pct": float(m.group(5)), "n_markers": int(m.group(6)),
    }


# ── per-accession worker ──────────────────────────────────────────────────────

def process_one(row: dict, genomes_dir: Path, busco_out: Path,
                busco_db: Path, cpus: int, thresholds: dict) -> dict:
    acc, name, kingdom, lineage, ftype = (
        row["accession"], row["organism_name"], row["kingdom"],
        row["busco_lineage"], row["fasta_type"],
    )
    result = {"accession": acc, "organism_name": name, "kingdom": kingdom,
              "busco_lineage": lineage, "complete_pct": "", "single_pct": "",
              "duplicated_pct": "", "fragmented_pct": "", "missing_pct": "",
              "n_markers": "", "status": ""}

    if ftype != "cds":
        result["status"] = "no_cds"
        return result

    dest_dir = genomes_dir / acc
    fnas = [f for f in dest_dir.glob("**/*.fna")
            if not f.name.endswith(("_clean.fna", "_combined.fna", ".tagged.fna"))]
    if not fnas:
        print(f"  [{_ts()}] {acc}  no CDS on disk — run kraken_db_search.py --download first",
              flush=True)
        result["status"] = "no_cds"
        return result

    combined  = dest_dir / f"{acc}_combined.fna"
    input_fna = concat_fnas(fnas, combined, acc)
    input_fna = clean_fasta_for_busco(input_fna, dest_dir / f"{acc}_clean.fna", acc)

    scores = run_busco(input_fna, lineage, acc, busco_out, busco_db, cpus)
    if scores is None:
        result["status"] = "busco_error"
        return result

    result.update(scores)
    threshold = thresholds.get(kingdom, 50.0)
    result["status"] = "pass" if scores["complete_pct"] >= threshold else "fail"
    return result


# ── scan (resumable, writes SCAN_CACHE incrementally) ────────────────────────

def run_scan(candidates: list, genomes_dir: Path, busco_out: Path, busco_db: Path,
            cpus: int, workers: int, thresholds: dict) -> None:
    busco_out.mkdir(parents=True, exist_ok=True)

    done_accs = set()
    if SCAN_CACHE.exists():
        with open(SCAN_CACHE, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                done_accs.add(row["accession"])
    print(f"Candidates: {len(candidates):,}  already scored: {len(done_accs):,}", flush=True)

    todo = [r for r in candidates if r["accession"] not in done_accs]
    if not todo:
        print("All accessions already scored.", flush=True)
        return

    write_header = not SCAN_CACHE.exists() or SCAN_CACHE.stat().st_size == 0
    fh = open(SCAN_CACHE, "a", newline="")
    writer = csv.DictWriter(fh, fieldnames=SCAN_COLS, delimiter="\t", extrasaction="ignore")
    if write_header:
        writer.writeheader()

    completed, n_exception, no_cds, busco_errors = 0, 0, 0, 0
    t_wall = time.time()

    def work(row):
        return process_one(row, genomes_dir, busco_out, busco_db, cpus, thresholds)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(work, row): row["accession"] for row in todo}
        for fut in as_completed(futs):
            acc = futs[fut]
            try:
                result = fut.result()
            except Exception as e:
                print(f"  [{_ts()}] ERROR {acc}: {e}", flush=True)
                n_exception += 1
                continue

            writer.writerow(result)
            fh.flush()

            # Delete BUSCO run dir immediately after scoring — Lustre inode limit
            # (~1M files/user on Setonix); each run writes ~5000 intermediate files
            # and we only need the short_summary, already parsed above.
            run_dir = busco_out / acc
            if run_dir.exists():
                shutil.rmtree(run_dir, ignore_errors=True)

            status = result["status"]
            pct = result.get("complete_pct", "")
            pct_s = f"{pct:.1f}%" if isinstance(pct, float) else "—"
            completed += 1
            if status == "no_cds":
                no_cds += 1
            elif status == "busco_error":
                busco_errors += 1

            total_done = completed + n_exception
            elapsed = time.time() - t_wall
            rate = total_done / elapsed if elapsed > 0 else 0
            eta = (len(todo) - total_done) / rate if rate > 0 else 0
            print(f"[{_ts()}] [{total_done}/{len(todo)}] {acc}  {status}  C={pct_s}"
                  f"  | rate={rate:.2f}/s  ETA={eta/3600:.1f}h", flush=True)

    fh.close()
    print(f"\n── Scan summary ─────────────────────────────────────────────────")
    print(f"  Processed:   {completed:,}")
    print(f"  No CDS:      {no_cds:,}")
    print(f"  BUSCO error: {busco_errors:,}")
    print(f"  Exceptions:  {n_exception:,}")


# ── finalize: merge candidates + scan cache, apply threshold + fallback ──────

def run_finalize(candidates: list, thresholds: dict) -> None:
    scores = {}
    if SCAN_CACHE.exists():
        with open(SCAN_CACHE, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                scores[row["accession"]] = row  # keep last (append-only cache)

    by_taxid = defaultdict(list)
    for cand in candidates:
        acc = cand["accession"]
        kingdom = cand["kingdom"]
        sc = scores.get(acc, {})
        status = sc.get("status", "not_run")
        try:
            pct = float(sc["complete_pct"]) if sc.get("complete_pct") else None
        except ValueError:
            pct = None
        threshold = thresholds.get(kingdom, 50.0)

        merged = dict(cand)
        merged["diversity_rank"]   = cand.get("selection_rank", "")
        merged["diversity_reason"] = cand.get("selection_reason", "")
        for col in ("complete_pct", "single_pct", "duplicated_pct", "fragmented_pct",
                    "missing_pct", "n_markers"):
            merged[col] = sc.get(col, "")
        merged["busco_status"] = (
            "pass" if status == "pass" and pct is not None and pct >= threshold
            else status if status in ("no_cds", "busco_error", "not_run")
            else "fail"
        )
        by_taxid[cand["taxid"]].append(merged)

    final_rows = []
    n_pass, n_fallback, n_no_option = 0, 0, 0
    for taxid, cands in by_taxid.items():
        passing = [c for c in cands if c["busco_status"] == "pass"]
        if passing:
            passing_ids = {id(c) for c in passing}
            for c in cands:
                c["selected"] = id(c) in passing_ids
                c["selected_reason"] = "busco_pass" if id(c) in passing_ids else ""
            final_rows.extend(cands)
            n_pass += len(passing)
        else:
            def fallback_key(c):
                try:
                    pct = float(c["complete_pct"]) if c["complete_pct"] else -1
                except ValueError:
                    pct = -1
                return (1 if c["fasta_type"] == "cds" else 0, pct)
            best = sorted(cands, key=fallback_key, reverse=True)[0]
            for c in cands:
                c["selected"] = (c is best)
            best["selected_reason"] = (
                f"fallback_below_threshold:{best['complete_pct']}%" if best["complete_pct"]
                else "fallback_no_score"
            )
            for c in cands:
                if c is not best:
                    c["selected_reason"] = ""
            final_rows.extend(cands)
            if best["complete_pct"]:
                n_fallback += 1
            else:
                n_no_option += 1

    final_rows.sort(key=lambda r: (r["kingdom"], r["source"], r["organism_name"].lower(),
                                   str(r["diversity_rank"])))

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(FINAL_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FINAL_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(final_rows)

    n_selected = sum(1 for r in final_rows if r["selected"])
    print(f"\n── Finalize summary ─────────────────────────────────────────────")
    print(f"  Total candidates:       {len(final_rows):>6,}")
    print(f"  Selected for build:     {n_selected:>6,}")
    print(f"    BUSCO pass:           {n_pass:>6,}")
    print(f"    Fallback (scored):    {n_fallback:>6,}  (below threshold, only option per taxid)")
    print(f"    Fallback (no score):  {n_no_option:>6,}  (busco_error / no_cds, only option)")
    print(f"  Thresholds: fungal >= {thresholds.get('fungal')}%,  oomycete >= {thresholds.get('oomycete')}%")
    print(f"\nOutput: {FINAL_TSV}")
    print(f"Next:   python kraken/kraken_db_build.py")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes-dir",    default=str(DEFAULT_GENOMES_DIR),
                    help="CDS directory populated by kraken_db_search.py --download")
    ap.add_argument("--busco-out",      default=str(DEFAULT_BUSCO_OUT),
                    help="Ephemeral directory for BUSCO run output (cleared per accession)")
    ap.add_argument("--busco-db-path",  default=str(DEFAULT_BUSCO_DB),
                    help="Path to pre-downloaded BUSCO lineage databases")
    ap.add_argument("--workers",        type=int, default=16)
    ap.add_argument("--cpus-per-busco", type=int, default=8)
    ap.add_argument("--fungi-threshold",    type=float, default=50.0)
    ap.add_argument("--oomycete-threshold", type=float, default=65.0)
    ap.add_argument("--finalize-only", action="store_true",
                    help="Skip scanning; just re-merge candidates + existing scan cache "
                         "(e.g. after changing thresholds)")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logs_base = OUT_DIR / "logs"
    log_dir   = make_log_dir(logs_base)
    log = _Tee(log_dir / "kraken_db_busco.log")
    link_latest(logs_base, log_dir / "kraken_db_busco.log")
    sys.stdout = log

    try:
        thresholds = {"fungal": args.fungi_threshold, "oomycete": args.oomycete_threshold}

        if not CANDIDATES_TSV.exists():
            print(f"ERROR: {CANDIDATES_TSV} not found — run kraken_db_search.py first")
            raise SystemExit(1)
        with open(CANDIDATES_TSV, newline="") as fh:
            candidates = list(csv.DictReader(fh, delimiter="\t"))

        if not args.finalize_only:
            run_scan(candidates, Path(args.genomes_dir), Path(args.busco_out),
                     Path(args.busco_db_path), args.cpus_per_busco, args.workers, thresholds)

        run_finalize(candidates, thresholds)
    finally:
        log.close()


if __name__ == "__main__":
    main()
