#!/usr/bin/env python3
"""
Kallisto MAL gate pilot — 45 runs × 3 organisms × 3 tiers.

Tiers:
  zero  — euk_pct = 0 (zero-kmer, negative control)
  low   — 0 < euk_pct < 1% (Dikarya-stall, main test group)
  high  — euk_pct >= 1%, library_detected=True (positive controls)

Reference assemblies (NCBI RefSeq RNA FASTA):
  pst — Puccinia striiformis f. sp. tritici  GCF_021901695.1
  pgt — Puccinia graminis f. sp. tritici     GCF_000149925.1
  por — Pyricularia oryzae 70-15             GCF_000002495.2

Run from crypt/: python scripts/pilot_kallisto.py
  --skip-download   if refs already downloaded
  --skip-index      if indices already built
  --runs N          process only first N runs (default: all)
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

PILOT_DIR  = Path("scripts/pilot")
REFS_DIR   = PILOT_DIR / "refs"
IDX_DIR    = PILOT_DIR / "indices"
READS_DIR  = PILOT_DIR / "reads"
QUANT_DIR  = PILOT_DIR / "quant"
RESULT_TSV = PILOT_DIR / "results.tsv"
RUNS_TSV   = Path("scripts/pilot_runs.tsv")

ASSEMBLIES = {
    "pst": "GCF_021901695.1",   # P. striiformis f. sp. tritici Pst134E36
    "pgt": "GCF_000149925.1",   # P. graminis f. sp. tritici CRL 75-36-700-3
    "por": "GCF_000002495.2",   # Pyricularia oryzae 70-15
}

N_READS = 100_000   # spots to stream per run


def run(cmd, **kwargs):
    print(f"  $ {' '.join(str(c) for c in cmd)}", flush=True)
    r = subprocess.run(cmd, **kwargs)
    if r.returncode != 0:
        print(f"  ERROR: exit code {r.returncode}", file=sys.stderr)
    return r.returncode == 0


def download_refs(skip=False):
    REFS_DIR.mkdir(parents=True, exist_ok=True)
    for org, acc in ASSEMBLIES.items():
        out_fa = REFS_DIR / f"{org}.rna.fa"
        if skip and out_fa.exists():
            print(f"  [skip] {out_fa} exists")
            continue
        zip_path = REFS_DIR / f"{org}.zip"
        print(f"\n── Downloading {org} ({acc}) ──")
        ok = run(["datasets", "download", "genome", "accession", acc,
                  "--include", "rna", "--filename", str(zip_path)])
        if not ok:
            print(f"  FAILED to download {acc}", file=sys.stderr)
            continue
        # extract rna.fna from zip
        with zipfile.ZipFile(zip_path) as zf:
            rna_names = [n for n in zf.namelist() if n.endswith("rna.fna")]
            if not rna_names:
                print(f"  ERROR: no rna.fna in {zip_path}", file=sys.stderr)
                continue
            with zf.open(rna_names[0]) as src, open(out_fa, "wb") as dst:
                shutil.copyfileobj(src, dst)
        zip_path.unlink()
        print(f"  → {out_fa}  ({out_fa.stat().st_size // 1_000_000} MB)")


def build_indices(skip=False):
    IDX_DIR.mkdir(parents=True, exist_ok=True)
    for org in ASSEMBLIES:
        idx = IDX_DIR / f"{org}.idx"
        if skip and idx.exists():
            print(f"  [skip] {idx} exists")
            continue
        fa = REFS_DIR / f"{org}.rna.fa"
        if not fa.exists():
            print(f"  MISSING {fa} — run without --skip-download first", file=sys.stderr)
            continue
        print(f"\n── Building kallisto index: {org} ──")
        run(["kallisto", "index", "-i", str(idx), str(fa)])


def _ena_fastq_urls(acc):
    """Return list of FASTQ FTP paths from ENA portal API (empty if not found)."""
    import urllib.request
    api = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={acc}&result=read_run&fields=fastq_ftp")
    try:
        with urllib.request.urlopen(api, timeout=15) as r:
            lines = r.read().decode().strip().split("\n")
        if len(lines) < 2:
            return []
        urls = [u.strip() for u in lines[1].split("\t")[1].split(";") if u.strip()]
        return urls  # like ["ftp.sra.ebi.ac.uk/vol1/fastq/ERR.../ERR..._1.fastq.gz", ...]
    except Exception:
        return []


def _stream_ena(ftp_path, out_fastq, n_reads):
    """Stream n_reads from an ENA FTP gzipped FASTQ path via curl | gunzip | head."""
    n_lines = n_reads * 4
    url = f"ftp://{ftp_path}"
    # head sends SIGPIPE when limit reached → curl exits 23 (write error); treat as ok
    cmd = f'curl --silent --fail "{url}" | gunzip -c | head -n {n_lines} > {out_fastq}'
    r = subprocess.run(cmd, shell=True)
    return r.returncode in (0, 23, 141) and Path(out_fastq).exists()


def stream_and_quant(run_row):
    """Fetch N_READS reads, run kallisto quant --single, return p_pseudoaligned."""
    acc      = run_row["Run"]
    org_key  = run_row["org_key"]
    idx      = IDX_DIR / f"{org_key}.idx"
    out_dir  = QUANT_DIR / acc
    reads_se = READS_DIR / f"{acc}.fastq"

    READS_DIR.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not idx.exists():
        print(f"  MISSING index {idx}", file=sys.stderr)
        return None

    # ── Skip if already done ──────────────────────────────────────────────────
    run_info = out_dir / "run_info.json"
    if run_info.exists():
        with open(run_info) as f:
            info = json.load(f)
        p = info.get("p_pseudoaligned", 0.0)
        print(f"  [cached] → {p:.1f}%")
        return p

    # ── Fetch reads via ENA FTP streaming (fast partial download) ─────────────
    print(f"  fetching {acc} (first {N_READS:,} reads) ...", flush=True)
    ena_urls = _ena_fastq_urls(acc)
    if not ena_urls:
        print(f"  no ENA URLs for {acc} — skipping", file=sys.stderr)
        return None

    # use _1 file only (forward reads); enough for detection
    r1_urls = [u for u in ena_urls if "_1.fastq" in u or (len(ena_urls) == 1)]
    target_url = r1_urls[0] if r1_urls else ena_urls[0]
    print(f"  streaming {target_url.split('/')[-1]} ...", flush=True)

    ok = _stream_ena(target_url, reads_se, N_READS)
    if not ok or not reads_se.exists() or reads_se.stat().st_size == 0:
        print(f"  ENA stream failed for {acc}", file=sys.stderr)
        reads_se.unlink(missing_ok=True)
        return None

    read_size_mb = reads_se.stat().st_size // 1_000_000
    print(f"  {read_size_mb} MB downloaded", flush=True)

    # ── Kallisto quant ────────────────────────────────────────────────────────
    print(f"  kallisto quant ({read_size_mb} MB) ...", flush=True)
    ok = run(["kallisto", "quant",
              "-i", str(idx),
              "-o", str(out_dir),
              "--single", "-l", "200", "-s", "20",
              str(reads_se)],
             capture_output=True)

    reads_se.unlink(missing_ok=True)   # free disk immediately

    if not ok:
        return None

    # ── Parse result ──────────────────────────────────────────────────────────
    run_info = out_dir / "run_info.json"
    if not run_info.exists():
        print(f"  Missing run_info.json for {acc}", file=sys.stderr)
        return None
    with open(run_info) as f:
        info = json.load(f)
    return info.get("p_pseudoaligned", 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-download", action="store_true")
    ap.add_argument("--skip-index",    action="store_true")
    ap.add_argument("--runs",          type=int, default=None,
                    help="Process only first N runs")
    args = ap.parse_args()

    PILOT_DIR.mkdir(parents=True, exist_ok=True)

    download_refs(skip=args.skip_download)
    build_indices(skip=args.skip_index)

    # ── Load run list ─────────────────────────────────────────────────────────
    with open(RUNS_TSV) as f:
        run_rows = list(csv.DictReader(f, delimiter="\t"))
    if args.runs:
        run_rows = run_rows[:args.runs]

    # ── Process runs ──────────────────────────────────────────────────────────
    QUANT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    for i, row in enumerate(run_rows, 1):
        print(f"\n[{i}/{len(run_rows)}] {row['Run']}  org={row['org_key']}  "
              f"tier={row['tier']}  stat_euk_pct={row['euk_pct']}")
        p = stream_and_quant(row)
        result = {**row, "kallisto_pct": f"{p:.3f}" if p is not None else "NA"}
        results.append(result)
        print(f"  → kallisto pseudoalignment: "
              f"{'NA' if p is None else f'{p:.1f}%'}")

    # ── Write results ─────────────────────────────────────────────────────────
    fields = list(run_rows[0].keys()) + ["kallisto_pct"]
    with open(RESULT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    print(f"\n── Summary ──────────────────────────────────────────────────────")
    print(f"{'Organism':<40} {'Tier':<6} {'STAT euk%':>10} {'Kallisto%':>10}")
    print("-" * 70)
    for row in results:
        print(f"{row['organism']:<40} {row['tier']:<6} "
              f"{row['euk_pct']:>10}  {row['kallisto_pct']:>10}")

    print(f"\nResults written: {RESULT_TSV}")


if __name__ == "__main__":
    main()
