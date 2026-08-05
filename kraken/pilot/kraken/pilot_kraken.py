#!/usr/bin/env python3
"""
Kraken2 MAL gate pilot — 45 runs × 3 organisms × 3 tiers.

Builds a mini Kraken2 database from the 3 transcriptome FASTAs already
present in kraken/pilot/kraken/refs/ (same references used by the kallisto pilot),
enabling a three-way comparison: STAT (euk_pct) vs kallisto vs Kraken2.

Organisms:
  pst — Puccinia striiformis f. sp. tritici  taxid 208830  refs/pst.rna.fa
  pgt — Puccinia graminis f. sp. tritici     taxid 208827  refs/pgt.rna.fa
  por — Pyricularia oryzae 70-15             taxid 318829  refs/por.rna.fa

Taxonomy note: uses kraken:taxid|TAXID| header format to avoid downloading
the 15 GB accession2taxid files. Only taxdump.tar.gz (~57 MB) is needed.

Using RNA refs (same as kallisto) keeps the local build fast and disk-light.
The Setonix production run will use full genomic FASTAs (kraken_build.py).

Run from crypt/ with kraken2 conda env active:
    conda activate kraken2
    python kraken/pilot/kraken/pilot_kraken.py
    python kraken/pilot/kraken/pilot_kraken.py --skip-build   # DB already built
    python kraken/pilot/kraken/pilot_kraken.py --runs 5       # test first 5 runs
"""

import argparse
import csv
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request
from pathlib import Path

PILOT_DIR    = Path("kraken/pilot/kraken")
REFS_DIR     = PILOT_DIR / "refs"
DB_DIR       = PILOT_DIR / "kraken_db"
RESULT_TSV   = PILOT_DIR / "kraken_results.tsv"
KALLISTO_TSV = PILOT_DIR / "results.tsv"

ORGANISMS = {
    "pst": {"ref": REFS_DIR / "pst.rna.fa", "taxid": 208830,
            "name": "Puccinia striiformis f. sp. tritici"},
    "pgt": {"ref": REFS_DIR / "pgt.rna.fa", "taxid": 208827,
            "name": "Puccinia graminis f. sp. tritici"},
    "por": {"ref": REFS_DIR / "por.rna.fa", "taxid": 318829,
            "name": "Pyricularia oryzae"},
}

N_READS = 500_000
CONFIDENCE = 0.1
KRAKEN_THREADS = 8

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


# ── Taxonomy setup ────────────────────────────────────────────────────────────

def download_taxonomy(db_dir: Path) -> None:
    """Download NCBI taxdump (names.dmp + nodes.dmp only, ~57 MB)."""
    tax_dir = db_dir / "taxonomy"
    if (tax_dir / "names.dmp").exists() and (tax_dir / "nodes.dmp").exists():
        print("  taxonomy already present, skipping")
        return

    tax_dir.mkdir(parents=True, exist_ok=True)
    archive = tax_dir / "taxdump.tar.gz"

    print(f"  downloading taxdump (~57 MB) …", flush=True)
    urllib.request.urlretrieve(TAXDUMP_URL, archive)

    print("  extracting names.dmp + nodes.dmp …", flush=True)
    with tarfile.open(archive) as tf:
        for member in tf.getmembers():
            if member.name in ("names.dmp", "nodes.dmp"):
                tf.extract(member, tax_dir)
    archive.unlink()
    print("  taxonomy ready")


# ── Reference tagging ────────────────────────────────────────────────────────

def tag_ref(org_key: str, src: Path, taxid: int) -> Path:
    """
    Copy the RNA FASTA into DB_DIR/tagged/ with kraken:taxid|TAXID| headers
    so Kraken2 bypasses accession2taxid lookup entirely.
    Returns path to tagged file; skips if already done.
    """
    tagged_dir = DB_DIR / "tagged"
    tagged_dir.mkdir(parents=True, exist_ok=True)
    out = tagged_dir / f"{org_key}.fa"

    if out.exists():
        print(f"  [{org_key}] already tagged, skipping")
        return out

    print(f"  [{org_key}] tagging {src.name} with taxid {taxid} …", flush=True)
    with open(src) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                rest = line[1:].rstrip()
                fout.write(f">kraken:taxid|{taxid}|{rest}\n")
            else:
                fout.write(line)
    print(f"    → {out}")
    return out


# ── Kraken2 build ─────────────────────────────────────────────────────────────

def build_database(threads: int) -> None:
    """Download taxonomy, tag refs, add to library, build mini Kraken2 DB."""
    DB_DIR.mkdir(parents=True, exist_ok=True)

    print("\n[1/3] Downloading taxonomy …")
    download_taxonomy(DB_DIR)

    print("\n[2/3] Tagging RNA refs + adding to library …")
    for org_key, org in ORGANISMS.items():
        if not org["ref"].exists():
            print(f"  ERROR: {org['ref']} not found — run pilot_kallisto.py first")
            sys.exit(1)
        tagged = tag_ref(org_key, org["ref"], org["taxid"])
        r = subprocess.run(
            ["kraken2-build", "--add-to-library", str(tagged), "--db", str(DB_DIR)],
            capture_output=True, text=True
        )
        if r.returncode != 0:
            print(f"  WARNING: add-to-library failed for {org_key}: {r.stderr[:120]}")

    print("\n[3/3] Building Kraken2 database …")
    r = subprocess.run(
        ["kraken2-build", "--build", "--db", str(DB_DIR),
         "--threads", str(threads)],
        text=True
    )
    if r.returncode != 0:
        raise RuntimeError("kraken2-build --build failed")
    print("  database ready")


# ── ENA fetch + classify ──────────────────────────────────────────────────────

def _ena_urls(run: str) -> list[str]:
    api = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={run}&result=read_run&fields=fastq_ftp,library_layout&format=tsv")
    try:
        req = urllib.request.Request(api, headers={"User-Agent": "crypt/pilot_kraken"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = resp.read().decode()
    except Exception as e:
        print(f"    ENA API error: {e}")
        return []
    lines = [l for l in body.strip().split("\n") if l]
    if len(lines) < 2:
        return []
    cols = lines[0].split("\t")
    vals = lines[1].split("\t")
    ftp_field = vals[cols.index("fastq_ftp")] if "fastq_ftp" in cols else vals[-1]
    return [f"ftp://{p.strip()}" for p in ftp_field.split(";") if p.strip()]


def _stream_reads(url: str, n: int, dest: Path) -> bool:
    cmd = ["bash", "-c",
           f'curl --silent --fail --max-time 300 "{url}" | gunzip -c | head -n {n * 4}']
    try:
        with open(dest, "w") as out:
            subprocess.run(cmd, stdout=out, timeout=360)
        return dest.exists() and dest.stat().st_size > 0
    except Exception:
        return False


def classify_run(run: str, layout: str, db_dir: Path,
                 tmp_dir: Path) -> dict:
    """Stream + classify one run. Returns {org_key: pct, pct_classified, n_reads, error}."""
    result = {"run": run, "error": None}

    urls = _ena_urls(run)
    if not urls:
        result["error"] = "no_ena_urls"
        return result

    run_tmp = tmp_dir / run
    run_tmp.mkdir(parents=True, exist_ok=True)

    try:
        is_paired = layout == "PAIRED" and len(urls) >= 2
        if is_paired:
            r1, r2 = run_tmp / "r1.fastq", run_tmp / "r2.fastq"
            if not (_stream_reads(urls[0], N_READS, r1) and
                    _stream_reads(urls[1], N_READS, r2)):
                result["error"] = "stream_failed"
                return result
            reads = [r1, r2]
        else:
            r1 = run_tmp / "r1.fastq"
            if not _stream_reads(urls[0], N_READS, r1):
                result["error"] = "stream_failed"
                return result
            reads = [r1]

        report = run_tmp / "report.txt"
        cmd = ["kraken2", "--db", str(db_dir),
               "--report", str(report),
               "--confidence", str(CONFIDENCE),
               "--threads", str(KRAKEN_THREADS),
               "--output", "/dev/null"]
        if is_paired:
            cmd += ["--paired", str(reads[0]), str(reads[1])]
        else:
            cmd.append(str(reads[0]))

        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if r.returncode != 0 or not report.exists():
            result["error"] = "kraken2_failed"
            return result

        # Parse report: find each organism by taxid
        taxid_pct: dict[int, float] = {}
        pct_classified = 0.0
        n_reads = 0
        with open(report) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                pct   = float(parts[0])
                reads_count = int(parts[1])
                rank  = parts[3].strip()
                taxid = int(parts[4].strip())

                if rank == "U":
                    # will compute classified from R line
                    pass
                elif rank == "R" and taxid == 1:
                    pct_classified = pct
                    n_reads_classified = reads_count

                if taxid in {o["taxid"] for o in ORGANISMS.values()}:
                    taxid_pct[taxid] = round(pct, 4)

        result["pct_classified"] = round(pct_classified, 4)
        for org_key, org in ORGANISMS.items():
            result[f"{org_key}_pct"] = taxid_pct.get(org["taxid"], 0.0)

    finally:
        for f in run_tmp.iterdir():
            f.unlink(missing_ok=True)
        run_tmp.rmdir()

    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-build", action="store_true",
                    help="Skip DB build (use existing kraken_db/)")
    ap.add_argument("--runs", type=int, default=None,
                    help="Process only first N runs (default: all 45)")
    ap.add_argument("--threads", type=int, default=KRAKEN_THREADS,
                    help=f"Kraken2 threads (default: {KRAKEN_THREADS})")
    args = ap.parse_args()

    # ── Build ──────────────────────────────────────────────────────────────────
    if not args.skip_build:
        build_database(args.threads)
    else:
        hash_file = DB_DIR / "hash.k2d"
        if not hash_file.exists():
            print(f"ERROR: --skip-build set but {hash_file} not found. Run without flag first.")
            sys.exit(1)
        print(f"Skipping build, using existing DB: {DB_DIR}")

    # ── Load pilot runs ────────────────────────────────────────────────────────
    if not KALLISTO_TSV.exists():
        print(f"ERROR: {KALLISTO_TSV} not found — run pilot_kallisto.py first")
        sys.exit(1)

    with open(KALLISTO_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        runs = list(reader)

    if args.runs:
        runs = runs[:args.runs]

    print(f"\nClassifying {len(runs)} pilot runs …\n")

    tmp_dir = Path(tempfile.gettempdir()) / "kraken_pilot"
    tmp_dir.mkdir(exist_ok=True)

    results = []
    t0 = time.time()

    for i, row in enumerate(runs, 1):
        run_id  = row["Run"]
        org_key = row["org_key"]
        layout  = "PAIRED"   # let URL count decide; classify_run treats 2-URL as PE

        print(f"[{i:>2}/{len(runs)}] {run_id} ({org_key}/{row['tier']}) …", flush=True)
        res = classify_run(run_id, layout, DB_DIR, tmp_dir)

        if res.get("error"):
            print(f"  ERROR: {res['error']}")
        else:
            target_pct = res.get(f"{org_key}_pct", 0.0)
            classified = res.get("pct_classified", 0.0)
            print(f"  classified={classified:.1f}%  {org_key}_pct={target_pct:.2f}%  "
                  f"kallisto={row['kallisto_pct']}%  STAT_euk={row['euk_pct']}%")

        results.append({**row, **res})

    # ── Write TSV ──────────────────────────────────────────────────────────────
    fieldnames = list(runs[0].keys()) + [
        "pst_pct", "pgt_pct", "por_pct", "pct_classified", "error"
    ]
    with open(RESULT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(results)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed/60:.1f} min — results: {RESULT_TSV}")

    # ── Quick summary ──────────────────────────────────────────────────────────
    print("\nTier × organism summary (mean target_pct):")
    from collections import defaultdict
    by_tier = defaultdict(list)
    for r in results:
        ok = r.get("org_key", "")
        pct = float(r.get(f"{ok}_pct", 0))
        by_tier[(r.get("org_key"), r.get("tier"))].append(pct)
    for (org, tier), vals in sorted(by_tier.items()):
        mean = sum(vals) / len(vals) if vals else 0
        print(f"  {org:>3}  {tier:<5}  mean kraken_pct = {mean:.2f}%  (n={len(vals)})")


if __name__ == "__main__":
    main()
