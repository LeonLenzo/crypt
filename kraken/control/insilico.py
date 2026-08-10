#!/usr/bin/env python3
"""
insilico.py — synthetic mixed-infection controls for Kraken2 DB validation.

Mixes reads from a known-composition source (SRA run or pre-generated FASTQ)
into a plant host background at specified ratios, runs Kraken2, and records
expected vs observed species proportions.  Directly calibrates KINGDOM_THRESHOLDS
and validates same-genus specificity.

Panel JSON (see --write-template):
  {
    "background": {"run": "SRR...", "species": "Triticum aestivum"},
    "mixtures": [
      {"name": "Botrytis cinerea",     "kingdom": "Fungi",    "run": "SRR..."},
      {"name": "Puccinia striiformis", "kingdom": "Fungi",
       "fastq": ["pst_sim1.fq.gz", "pst_sim2.fq.gz"]}
    ],
    "ratios":     [0.1, 0.5, 1.0, 5.0, 10.0, 25.0],
    "n_total":    500000,
    "replicates": 1
  }

For obligate biotrophs without culture RNA-seq (e.g. PST), simulate reads from the
reference genome first:
  art_illumina -ss HS25 -i pst_cyr32.fa -l 150 -f 15 -m 200 -s 10 -p -o pst_sim
  gzip pst_sim1.fq pst_sim2.fq

Note: this script uses R1 only for mixing.  Paired reads within each source are
used for acquisition (parallel streaming), but classification uses R1 as single-end
to keep mixing arithmetic simple.  This is conservative; PE classification would
give equal or higher observed percentages.

Outputs (under --out-dir):
  results.tsv        — per mixture row: expected_pct, observed_pct, detected, ...
  control_rows.tsv   — stratum-J rows to append to control_runs.tsv

Usage:
  python kraken/control/insilico.py --write-template
  # fill in panel.json, then:
  python kraken/control/insilico.py --panel panel.json --db /scratch/.../db_pathogens
  # append to control_runs.tsv when satisfied:
  # tail -n +2 kraken/control/output/data/insilico/control_rows.tsv >> kraken/control/output/data/control_runs.tsv
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from math import ceil
from pathlib import Path

ROOT        = Path(__file__).resolve().parents[2]
OUT_DIR_DEF = ROOT / "kraken/control/output/data/insilico"
CONTROL_TSV = ROOT / "kraken/control/output/data/control_runs.tsv"
SCRATCH_DEF = Path(os.environ.get("MYSCRATCH", tempfile.gettempdir())) / "insilico_tmp"

CONFIDENCE     = 0.15
MIN_HIT_GROUPS = 3
KRAKEN_THREADS = 4

# Detection threshold per kingdom — mirrors KINGDOM_THRESHOLDS in stat/filter_runs.py
KINGDOM_THRESHOLD = {"Fungi": 0.5, "Oomycota": 0.5, "Nematoda": 1.0}

CONTROL_COLS = [
    "stratum", "Run", "mode", "BioSample", "BioProject",
    "named_host", "host_pct", "euk_pct", "fungi_pct", "oomycete_pct",
    "stat_pathogens", "co_infection_flag", "same_genus_secondary",
    "llm_treatment", "llm_study_setting", "expected",
]

RESULT_COLS = [
    "run_id", "species", "kingdom", "source_type", "source_id",
    "background_run", "background_species",
    "ratio_pct", "observed_pct", "detected",
    "detection_threshold", "replicate",
    "n_reads_total", "pct_classified", "top_hits",
]

TEMPLATE = {
    "_note": (
        "Replace SRR_PLACEHOLDER with real SRA accessions from pure-culture RNA-seq. "
        "For obligate biotrophs, use 'fastq' key with ART-simulated reads instead of 'run'."
    ),
    "background": {"run": "SRR_PLACEHOLDER", "species": "Triticum aestivum"},
    "mixtures": [
        {"name": "Botrytis cinerea",        "kingdom": "Fungi",    "run": "SRR_PLACEHOLDER"},
        {"name": "Sclerotinia sclerotiorum", "kingdom": "Fungi",    "run": "SRR_PLACEHOLDER"},
        {"name": "Fusarium oxysporum",       "kingdom": "Fungi",    "run": "SRR_PLACEHOLDER"},
        {"name": "Fusarium graminearum",     "kingdom": "Fungi",    "run": "SRR_PLACEHOLDER"},
        {"name": "Phytophthora infestans",   "kingdom": "Oomycota", "run": "SRR_PLACEHOLDER"},
        {
            "name":    "Puccinia striiformis",
            "kingdom": "Fungi",
            "fastq":   ["pst_sim1.fq.gz", "pst_sim2.fq.gz"],
        },
    ],
    "ratios":     [0.1, 0.5, 1.0, 5.0, 10.0, 25.0],
    "n_total":    500_000,
    "replicates": 1,
}


# ── ENA FTP helpers ────────────────────────────────────────────────────────────

def _ena_urls(run: str) -> list[str]:
    """Return ENA FTP FASTQ URLs for run. Empty list if unavailable."""
    import urllib.request
    api = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={run}&result=read_run&fields=fastq_ftp&format=tsv")
    req = urllib.request.Request(api, headers={"User-Agent": "crypt/insilico"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            body = r.read().decode()
    except Exception:
        return []
    lines = [l for l in body.strip().split("\n") if l]
    if len(lines) < 2:
        return []
    return [f"ftp://{p.strip()}" for p in lines[1].split("\t")[-1].split(";") if p.strip()]


def _stream_r1(url: str, n_reads: int, dest: Path) -> bool:
    """Stream n_reads records from gzipped ENA URL into dest (plain FASTQ, R1 only)."""
    cmd = ["bash", "-c",
           f'curl --silent --fail --max-time 600 "{url}" | gunzip -c | head -n {n_reads * 4}']
    try:
        with open(dest, "w") as out:
            subprocess.run(cmd, stdout=out, timeout=660, check=False)
        return dest.exists() and dest.stat().st_size > 0
    except Exception:
        return False


# ── Read acquisition ───────────────────────────────────────────────────────────

def _acquire_sra(run: str, n_reads: int, cache: Path) -> Path | None:
    """Download R1 from ENA into cache/{run}_r1.fastq. Returns path or None."""
    dest = cache / f"{run}_r1.fastq"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"    reusing {dest.name}")
        return dest
    urls = _ena_urls(run)
    if not urls:
        print(f"    ERROR: no ENA URLs for {run}")
        return None
    print(f"    streaming {urls[0]} → {dest.name} ({n_reads:,} reads)...", flush=True)
    ok = _stream_r1(urls[0], n_reads, dest)
    return dest if ok else None


def _acquire_fastq(fastq_paths: list[str], n_reads: int, cache: Path) -> Path | None:
    """Load R1 from a local FASTQ (plain or .gz). Returns path or None."""
    src = Path(fastq_paths[0])
    if not src.exists():
        print(f"    ERROR: {src} not found")
        return None
    dest = cache / "local_r1.fastq"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"    reusing {dest.name}")
        return dest
    if str(src).endswith(".gz"):
        cmd = ["bash", "-c", f'gunzip -c "{src}" | head -n {n_reads * 4}']
    else:
        cmd = ["bash", "-c", f'head -n {n_reads * 4} "{src}"']
    print(f"    loading {src.name} → {dest.name} ({n_reads:,} reads)...", flush=True)
    with open(dest, "w") as f:
        subprocess.run(cmd, stdout=f, timeout=300, check=False)
    return dest if (dest.exists() and dest.stat().st_size > 0) else None


# ── FASTQ record utilities ─────────────────────────────────────────────────────

def _count_records(path: Path) -> int:
    n = sum(1 for _ in open(path))
    return n // 4


def _take_records(src: Path, offset: int, n: int, dest: Path) -> int:
    """Write records [offset, offset+n) from src to dest. Returns records written."""
    start_line = offset * 4
    end_line   = start_line + n * 4
    written = 0
    with open(src) as fin, open(dest, "w") as fout:
        for i, line in enumerate(fin):
            if i < start_line:
                continue
            if i >= end_line:
                break
            fout.write(line)
            if (i + 1 - start_line) % 4 == 0:
                written += 1
    return written


# ── Kraken2 ────────────────────────────────────────────────────────────────────

def _run_kraken2(db: Path, r1: Path, report: Path, threads: int) -> bool:
    cmd = [
        "kraken2", "--db", str(db),
        "--report", str(report),
        "--confidence",        str(CONFIDENCE),
        "--minimum-hit-groups", str(MIN_HIT_GROUPS),
        "--threads", str(threads),
        "--output", "/dev/null",
        str(r1),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        return r.returncode == 0
    except Exception:
        return False


def _parse_report(report: Path) -> dict:
    species, pct_unclassified, n_reads = [], 0.0, 0
    try:
        with open(report) as f:
            lines = f.readlines()
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct, rank = float(parts[0]), parts[3].strip()
            name = parts[5].strip()
            if rank == "U":
                pct_unclassified = pct
            elif rank == "S" and pct > 0:
                species.append({"name": name, "pct": round(pct, 4)})
        u_line = next((l for l in lines if "\tU\t" in l), None)
        r_line = next((l for l in lines if "\tR\t1\t" in l), None)
        u_n = int(u_line.split("\t")[1]) if u_line else 0
        r_n = int(r_line.split("\t")[1]) if r_line else 0
        n_reads = u_n + r_n
    except Exception:
        pass
    return {
        "pct_classified": round(100.0 - pct_unclassified, 4),
        "n_reads":        n_reads,
        "species":        sorted(species, key=lambda x: -x["pct"]),
    }


def _observed_pct(species: list[dict], target: str) -> float:
    """Sum pct for all report entries where every word in target appears in the entry name."""
    words = set(target.lower().split())
    total = 0.0
    for s in species:
        if words.issubset(set(s["name"].lower().split())):
            total += s["pct"]
    return round(total, 4)


def _slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel",          metavar="PATH", help="Panel JSON config")
    ap.add_argument("--db",             metavar="PATH", help="Kraken2 database directory")
    ap.add_argument("--out-dir",        default=str(OUT_DIR_DEF))
    ap.add_argument("--tmp-dir",        default=str(SCRATCH_DEF))
    ap.add_argument("--kraken-threads", type=int, default=KRAKEN_THREADS)
    ap.add_argument("--write-template", action="store_true",
                    help="Write starter panel.json to cwd and exit")
    ap.add_argument("--dry-run",        action="store_true",
                    help="Print plan without downloading or classifying")
    args = ap.parse_args()

    if args.write_template:
        out = Path("panel.json")
        out.write_text(json.dumps(TEMPLATE, indent=2))
        print(f"Wrote {out.resolve()}")
        print("Fill in SRR accessions, then run:")
        print("  python kraken/control/insilico.py --panel panel.json --db /path/to/db")
        return

    if not args.panel:
        ap.error("--panel is required (or use --write-template)")
    if not args.db and not args.dry_run:
        ap.error("--db is required")

    with open(args.panel) as f:
        panel = json.load(f)

    ratios     = panel["ratios"]
    n_total    = panel["n_total"]
    replicates = panel.get("replicates", 1)
    bg_cfg     = panel["background"]
    max_ratio  = max(ratios)

    # reads needed from each source for all replicates
    n_path_needed = ceil(n_total * max_ratio / 100) * replicates + 500

    out_dir = Path(args.out_dir);  out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(args.tmp_dir);  tmp_dir.mkdir(parents=True, exist_ok=True)
    db_dir  = Path(args.db) if args.db else None

    print(f"Panel   : {len(panel['mixtures'])} species, {len(ratios)} ratios, "
          f"{replicates} replicate(s), {n_total:,} reads/mixture")
    if not args.dry_run:
        print(f"DB      : {db_dir}")
    print()

    # ── Acquire background reads (once) ───────────────────────────────────────
    bg_cache = tmp_dir / "background"
    bg_cache.mkdir(exist_ok=True)
    print(f"[background] {bg_cfg['species']} ({bg_cfg['run']})")
    if args.dry_run:
        bg_r1 = None
    else:
        bg_r1 = _acquire_sra(bg_cfg["run"], n_total, bg_cache)
        if bg_r1 is None:
            sys.exit("ERROR: failed to acquire background reads")
        n_bg_available = _count_records(bg_r1)
        print(f"    {n_bg_available:,} background records available")
    print()

    results, control_rows = [], []

    # ── Process each species ──────────────────────────────────────────────────
    for mix_cfg in panel["mixtures"]:
        sp_name = mix_cfg["name"]
        kingdom = mix_cfg["kingdom"]
        sp_slug = _slug(sp_name)
        threshold = KINGDOM_THRESHOLD.get(kingdom, 0.5)

        sp_cache = tmp_dir / sp_slug
        sp_cache.mkdir(exist_ok=True)

        print(f"[{sp_name}] ({kingdom}, detection threshold {threshold}%)")

        if args.dry_run:
            for ratio in ratios:
                n_path = round(n_total * ratio / 100)
                print(f"  would mix: {n_path:,} pathogen + {n_total - n_path:,} background @ {ratio}%")
            print()
            continue

        # acquire source reads
        if "fastq" in mix_cfg:
            src_r1   = _acquire_fastq(mix_cfg["fastq"], n_path_needed, sp_cache)
            src_type = "fastq"
            src_id   = mix_cfg["fastq"][0]
        else:
            src_r1   = _acquire_sra(mix_cfg["run"], n_path_needed, sp_cache)
            src_type = "sra"
            src_id   = mix_cfg["run"]

        if src_r1 is None:
            print(f"  ERROR: could not acquire reads — skipping\n")
            continue

        n_src_available = _count_records(src_r1)
        print(f"    {n_src_available:,} source records available (need {n_path_needed:,})")
        if n_src_available < n_path_needed:
            actual_reps = n_src_available // ceil(n_total * max_ratio / 100)
            if actual_reps < replicates:
                print(f"    WARNING: insufficient reads for {replicates} replicates; "
                      f"will do {max(1, actual_reps)}")
                replicates = max(1, actual_reps)

        for rep in range(replicates):
            path_offset = rep * ceil(n_total * max_ratio / 100)

            for ratio in ratios:
                n_path = round(n_total * ratio / 100)
                n_bg   = n_total - n_path
                if n_path == 0:
                    continue

                run_id  = f"INSILICO_{sp_slug}_{ratio:.1f}pct_rep{rep + 1}"
                mix_dir = sp_cache / run_id
                mix_dir.mkdir(exist_ok=True)

                path_tmp = mix_dir / "path.fastq"
                bg_tmp   = mix_dir / "bg.fastq"
                mix_r1   = mix_dir / "mix.fastq"
                report   = mix_dir / "report.txt"

                n_got = _take_records(src_r1, path_offset, n_path, path_tmp)
                _take_records(bg_r1, 0, n_bg, bg_tmp)

                # concatenate: pathogen first, then background
                with open(mix_r1, "w") as out:
                    for src in [path_tmp, bg_tmp]:
                        with open(src) as fin:
                            out.write(fin.read())
                path_tmp.unlink(missing_ok=True)
                bg_tmp.unlink(missing_ok=True)

                print(f"  {run_id}: {n_got} path + {n_bg} bg → kraken2...", flush=True)
                ok = _run_kraken2(db_dir, mix_r1, report, args.kraken_threads)
                mix_r1.unlink(missing_ok=True)

                if not ok:
                    print(f"    ERROR: kraken2 failed")
                    continue

                parsed       = _parse_report(report)
                observed_pct = _observed_pct(parsed["species"], sp_name)
                detected     = observed_pct >= threshold

                print(f"    expected={ratio:.1f}%  observed={observed_pct:.3f}%  "
                      f"detected={detected}  classified={parsed['pct_classified']:.1f}%")

                results.append({
                    "run_id":            run_id,
                    "species":           sp_name,
                    "kingdom":           kingdom,
                    "source_type":       src_type,
                    "source_id":         src_id,
                    "background_run":    bg_cfg["run"],
                    "background_species": bg_cfg["species"],
                    "ratio_pct":         ratio,
                    "observed_pct":      observed_pct,
                    "detected":          detected,
                    "detection_threshold": threshold,
                    "replicate":         rep + 1,
                    "n_reads_total":     parsed["n_reads"],
                    "pct_classified":    parsed["pct_classified"],
                    "top_hits":          json.dumps(parsed["species"][:5]),
                })

                control_rows.append({
                    "stratum":              "J",
                    "Run":                  run_id,
                    "mode":                 "insilico",
                    "BioSample":            "",
                    "BioProject":           "",
                    "named_host":           bg_cfg["species"],
                    "host_pct":             "",
                    "euk_pct":              ratio,
                    "fungi_pct":            ratio if kingdom == "Fungi" else 0.0,
                    "oomycete_pct":         ratio if kingdom == "Oomycota" else 0.0,
                    "stat_pathogens":       sp_name,
                    "co_infection_flag":    "insilico",
                    "same_genus_secondary": "False",
                    "llm_treatment":        "insilico",
                    "llm_study_setting":    "insilico",
                    "expected":             f"{sp_name} {ratio:.1f}% ({kingdom})",
                })

        print()

    # ── Write outputs ─────────────────────────────────────────────────────────
    res_path = out_dir / "results.tsv"
    with open(res_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=RESULT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {res_path}  ({len(results)} rows)")

    ctrl_path = out_dir / "control_rows.tsv"
    with open(ctrl_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CONTROL_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(control_rows)
    print(f"Wrote {ctrl_path}  ({len(control_rows)} rows)")

    if not args.dry_run and control_rows:
        print(f"\nTo append to control_runs.tsv:")
        print(f"  tail -n +2 {ctrl_path} >> {CONTROL_TSV}")


if __name__ == "__main__":
    main()
