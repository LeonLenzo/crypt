#!/usr/bin/env python3
"""fetch_genomes.py — download GENOMIC FASTA for the benchmark sample.

Step 2 of the predict evaluation. Reads `data/silver/predict_sample.tsv`, writes under
`data/bronze/genomes/<accession>/`, and touches nothing else.

## Why this is a separate download and not a reuse

db_v4 holds **CDS only**, deliberately: the cohort is 100% RNA-seq, so the transcriptome is
what reads come from and genomic references are wrong here rather than merely bigger. But a
gene predictor takes a genome — you cannot predict genes from the CDS you are trying to
produce. So this is new bytes, not a reconciliation, and there is no Acacia pass because these
accessions were never downloaded in genomic form.

Scale check before running: fungal genomes are ~35-60 Mbp, so ~120 assemblies is ~7 GB.
That is small enough to redo, which is the point of benchmarking on a sample rather than
committing to all 2,138 unannotated species.

**Presence is a non-empty file above a floor**, the same rule as `02_build/fetch_cds.py`, for
the same reason: /scratch's purge removes files and leaves directories, so a directory check
counts an empty skeleton as held. The floor here is 1 MB rather than 100 kB because a genomic
assembly is an order of magnitude larger than a CDS set.

Run on Setonix inside a job, never on a login node:
    sbatch 02_classify/slurm/predict_fetch.slurm
    python 02_classify/01_refs/predict/fetch_genomes.py --inventory
"""

import argparse
import csv
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

MODULE = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(MODULE))
sys.path.insert(0, str(MODULE.parent))

from _layout import GENOMES, LOGS, PREDICT_SAMPLE, SILVER  # noqa: E402
from _util import _Tee, atomic_write, link_latest, make_log_dir  # noqa: E402

MIN_FASTA_BYTES = 1_000_000
DERIVED_SUFFIXES = ("_clean.fna", "_combined.fna", ".tagged.fna")
GENOME_MANIFEST = SILVER / "predict_genomes.tsv"


def genome_file(acc_dir: Path):
    """Largest non-derived .fna at or above the floor, else None."""
    if not acc_dir.is_dir():
        return None
    c = [f for f in acc_dir.glob("**/*.fna")
         if f.is_file() and not f.name.endswith(DERIVED_SUFFIXES)
         and f.stat().st_size >= MIN_FASTA_BYTES]
    return max(c, key=lambda f: f.stat().st_size) if c else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample", default=None, help=f"default {PREDICT_SAMPLE}")
    ap.add_argument("--genomes-dir", default=None, help=f"default {GENOMES}")
    ap.add_argument("--workers", type=int, default=8,
                    help="8, not 16: measured on the CDS pull, past ~8 NCBI throttles and the "
                         "retries cost more than the concurrency returns")
    ap.add_argument("--arm", default=None, choices=("annotated", "unannotated"),
                    help="restrict to one arm")
    ap.add_argument("--limit", type=int, default=0,
                    help="fetch at most N (for an install probe, which needs one genome)")
    ap.add_argument("--inventory", action="store_true", help="report and exit, move nothing")
    args = ap.parse_args()

    log_dir = make_log_dir(LOGS)
    log = _Tee(log_dir / "predict_fetch_genomes.log")
    link_latest(LOGS, log_dir)
    try:
        sample_path = Path(args.sample) if args.sample else PREDICT_SAMPLE
        gdir = Path(args.genomes_dir) if args.genomes_dir else GENOMES
        if not sample_path.exists():
            sys.exit(f"no sample at {sample_path} — run 01_refs/predict/sample.py first")
        rows = list(csv.DictReader(open(sample_path, newline=""), delimiter="\t"))
        if args.arm:
            rows = [r for r in rows if r["arm"] == args.arm]
        if args.limit:
            rows = rows[:args.limit]

        have = {r["accession"]: f for r in rows
                if (f := genome_file(gdir / r["accession"]))}
        todo = [r for r in rows if r["accession"] not in have]
        print(f"sample: {len(rows)} assemblies  ({sample_path})")
        print(f"  genomic FASTA already on disk: {len(have)}")
        print(f"  to download:                   {len(todo)}")
        skeletons = [r["accession"] for r in todo if (gdir / r["accession"]).is_dir()]
        if skeletons:
            print(f"  of which EMPTY directories (purged, not held): {len(skeletons)}")
        if args.inventory:
            print("\n--inventory: nothing moved.")
            return

        if todo:
            from _assembly import download_cds  # noqa: E402  late: NCBI-only path
            gdir.mkdir(parents=True, exist_ok=True)
            print(f"\nDownloading {len(todo)} genomic assemblies ({args.workers} workers)")
            done = 0
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                futs = {ex.submit(download_cds, r["accession"], gdir / r["accession"],
                                  "genome"): r for r in todo}
                for fut in as_completed(futs):
                    r = futs[fut]
                    try:
                        fut.result()
                    except Exception as exc:                      # noqa: BLE001
                        print(f"  {r['accession']}  FAILED: {type(exc).__name__}: {exc}")
                    f = genome_file(gdir / r["accession"])
                    if f:
                        have[r["accession"]] = f
                    done += 1
                    if done % 20 == 0 or done == len(todo):
                        print(f"  [{done}/{len(todo)}]  ok={len(have)}", flush=True)

        cols = ["arm", "accession", "taxid", "organism_name", "assembly_level",
                "has_truth", "focal", "bytes", "path"]
        missing = []
        with atomic_write(GENOME_MANIFEST) as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(cols)
            for r in rows:
                f = have.get(r["accession"])
                if not f:
                    missing.append(r)
                w.writerow([r["arm"], r["accession"], r["taxid"], r["organism_name"],
                            r["assembly_level"], r["has_truth"], r["focal"],
                            f.stat().st_size if f else 0,
                            str(f.relative_to(gdir)) if f else ""])
        print(f"\nmanifest → {GENOME_MANIFEST}")
        tally = Counter(r["arm"] for r in rows if r["accession"] in have)
        print("── summary ──────────────────────────────────────────────")
        for k, v in sorted(tally.items()):
            print(f"  {k:12s} {v:4d}")
        total_b = sum(f.stat().st_size for f in have.values())
        print(f"  {'total':12s} {len(have):4d} of {len(rows)}  ({total_b/1e9:.1f} GB)")

        if missing:
            print(f"\n{len(missing)} assemblies have no genomic FASTA:")
            for r in missing[:10]:
                print(f"  {r['accession']:18s} {r['organism_name']}")
            print("  Re-run to retry; the download is resumable. A persistent failure is a "
                  "real NCBI withdrawal and the sample should be re-drawn without it rather "
                  "than silently shrinking one stratum.")
            sys.exit(1)
    finally:
        log.close()


if __name__ == "__main__":
    main()
