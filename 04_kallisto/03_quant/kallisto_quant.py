#!/usr/bin/env python3
"""kallisto_quant.py — pseudoalign each selected run against its host index.

Two passes per run, because one cannot give both facts 04_confirm needs:

  quant  kallisto quant -> abundance.tsv (est_counts, tpm per transcript). The EM's
         apportionment, which is what aggregates to species abundance and gene breadth.
  bus    kallisto bus -x Bulk -> output.bus + matrix.ec. The equivalence class each read fell
         in, which is the only way to split uniquely-assigned from shared reads.

The second pass is not optional and not redundant. `kallisto quant` in 0.50.1 has no
equivalence-class output at all: `--dump-eq` and `kallisto pseudo` were both removed after
0.46, and `--pseudobam` is gone from quant. Verified against the installed binary, not assumed.
`quant-tcc` could in principle collapse this to one pass (bus -> TCC -> quant-tcc), but that
means hand-rolling kallisto's TCC format, and a silent format mismatch would corrupt the
abundances that everything downstream rests on. Two passes is the cheap, safe option: kallisto
is fast and the selected set is only 0.49 TB.

bustools is not installed and not in Pawsey's spack repo, so 04_confirm parses output.bus
directly. The format is trivial and stable: a header then fixed 32-byte records.

Reads come from 03_kraken/03_select/data/reads on Setonix scratch, under a 21-day atime purge.
Run freshen.py --report first if it has been a while.

    sbatch 04_kallisto/slurm/kallisto_quant.slurm
    python 04_kallisto/03_quant/kallisto_quant.py --list
    python 04_kallisto/03_quant/kallisto_quant.py --run SRR123456
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import IDX_DIR, QUANT_DIR, READS_ROOT, STRATA, slug, step_log

# Required by kallisto for single-end data, which cannot estimate the fragment distribution
# from the reads. Matches the pilot in 03_kraken/utilities/pilot/kallisto/run_kallisto.py.
SINGLE_FRAG_LEN, SINGLE_FRAG_SD = 200, 20


def targets(host=None, run=None):
    """[(host, run, layout, [abs read paths])], one entry per run.

    Deduplicated by run: a run carries every detection it holds, so quant is per run, not per
    detection. 174 runs for 200 detections in the current stratified set.
    """
    seen, out = set(), []
    for r in csv.DictReader(open(STRATA), delimiter="\t"):
        if r["run"] in seen or (host and r["host"] != host) or (run and r["run"] != run):
            continue
        seen.add(r["run"])
        reads = [READS_ROOT / p for p in r["read_files"].split(",") if p]
        out.append((r["host"], r["run"], r["layout"], reads))
    return out


def _done(out_dir, want_bus):
    """True if this run's output is complete AND readable.

    Integrity, not existence. The 2026-09-14 truncated-gzip incident came from a resume check
    that only asked whether a file was present, so a half-written file counted as success and
    its source was deleted. Here: run_info.json must parse and carry n_processed, and the bus
    outputs must be non-empty.
    """
    info = out_dir / "run_info.json"
    if not info.exists():
        return False
    try:
        if json.load(open(info)).get("n_processed") is None:
            return False
    except (json.JSONDecodeError, OSError):
        return False
    if not (out_dir / "abundance.tsv").exists():
        return False
    if want_bus:
        for f in ("output.bus", "matrix.ec", "transcripts.txt"):
            p = out_dir / "bus" / f
            if not p.exists() or p.stat().st_size == 0:
                return False
    return True


def _run(cmd):
    print(f"  $ {' '.join(str(c) for c in cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        tail = (r.stderr or r.stdout or "").strip().splitlines()[-4:]
        print(f"  !! exit {r.returncode}: {' | '.join(tail)}", file=sys.stderr)
    return r.returncode == 0


def quant_one(host, run, layout, reads, threads=8, with_bus=True, bootstrap=0, force=False):
    idx = IDX_DIR / f"{slug(host)}.idx"
    out_dir = QUANT_DIR / slug(host) / run
    if not idx.exists():
        print(f"  {run}: MISSING index {idx.name}; build it first (02_build)", file=sys.stderr)
        return False
    missing = [p for p in reads if not p.exists()]
    if missing:
        print(f"  {run}: reads gone from scratch ({missing[0].name}); purged? "
              f"re-pull from run_list.tsv", file=sys.stderr)
        return False
    if not force and _done(out_dir, with_bus):
        print(f"  {run}: already quanted, skipping")
        return True

    out_dir.mkdir(parents=True, exist_ok=True)
    layout_args = ([] if layout == "paired"
                   else ["--single", "-l", str(SINGLE_FRAG_LEN), "-s", str(SINGLE_FRAG_SD)])

    print(f"{run} ({host}, {layout}, {len(reads)} files)")
    cmd = ["kallisto", "quant", "-i", str(idx), "-o", str(out_dir), "-t", str(threads)]
    if bootstrap:
        cmd += ["-b", str(bootstrap)]
    cmd += layout_args + [str(p) for p in reads]
    if not _run(cmd):
        return False

    try:
        info = json.load(open(out_dir / "run_info.json"))
        print(f"  quant: {info.get('n_processed', 0):,} reads, "
              f"{info.get('p_pseudoaligned', 0)}% pseudoaligned")
    except (json.JSONDecodeError, OSError):
        print("  !! run_info.json unreadable after quant", file=sys.stderr)
        return False

    if with_bus:
        bus_dir = out_dir / "bus"
        bus_dir.mkdir(exist_ok=True)
        cmd = ["kallisto", "bus", "-i", str(idx), "-o", str(bus_dir),
               "-x", "Bulk", "-t", str(threads), "--num"]
        if layout == "paired":
            cmd.append("--paired")
        cmd += [str(p) for p in reads]
        if not _run(cmd):
            return False
        print(f"  bus: {(bus_dir / 'output.bus').stat().st_size / 1e6:.1f} MB, "
              f"{sum(1 for _ in open(bus_dir / 'matrix.ec')):,} equivalence classes")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", help="restrict to one host")
    ap.add_argument("--run", help="one run accession")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bootstrap", type=int, default=0,
                    help="bootstrap samples; 04_confirm's three facts do not need them, so 0")
    ap.add_argument("--no-bus", action="store_true",
                    help="skip the equivalence-class pass (then unique reads are unavailable)")
    ap.add_argument("--force", action="store_true", help="re-quant even if output is complete")
    ap.add_argument("--list", action="store_true", help="show what would be quanted and exit")
    args = ap.parse_args()

    todo = targets(args.host, args.run)
    if args.list:
        print(f"{len(todo)} runs over {len({h for h, _, _, _ in todo})} hosts")
        for host, run, layout, reads in todo:
            idx = IDX_DIR / f"{slug(host)}.idx"
            print(f"  {run:12s} {host:24s} {layout:7s} "
                  f"idx={'ok' if idx.exists() else 'MISSING':7s} "
                  f"reads={'ok' if all(p.exists() for p in reads) else 'MISSING'}")
        return

    log = step_log(__file__, "kallisto_quant")
    try:
        ok = sum(quant_one(h, r, l, rd, args.threads, not args.no_bus,
                           args.bootstrap, args.force)
                 for h, r, l, rd in todo)
        print(f"\n{ok}/{len(todo)} runs quanted")
    finally:
        log.close()
    sys.exit(0 if ok == len(todo) else 1)


if __name__ == "__main__":
    main()
