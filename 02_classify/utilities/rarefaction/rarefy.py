#!/usr/bin/env python3
"""rarefy.py — at what sequencing depth does a detection stop being detectable?

The question this answers is not academic. The cohort's libraries span 700x in the quantity
detection actually depends on, and the spread is CONFOUNDED WITH CROP:

    PRJDB7234     1,557 runs  SINGLE  51 nt  0.057 Gb  ->  0.019 Gkmer   Rice
    PRJNA726716     104 runs  SINGLE  51 nt  0.109 Gb  ->  0.036 Gkmer   Rice
    PRJEB39201      538 runs  PAIRED 202 nt  4.533 Gb  ->  3.007 Gkmer   Wheat
    PRJNA1181034     12 runs  PAIRED 300 nt 17.123 Gb  -> 13.242 Gkmer   Wheat

Those two rice projects are 1,661 of the 2,784 rice runs, so 60% of rice sits at the bottom of
that range. If rice returns a low co-infection rate, depth is a sufficient explanation and we
could not tell it apart from biology. A depth cut-off chosen by eye would also BE a rice
filter: at 0.2 Gb it drops 1,551 rice runs and 18 maize. So the threshold has to be measured,
not guessed, and this is the measurement.

## Gkmer, not Gb, and not reads

At k=35 a read of length L yields L-34 k-mers, so k-mer yield PER BASE is (L-34)/L: 0.33 at
51 nt, 0.66 at 101 nt, 0.77 at 150 nt. Read length is therefore worth a 2.3x factor on its own,
independent of depth, and a library must be ranked on bases x yield. That composite is `Gkmer`
and it is what the rungs are set in. Ranking on layout instead gets it plainly wrong: a
paired-only filter would discard 301 single-end runs deeper than 2 Gb while keeping 311 paired
runs under 1 Gb.

## Design

FIVE RUNS, 22 DETECTIONS. One Kraken2 report carries every taxon in the sample, so rarefying
one run measures all of its detections at once. The five were chosen to be deep enough for 170x
headroom and to span coverage from 1.0% to 65%, with 7 detections under 5% -- those sit nearest
the floor and cross it first, so they carry most of the information.

TWO ARMS, matched on Gkmer so the only difference is library SHAPE:

    native    subsample only, keeping the run's own read length and pairing
    shallow   R1 only, trimmed to 51 nt, then subsampled to the SAME Gkmer

If `shallow` loses detections that `native` keeps at equal Gkmer, then read length matters
beyond k-mer count -- which is expected, because `--minimum-hit-groups 3` is harder to satisfy
with 17 k-mers per read than with 117, and a short read gives the LCA less to work with. That
is exactly PRJDB7234's situation and is the thing a Gb-only analysis would miss.

SEED REPLICATES at the two shallowest rungs, where subsampling variance is largest and a single
draw could put a detection either side of the floor by luck.

## Subsampling

Bernoulli selection on the record index in awk: one streaming pass, no memory, and because both
mates are processed by the same awk build with the same seed and one rand() call per record, the
same records are drawn from both files. There is no seqtk or seqkit on Setonix and this needs
neither.

It is NOT `head -n`. Reads in an SRA-derived fastq are in flowcell order, so the first N reads
are a corner of the flowcell rather than a sample of the library, and tile-edge quality effects
would be read as a depth effect.

## Cost

Kraken2 reloads db_v4's 41 GB hash per invocation, and at these depths the LOAD dominates the
classification. So the hash is staged once into /dev/shm and every worker mmaps that one copy.
CLAUDE.md says do not use --memory-mapping, measured at 13 minutes for zero reports; that
finding was mmap against LUSTRE, where each page fault becomes a filesystem RPC. /dev/shm is
RAM and is the opposite case. The note stands for the case it measured.

    python 02_classify/utilities/rarefaction/rarefy.py --select
    python 02_classify/utilities/rarefaction/rarefy.py --task 3
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _layout import DB, DB_INSPECT, GOLD, READS, READS_MANIFEST, SILVER  # noqa: E402
from _reports import coverage, load_inspect, parse_report  # noqa: E402

HERE = Path(__file__).resolve().parent
WORKLIST = SILVER / "rarefy_worklist.tsv"
CURVE = GOLD / "rarefaction_curve.tsv"
K = 35

# Rungs in Gkmer, set at the REAL project medians so a result reads directly as a statement
# about a project rather than about an abstract depth.
RUNGS = [0.019, 0.036, 0.10, 0.21, 0.50, 1.00, 2.00, 3.00]
SHALLOW_LEN = 51          # PRJDB7234's read length
REPLICATE_RUNGS = 2       # how many of the lowest rungs get seed replicates
N_SEEDS = 3
SEED0 = 20261019

# db_v4's operating point, unchanged from the published db_v3 settings so the curve is about
# depth and nothing else.
CONFIDENCE = "0.15"
MIN_HIT_GROUPS = "3"

WORK_COLS = ["task", "run", "arm", "rung_gkmer", "seed", "native_gkmer", "frac",
             "layout", "read_len", "bases"]
CURVE_COLS = ["run", "arm", "rung_gkmer", "seed", "taxid", "species", "rank",
              "reads_clade", "distinct", "db_clade", "coverage", "native_coverage"]


def kmers_per_base(read_len: float, layout: str) -> float:
    """(L-K+1)/L, with L the length of ONE read.

    ENA's read_len is bases per spot, so for a paired library it is both mates and the
    per-read length is half. Getting this wrong doubles the apparent yield of every paired
    library, which would make the shallow single-end projects look relatively better than
    they are.
    """
    L = read_len / 2 if layout == "PAIRED" else read_len
    return max(0.0, L - K + 1) / L if L > K else 0.0


def select(args) -> None:
    """Write the worklist: one row per (run, arm, rung, seed) to classify."""
    depth = {r["Run"]: r for r in csv.DictReader(open(SILVER / "run_depth.tsv"),
                                                 delimiter="\t")}
    runs = [r.strip() for r in args.runs.split(",") if r.strip()]
    rows, task = [], 0
    for run in runs:
        d = depth.get(run)
        if not d:
            sys.exit(f"{run} has no depth record. Run fetch_depth.py first.")
        bases, rlen, lay = int(d["bases"]), float(d["read_len"]), d["layout"]
        native_native = bases * kmers_per_base(rlen, lay) / 1e9
        # The shallow arm keeps R1 only and trims to 51 nt, so its yield per ORIGINAL base is
        # halved for pairing and recomputed for the shorter read.
        # The shallow arm keeps ONE mate and truncates it to 51 nt, so of the `rlen` bases
        # per spot it retains min(51, per_read). For a 2x75 library that is 51/150 = 34% of
        # the bases; for a library already at 51 nt single-end it is all of them.
        per_read = rlen / 2 if lay == "PAIRED" else rlen
        kept_len = min(SHALLOW_LEN, per_read)
        native_shallow = bases * (kept_len / rlen) * kmers_per_base(kept_len, "SINGLE") / 1e9
        for arm, native in (("native", native_native), ("shallow", native_shallow)):
            for i, rung in enumerate(RUNGS):
                if rung > native * 0.98:
                    continue              # cannot rarefy UP; skip rather than fake it
                seeds = range(N_SEEDS) if i < REPLICATE_RUNGS else range(1)
                for s in seeds:
                    task += 1
                    rows.append({"task": task, "run": run, "arm": arm,
                                 "rung_gkmer": rung, "seed": SEED0 + s,
                                 "native_gkmer": f"{native:.4f}",
                                 "frac": f"{rung / native:.6f}",
                                 "layout": lay, "read_len": f"{rlen:.1f}",
                                 "bases": bases})
    SILVER.mkdir(parents=True, exist_ok=True)
    tmp = WORKLIST.with_suffix(".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=WORK_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(WORKLIST)
    print(f"-> {WORKLIST.name}: {len(rows)} classifications")
    for run in runs:
        rr = [r for r in rows if r["run"] == run]
        nat = [r for r in rr if r["arm"] == "native"]
        print(f"  {run:13s} native {float(nat[0]['native_gkmer']) if nat else 0:5.2f} Gkmer, "
              f"{len(rr)} tasks")
    print(f"\nsbatch --array=1-{len(rows)} 02_classify/slurm/rarefaction.slurm")


def subsample(src: Path, dst: Path, frac: float, seed: int, trim: int | None) -> None:
    """Bernoulli-sample whole records, optionally truncating sequence and quality."""
    trim_expr = (f'if(n==2||n==4) $0=substr($0,1,{trim}); ' if trim else "")
    prog = (f'BEGIN{{srand({seed})}} '
            f'(NR-1)%4==0{{keep=(rand()<{frac})}} '
            f'keep{{n=((NR-1)%4)+1; {trim_expr}print}}')
    with open(dst, "wb") as out:
        z = subprocess.Popen(["gzip", "-cd", str(src)], stdout=subprocess.PIPE)
        a = subprocess.Popen(["awk", prog], stdin=z.stdout, stdout=out)
        z.stdout.close()
        a.communicate()
        z.wait()
        if a.returncode != 0:
            raise RuntimeError(f"awk failed on {src}")


def classify(run: str, files: list, db: Path, report: Path, threads: int, mmap: bool) -> None:
    cmd = ["kraken2", "--db", str(db), "--threads", str(threads),
           "--confidence", CONFIDENCE, "--minimum-hit-groups", MIN_HIT_GROUPS,
           # Without this the report has 6 columns and no distinct count, so coverage -- the
           # entire output of this experiment -- would be None for every row.
           "--report-minimizer-data", "--report", str(report), "--output", "/dev/null"]
    if mmap:
        cmd.append("--memory-mapping")
    if len(files) == 2:
        cmd.append("--paired")
    cmd += [str(f) for f in files]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not report.exists():
        raise RuntimeError(f"kraken2 failed for {run}: {r.stderr[-800:]}")


def run_task(args) -> None:
    work = {int(r["task"]): r for r in csv.DictReader(open(WORKLIST), delimiter="\t")}
    t = work.get(args.task)
    if not t:
        sys.exit(f"no task {args.task} in {WORKLIST.name} ({len(work)} tasks)")
    man = {r["Run"]: r for r in csv.DictReader(open(READS_MANIFEST), delimiter="\t")}
    m = man.get(t["run"])
    if not m:
        sys.exit(f"{t['run']} is not in {READS_MANIFEST.name}")
    src_dir = READS / m["BioProject"]
    files = [src_dir / n for n in m["files"].split(";")]
    missing = [f for f in files if not f.exists()]
    if missing:
        sys.exit(f"missing reads: {missing}")

    # The shallow arm is R1 only, so a paired run contributes one file.
    if t["arm"] == "shallow":
        files = [f for f in files if "_2.fastq" not in f.name][:1]
        trim = SHALLOW_LEN
    else:
        trim = None

    tmp = Path(args.tmp_dir) / f"rarefy_{args.task}"
    tmp.mkdir(parents=True, exist_ok=True)
    frac = min(1.0, float(t["frac"]))
    subs = []
    t0 = time.time()
    for f in files:
        d = tmp / f.name.replace(".fastq.gz", ".fastq")
        subsample(f, d, frac, int(t["seed"]), trim)
        subs.append(d)
    sub_bytes = sum(p.stat().st_size for p in subs)
    print(f"  subsampled {len(subs)} file(s), frac={frac:.5f}, "
          f"{sub_bytes/1e6:.1f} MB in {time.time()-t0:.0f}s", flush=True)

    report = tmp / "report.tsv"
    classify(t["run"], subs, Path(args.db), report, args.threads, args.mmap)

    insp = load_inspect(Path(args.inspect))
    out = []
    for rec in parse_report(report):
        if rec["rank"] != "S" or not rec["reads_clade"]:
            continue
        cov = coverage(rec["distinct"], rec["taxid"], insp)
        if cov is None:
            continue
        out.append({"run": t["run"], "arm": t["arm"], "rung_gkmer": t["rung_gkmer"],
                    "seed": t["seed"], "taxid": rec["taxid"], "species": rec["name"],
                    "rank": rec["rank"], "reads_clade": rec["reads_clade"],
                    "distinct": rec["distinct"],
                    "db_clade": insp[rec["taxid"]]["clade"],
                    "coverage": f"{cov:.6f}", "native_coverage": ""})
    part = SILVER / "rarefy_parts"
    part.mkdir(parents=True, exist_ok=True)
    pf = part / f"task_{args.task:04d}.tsv"
    with open(pf, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CURVE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    for p in subs:
        p.unlink(missing_ok=True)
    report.unlink(missing_ok=True)
    print(f"-> {pf.name}: {len(out)} species rows, total {time.time()-t0:.0f}s")


def collect(args) -> None:
    part = SILVER / "rarefy_parts"
    rows = []
    for p in sorted(part.glob("task_*.tsv")):
        rows += list(csv.DictReader(open(p), delimiter="\t"))
    if not rows:
        sys.exit(f"no parts in {part}")
    GOLD.mkdir(parents=True, exist_ok=True)
    tmp = CURVE.with_suffix(".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CURVE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(CURVE)
    print(f"-> {CURVE.name}: {len(rows)} rows from {len(list(part.glob('task_*.tsv')))} tasks")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--select", action="store_true", help="write the worklist")
    ap.add_argument("--runs", default="ERR15005517,ERR15005540,ERR3170251,"
                                      "ERR1607752,ERR1607774",
                    help="runs to rarefy (chosen for depth headroom and coverage spread)")
    ap.add_argument("--task", type=int, help="run one worklist task")
    ap.add_argument("--collect", action="store_true", help="merge the parts")
    ap.add_argument("--db", default=str(DB / "db_v4"))
    ap.add_argument("--inspect", default=str(DB_INSPECT))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--mmap", action="store_true",
                    help="--memory-mapping. Correct ONLY against /dev/shm, never Lustre")
    ap.add_argument("--tmp-dir", default=os.environ.get("MYSCRATCH", "/tmp") + "/rarefy")
    args = ap.parse_args()

    if args.select:
        select(args)
    elif args.collect:
        collect(args)
    elif args.task:
        run_task(args)
    else:
        ap.error("one of --select, --task or --collect")


if __name__ == "__main__":
    main()
