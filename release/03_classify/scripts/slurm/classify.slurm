#!/bin/bash
#SBATCH --job-name=classify
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=release/03_classify/scripts/slurm/logs/classify_%j.out
#SBATCH --error=release/03_classify/scripts/slurm/logs/classify_%j.err
#
# Classify one or more BioProjects against one arm of the ladder.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/03_classify/scripts/slurm/classify.slurm --arm db1p --project PRJEB39201
#   sbatch release/03_classify/scripts/slurm/classify.slurm --arm db2 --project PRJNA1 PRJNA2
#
# ARGUMENTS ARE FLAGS, NOT ENVIRONMENT VARIABLES. Setonix defaults sbatch to --export=NONE, so
# `ARM=db1p sbatch ...` reaches the job with ARM unset; a build job died in 9 s on exactly that.
# Everything here is forwarded to classify.py verbatim, so the flags are the script's own and
# there is one place to read them.
#
# --mem=100G, and the request is the CPU BILL. Memory bills about 1.9 GB per physical core, so
# 100G bills roughly 53 of the node's 64 cores whether they are busy or not. It is sized to the
# hash and no further: 57.5 GB (db1p) or 59.8 GB (db2) staged into /dev/shm is charged to the
# job's cgroup, and the workers' working sets sit on top. db_v4's equivalent job ran at 65.5 GB
# against a 41 GB hash, so the headroom scales to about 85 GB here. A request sized for db_v4
# does not fit.
#
# GROUP THE SMALL PROJECTS. Staging the hash costs a few minutes and is per JOB, not per
# project, so a 4-run project in its own job spends the stage on two minutes of work. Pass
# several accessions to one job.
#
# --memory-mapping is used here and only because the hash is in /dev/shm, which is RAM. Against
# Lustre it is catastrophic and measured: 13 minutes, zero reports, kraken2 at 1.6% CPU, because
# every hash lookup becomes a filesystem RPC. Never pass --mmap without the staging below.
#
# 8 workers x 8 threads = 64 concurrent threads, matching the 64 PHYSICAL cores (sacct reports
# 128, which is the SMT logical count). Workers are cheap to add because they mmap ONE shared
# /dev/shm copy of the hash rather than loading their own.
#
# -e is deliberately NOT set: the exit code is captured, the manifest is published, and the
# cleanup trap removes the staged hash. Under -e a failed run would skip all three.

set -uo pipefail
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p release/03_classify/scripts/slurm/logs
module load python/3.11.6

# kraken2 PINNED BY PATH, never `spack load`: spack resolves against fl3 rather than the job's
# account, and in a non-interactive shell `spack` is not on PATH after `module load spack` at
# all. zen3 for compute nodes, zen2 as the fallback.
K2_ZEN3="/software/projects/pawsey1168/llenzo/setonix/2025.08/software/linux-sles15-zen3/gcc-14.2.0/kraken2-2.1.2-qk2lek64msmwg57n2fowuf5b7etidqaj/bin"
K2_ZEN2="/software/projects/pawsey1168/llenzo/setonix/2025.08/software/linux-sles15-zen2/gcc-14.2.0/kraken2-2.1.2-errlhphjipv3gwjwefhruqqogybieteb/bin"
if "$K2_ZEN3/kraken2" --version >/dev/null 2>&1; then export PATH="$K2_ZEN3:$PATH"
elif "$K2_ZEN2/kraken2" --version >/dev/null 2>&1; then
    echo "note: zen3 kraken2 would not run here, falling back to zen2"
    export PATH="$K2_ZEN2:$PATH"
else echo "ABORT: no runnable kraken2 found at either pinned path."; exit 1; fi
echo "kraken2: $(which kraken2)  $(kraken2 --version 2>&1 | head -1)"

# --arm decides which hash is staged, so it is read back out of the forwarded arguments.
ARM=""
prev=""
for a in "$@"; do
    [ "$prev" = "--arm" ] && ARM="$a"
    case "$a" in --arm=*) ARM="${a#--arm=}" ;; esac
    prev="$a"
done
case "$ARM" in
  db1p|db2) ;;
  *) echo "usage: sbatch classify.slurm --arm {db1p|db2} --project PRJ [PRJ ...]"; exit 2 ;;
esac

DB="release/02_reference/data/db/$ARM"
SHM="/dev/shm/${ARM}_${SLURM_JOB_ID}"
[ -f "$DB/hash.k2d" ] || { echo "ABORT: no index at $DB"; exit 1; }

echo "host     $(hostname)"
echo "arm      $ARM"
echo "started  $(date -Is)"

SECONDS=0
mkdir -p "$SHM"
cleanup() { rm -rf "$SHM"; }
trap cleanup EXIT
cp "$DB"/hash.k2d "$DB"/opts.k2d "$DB"/taxo.k2d "$SHM"/ || { echo "shm stage FAILED"; exit 1; }
echo "staged $(du -sh --apparent-size "$SHM" | cut -f1) into RAM in ${SECONDS}s"

# The tmp dir is PER JOB, and that is a correctness requirement rather than tidiness: two
# classify jobs overlapping on one run would otherwise have two kraken2 processes writing one
# .partial file and both renaming it over the real report. The result is a full-size report of
# interleaved output, which the size check cannot catch.
SECONDS=0
python3 -u release/03_classify/scripts/classify.py \
    --db "$SHM" --mmap \
    --workers "${WORKERS:-8}" --threads "${THREADS:-8}" \
    --tmp-dir "$MYSCRATCH/kraken_tmp/$SLURM_JOB_ID" "$@"
rc=$?

# The manifest describes reports only this machine holds, so this machine commits it.
tools/repo_sync.sh publish \
    "03_classify: ${ARM} reports manifest (job ${SLURM_JOB_ID:-?})" \
    "release/03_classify/data/raw/reports_manifest_${ARM}.tsv"
prc=$?
[ "$rc" -eq 0 ] && rc=$prc

echo "exit     $rc after ${SECONDS}s"
echo "finished $(date -Is)"
exit $rc
