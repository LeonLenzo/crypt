#!/bin/bash
#SBATCH --job-name=build_db
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=12:00:00
#SBATCH --output=release/02_reference/scripts/slurm/logs/build_db_%j.out
#SBATCH --error=release/02_reference/scripts/slurm/logs/build_db_%j.err
#SBATCH --dependency=singleton
#
# Build one Kraken2 index from the CDS on disk. Fetches nothing.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/02_reference/scripts/slurm/build_db.slurm --name db1p
#   sbatch ... --name db2 --extra-cds data/curation/predicted_cds.tsv
#   sbatch ... --name db1p --build-only
#
# ARGUMENTS ARE FLAGS, NOT ENVIRONMENT VARIABLES, and that is not a style choice: Setonix
# defaults sbatch to --export=NONE, so `NAME=db1p sbatch ...` reaches the job with NAME unset.
# Job 53297249 died in 9 s on exactly that. Everything here is forwarded to build_db.py
# verbatim, so the flags are the script's own and there is one place to read them.
#
# PARTITION=work and a whole node: db_v4 took 5 h 40 m and peaked near 106 GB of hash for
# ~10 B minimizers, and --build is the only step that needs a compute node at all. The wall is
# 12 h rather than 24 h because the proven run fits twice over and a shorter wall queues
# sooner; --build-only resumes without redoing the library if it ever does run out.
#
# --dependency=singleton so two arms of the ladder cannot build into the same library at once.
# -e is deliberately NOT set: the exit code is captured and reported, and under -e a failed
# build would exit before the publish hook and before rc is read.
set -uo pipefail
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p release/02_reference/scripts/slurm/logs
module load python/3.11.6

# kraken2 PINNED BY PATH, never `spack load`: spack resolves against fl3 rather than the
# job's account, and in a non-interactive shell `spack` is not on PATH after `module load
# spack` at all. zen3 for compute nodes, zen2 as the fallback.
K2_ZEN3="/software/projects/pawsey1168/llenzo/setonix/2025.08/software/linux-sles15-zen3/gcc-14.2.0/kraken2-2.1.2-qk2lek64msmwg57n2fowuf5b7etidqaj/bin"
K2_ZEN2="/software/projects/pawsey1168/llenzo/setonix/2025.08/software/linux-sles15-zen2/gcc-14.2.0/kraken2-2.1.2-errlhphjipv3gwjwefhruqqogybieteb/bin"
if "$K2_ZEN3/kraken2-build" --version >/dev/null 2>&1; then export PATH="$K2_ZEN3:$PATH"
elif "$K2_ZEN2/kraken2-build" --version >/dev/null 2>&1; then
    echo "note: zen3 kraken2 would not run here, falling back to zen2"
    export PATH="$K2_ZEN2:$PATH"
else echo "ABORT: no runnable kraken2-build found at either pinned path."; exit 1; fi
echo "kraken2-build: $(which kraken2-build)  $(kraken2-build --version 2>&1 | head -1)"

# --name is required, and is also the inspect table's name, so read it back out of the args.
NAME=""
prev=""
for a in "$@"; do
    [ "$prev" = "--name" ] && NAME="$a"
    case "$a" in --name=*) NAME="${a#--name=}" ;; esac
    prev="$a"
done
if [ -z "$NAME" ]; then
    echo "usage: sbatch build_db.slurm --name <db> [--tier a,b] [--extra-cds F] [--build-only]"
    exit 2
fi

echo "START $(date -Is)  node $(hostname)  mem $(free -g | awk '/^Mem:/{print $2"G"}')"
SECONDS=0
python3 release/02_reference/scripts/build_db.py --threads "${THREADS:-128}" "$@"
rc=$?

# db_inspect is the coverage denominator every detection figure divides by, and it describes
# an index only this machine holds, so this machine commits it.
tools/repo_sync.sh publish \
    "02_reference: db_inspect for ${NAME} (job ${SLURM_JOB_ID:-?})" \
    "release/02_reference/data/outputs/db_inspect_${NAME}.tsv"
prc=$?
[ "$rc" -eq 0 ] && rc=$prc
echo "END $(date -Is) after ${SECONDS}s  rc=$rc"
exit $rc
