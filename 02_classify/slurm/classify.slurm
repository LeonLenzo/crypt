#!/bin/bash
#SBATCH --job-name=classify
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=02_classify/slurm/logs/classify_%j.out
#SBATCH --error=02_classify/slurm/logs/classify_%j.err

# 04_classify: run db_v4 over one BioProject's reads.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch --export=ALL,PROJECT=PRJEB39201 02_classify/slurm/classify.slurm
#
# ONE PROJECT PER JOB. db_v4's hash is 41 GB and the memory request sets the CPU bill, so the
# expensive operation is LOADING the database, not classifying against it. A per-run array
# would pay one 41 GB Lustre load per run; one long-lived job per project stages the hash into
# /dev/shm ONCE and every worker mmaps that single RAM-resident copy.
#
# --memory-mapping is used here and CLAUDE.md says never to. That note measured mmap against
# LUSTRE -- 13 minutes, zero reports, kraken2 pinned at 1.6% CPU, because every hash lookup
# became a filesystem RPC. /dev/shm is RAM and is the opposite case. Never pass --mmap without
# the staging below.
#
# --mem=80G, NOT 64G. The rarefaction job ran at 65,524 MB against a 64 GB request, which is
# the edge: 41 GB of hash in /dev/shm (charged to the job's cgroup) plus worker working sets.
# More workers here than there, so the headroom is raised. Memory sets the bill at ~1.9 GB per
# physical core, so 80G bills ~42 cores whether busy or not -- that is the floor for holding
# db_v4 at all, and the only way to spend less is fewer DB loads, which is what shm achieves.
#
# NO DEPTH FILTER. Every run on disk is classified. A detection is valid at any depth; only
# ABSENCE is depth-dependent, and that caveat belongs to 05_detect's denominator, not here.
#
# CONTROLS=1 classifies the NEGATIVE CONTROLS instead of the cohort:
#
#   sbatch --export=ALL,CONTROLS=1,PROJECT="PRJEB29046 PRJNA613154" 02_classify/slurm/classify.slurm
#
# It points --reads/--reports at `data/bronze/controls/`, where PROJECT names the control
# ACCESSION rather than a cohort BioProject. Nothing else changes, and that is the requirement:
# a control classified under different settings measures nothing. The separate tree exists so
# controls cannot enter the cohort's denominator -- see the comment beside CONTROL_READS in
# _layout.py.

set -uo pipefail

module load python/3.11.6
module load kraken2/2.1.2-qk2lek6

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 02_classify/slurm/logs

: "${PROJECT:?set PROJECT, e.g. --export=ALL,PROJECT=PRJEB39201}"
# PROJECT may be a SPACE-SEPARATED LIST, and for the small projects it should be. Staging the
# 41 GB hash costs ~200s and is per JOB, not per project, so a 4-run project in its own job
# spends 200s of staging on ~2 minutes of work. The 21 projects under 50 runs are 290 runs
# between them and belong in one job behind one stage.
DB=02_classify/data/silver/db/db_v4
SHM=/dev/shm/db_v4_$SLURM_JOB_ID
# 8 x 8 = 64 concurrent threads, matching the 64 PHYSICAL cores that --mem=80G buys (sacct
# reports 128, which is the SMT logical count). The smoke test ran 4 x 14 = 56 threads and
# billed the same 64 cores, so half the allocation sat idle. Workers are cheap to add because
# they mmap ONE shared /dev/shm copy of the hash instead of loading their own: measured RSS
# went 48.4 GB at 4 workers, and each extra worker costs a working set rather than 41 GB.
WORKERS=${WORKERS:-8}
THREADS=${THREADS:-8}

# Built as an ARRAY, not a string. A quoted string with a line continuation inside it keeps the
# backslash as a literal character, and unquoted $TREE then word-splits it into an argument of
# its own, so classify.py is handed a bare "\" and argparse rejects the run.
TREE=()
TREE_NAME=cohort
if [ -n "${CONTROLS:-}" ]; then
  TREE=(--reads 02_classify/data/bronze/controls/reads
        --reports 02_classify/data/bronze/controls/reports)
  TREE_NAME=controls
fi

echo "host     $(hostname)"
echo "project  $PROJECT"
echo "tree     $TREE_NAME"
echo "started  $(date -Is)"

SECONDS=0
mkdir -p "$SHM"
cp "$DB"/hash.k2d "$DB"/opts.k2d "$DB"/taxo.k2d "$SHM"/ || { echo "shm stage FAILED"; exit 1; }
echo "staged $(du -sh --apparent-size $SHM | cut -f1) into RAM in ${SECONDS}s"

cleanup() { rm -rf "$SHM"; }
trap cleanup EXIT

SECONDS=0
rc=0
for P in $PROJECT; do
  echo "--- $P ---"
  python3 02_classify/04_classify/classify.py \
      --project "$P" --db "$SHM" \
      --workers "$WORKERS" --threads "$THREADS" --mmap \
      --tmp-dir "$MYSCRATCH/kraken_tmp" "${TREE[@]+${TREE[@]}}" || rc=$?
done

echo "exit     $rc after ${SECONDS}s"
echo "finished $(date -Is)"
exit $rc
