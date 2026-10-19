#!/bin/bash
#SBATCH --job-name=migreads
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=02_classify/slurm/logs/migreads_%j.out
#SBATCH --error=02_classify/slurm/logs/migreads_%j.err

# 03_cohort: move the flat 03_kraken reads tree into reads/<BioProject>/.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch 02_classify/slurm/migrate_reads.slurm
#
# This is a JOB and not a login-node command even though it computes nothing, because 2,362
# renames is real work on a SHARED Lustre metadata server and that is exactly the kind of thing
# the no-compute-on-login rule is protecting. It is pure metadata: no bytes are read or
# written, so 4.05 TiB relocates in seconds.
#
# --mem=4G, NOT because the script needs 4 GB (it needs ~50 MB) but because the MEMORY request
# sets the CPU bill on Setonix at ~1.9 GB/core, so this asks for ~2 cores' worth and nothing
# more. A --mem=120G habit billed 134 cores for a job using 30 GB once already.
#
# Idempotent: every move is an independent rename, and a destination that already exists aborts
# rather than overwriting. Rerun after any interruption and it continues from where it stopped.

set -uo pipefail

module load python/3.11.6

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 02_classify/slurm/logs

echo "host      $(hostname)"
echo "started   $(date -Is)"
SECONDS=0

python3 02_classify/03_cohort/migrate_reads.py --apply --clean-orphans
rc=$?

echo "exit      $rc"
echo "elapsed   ${SECONDS}s"
echo "finished  $(date -Is)"
exit $rc
