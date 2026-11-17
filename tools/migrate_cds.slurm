#!/bin/bash
#SBATCH --job-name=migrate_cds
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=tools/logs/migrate_cds_%j.out
#SBATCH --error=tools/logs/migrate_cds_%j.err
#
# One-off: move the fungal CDS tree into release/02_reference (see migrate_cds_to_release.py).
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch tools/migrate_cds.slurm
#
# PARTITION=copy. Nothing here computes: it is 5,302 renames within one Lustre filesystem
# plus a directory scan, which is metadata work and belongs off the login node. --mem=4G
# because the memory request sets the CPU bill and this needs almost none.
#
# Renames on the same filesystem do not move bytes, so the wall is generous rather than
# needed. Idempotent: an accession already at its destination is skipped, so a killed job
# resumes by resubmitting.
set -uo pipefail
module load python/3.11.6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p tools/logs
echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
python3 tools/migrate_cds_to_release.py
rc=$?
echo "exit     $rc after ${SECONDS}s"

# The manifest describes bytes only this machine holds, so this machine commits it.
tools/repo_sync.sh publish \
    "02_reference: CDS manifest after the migration into release/ (job ${SLURM_JOB_ID:-?})" \
    release/02_reference/data/raw/cds_fetch_manifest.tsv
prc=$?
[ "$rc" -eq 0 ] && rc=$prc
echo "finished $(date -Is)"
exit $rc
