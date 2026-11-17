#!/bin/bash
#SBATCH --job-name=fetch_cds
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=release/02_reference/scripts/slurm/logs/fetch_cds_%j.out
#SBATCH --error=release/02_reference/scripts/slurm/logs/fetch_cds_%j.err
#
# Fetch the fungal CDS a build tier still lacks, into release/02_reference/data/raw/cds/.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   TIER=db1_congener sbatch release/02_reference/scripts/slurm/fetch_cds.slurm
#
# PARTITION=copy: network and filesystem I/O, nothing computed. --mem=8G sets the CPU bill
# low; `datasets` needs little.
#
# WORKERS=8 and no higher. The recurring failure against NCBI here is a dropped connection
# mid-package, not rate-limiting, and 8 parallel single-accession packages retry cleanly where
# one 700-accession package restarts whole.
#
# Idempotent: an accession whose cds_from_genomic.fna is present and non-empty is skipped, so
# a killed job resumes by resubmitting. Most of the fungal set is already on disk from the
# pre-release fetch (tools/migrate_cds_to_release.py moved it), so this closes a gap rather
# than fetching 92 GB.
set -uo pipefail
module load python/3.11.6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p release/02_reference/scripts/slurm/logs
echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
python3 release/02_reference/scripts/fetch_cds.py \
    --tier "${TIER:-db0_baseline,db1_congener}" --workers "${WORKERS:-8}"
rc=$?
echo "exit     $rc after ${SECONDS}s"

# The manifest describes bytes only this machine holds, so this machine commits it.
tools/repo_sync.sh publish \
    "02_reference: CDS fetch manifest (job ${SLURM_JOB_ID:-?})" \
    release/02_reference/data/raw/cds_fetch_manifest.tsv
prc=$?
[ "$rc" -eq 0 ] && rc=$prc
echo "finished $(date -Is)"
exit $rc
