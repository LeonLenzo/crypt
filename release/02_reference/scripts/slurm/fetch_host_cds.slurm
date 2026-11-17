#!/bin/bash
#SBATCH --job-name=fetch_host_cds
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=release/02_reference/scripts/slurm/logs/fetch_host_cds_%j.out
#SBATCH --error=release/02_reference/scripts/slurm/logs/fetch_host_cds_%j.err
#
# Fetch the host CDS from Ensembl Plants into release/02_reference/data/raw/host_cds/.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/02_reference/scripts/slurm/fetch_host_cds.slurm
#
# PARTITION=copy, like every fetch in this project: pure network and filesystem I/O, nothing
# decompressed and no core doing arithmetic. --mem=4G because the memory request sets the CPU
# bill (see the compute-accounting note) and a single-stream download needs almost none.
#
# 21 files, ~1.3 GB gz: 19 Triticum aestivum cultivars, Hordeum vulgare MorexV3 and
# Arabidopsis thaliana TAIR10. ONE REQUEST AT A TIME with a pause between, because EBI shares
# a 2 req/s limiter across EuropePMC and ENA and this project has already collapsed one EBI
# job by running 48 parallel streams at it. Minutes, not hours; the 4 h wall is for a bad day
# at the other end.
#
# Idempotent: a file already on disk with a non-zero size is skipped, so a killed job resumes
# by resubmitting. There is no HTTP Range resume, so a .part killed at the wall restarts from
# zero -- which for 70 MB files costs nothing worth engineering around.
set -uo pipefail
module load python/3.11.6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p release/02_reference/scripts/slurm/logs
echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
python3 release/02_reference/scripts/fetch_host_cds.py
rc=$?
echo "exit     $rc after ${SECONDS}s"

# The manifest describes bytes only this machine holds, so this machine commits it.
tools/repo_sync.sh publish \
    "02_reference: host CDS manifest (job ${SLURM_JOB_ID:-?})" \
    release/02_reference/data/raw/host_cds_manifest.tsv
prc=$?
[ "$rc" -eq 0 ] && rc=$prc
echo "finished $(date -Is)"
exit $rc
