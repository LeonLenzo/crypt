#!/bin/bash
#SBATCH --job-name=scan_reads
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/release/01_cohort/slurm/logs/scan_reads_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/release/01_cohort/slurm/logs/scan_reads_%j.err

# Rewrite reads_manifest.tsv from a DISK SCAN, fetching nothing.
#
# Needed because the manifest used to be written only after the fetch pool drained, so every
# job killed at its 23 h wall left the manifest describing the previous run: on 2026-11-17 it
# was dated 2026-11-13 and still showed the 48-worker collapse while 875 runs sat on disk.
# fetch_reads.py now writes before, during and after a fetch, but a scan-only pass is still
# the way to get a current manifest WITHOUT waiting on a fetch, and it is safe to run while a
# fetch job is in flight: presence is read off the filesystem and nothing is written but the
# manifest, through a temp-file rename.
set -uo pipefail
cd /scratch/pawsey1168/llenzo/crypt/release/01_cohort
module load python/3.11.6
echo "START $(date -Is)"
python -u scripts/fetch_reads.py --reads-dir raw/reads
echo "END $(date -Is)  rc=$?"
