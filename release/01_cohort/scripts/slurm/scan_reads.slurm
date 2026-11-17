#!/bin/bash
#SBATCH --job-name=scan_reads
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=release/01_cohort/scripts/slurm/logs/scan_reads_%j.out
#SBATCH --error=release/01_cohort/scripts/slurm/logs/scan_reads_%j.err

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
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p release/01_cohort/scripts/slurm/logs
module load python/3.11.6
tools/repo_sync.sh require-current || exit 1
echo "START $(date -Is)"
python3 -u release/01_cohort/scripts/fetch_reads.py --reads-dir release/01_cohort/raw/reads
rc=$?
tools/repo_sync.sh publish \
    "01_cohort: reads manifest from a disk scan (job ${SLURM_JOB_ID:-?})" \
    release/01_cohort/data/raw/reads_manifest.tsv
prc=$?
[ "$rc" -eq 0 ] && rc=$prc
echo "END $(date -Is)  rc=$rc"
exit $rc
