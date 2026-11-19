#!/bin/bash
#SBATCH --job-name=refetch_failed
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/release/01_cohort/scripts/slurm/logs/refetch_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/release/01_cohort/scripts/slurm/logs/refetch_%j.err
#
# Re-fetch the runs job 53279097 could not land, at LOW concurrency.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/01_cohort/scripts/slurm/refetch_failed.slurm
#
# WHY 4 WORKERS AND NOT 24. 53279097 finished 1,533 of 1,560 runs and left 27, every one of
# them a transport symptom: 21 md5 mismatch, 3 remote-end-closed, 2 HTTP 403, 1 connection
# refused. None is a data problem. Seven were re-requested by hand on 2026-11-18 (job
# 53347244) one at a time: three came back byte-exact and checksum-exact, four got HTTP 403
# after 77 minutes of sustained transfer, and HEAD-testing every URL afterwards returned 200
# on all of them. So the files are intact at EBI and the failures are EBI throttling Pawsey's
# egress under sustained load. 24 streams is what reaches that threshold fastest.
#
# This is ~70 GB over 49 files. At the ~1 MB/s a single throttled stream gave, four streams
# is a few hours against the 20 h wall, with room for the backoff retries to do their work.
#
# NOTHING COMPLETE IS RE-REQUESTED. fetch_reads.py checks the disk, not the manifest, so this
# touches only the 27 incomplete runs; 22 of them already hold one mate of two and are queued
# first. Re-running it after a partial result is safe and resumes.
#
# If 4 workers still trips the throttle, the next step is NCBI rather than EBI: a different
# host and so a different limit. It is wrong at cohort scale (~10,000 SU for 14 TB, measured)
# and perfectly reasonable for 27 runs.
set -uo pipefail
module load python/3.11.6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p release/01_cohort/scripts/slurm/logs
tools/repo_sync.sh require-current || exit 1
echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
python3 -u release/01_cohort/scripts/fetch_reads.py \
    --reads-dir release/01_cohort/raw/reads --fetch --workers "${WORKERS:-4}"
rc=$?
echo "python   rc=$rc after ${SECONDS}s"

# The manifest describes bytes only this machine holds, so this machine commits it.
tools/repo_sync.sh publish \
    "01_cohort: reads manifest after refetch job ${SLURM_JOB_ID:-?}" \
    release/01_cohort/data/raw/reads_manifest.tsv

# fetch_reads.py counts failures and still exits 0, so the job's own exit code has to carry
# them or a half-finished fetch reports success. Same guard as fetch_wheat_reads.slurm.
failed=$(awk -F'\t' 'NR>1 && $NF ~ /^failed/ {n++} END {print n+0}' \
    release/01_cohort/data/raw/reads_manifest.tsv)
echo "failed_runs=$failed"
echo "finished $(date -Is)"
[ "$rc" -eq 0 ] && [ "$failed" -gt 0 ] && { echo "FAILING the job: rc=0 failed_runs=$failed"; exit 1; }
exit $rc
