#!/bin/bash
#SBATCH --job-name=krsel_full
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=kraken/slurm/logs/krsel_full_%j.out
#SBATCH --error=kraken/slurm/logs/krsel_full_%j.err

# Full-cohort read download for submodule 2 step 1/3. Replaces the ad-hoc
# `sbatch --wrap=...` used for the 2026-09-08 chain (jobs 48206159-48206172),
# which left no script in the repo and wrote its logs to crypt/ root.
#
# The `work` partition caps a single job at 24h and the full pull needs several
# days, so this is designed to be submitted as an --dependency=afterany chain:
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   prev=""
#   for i in $(seq 1 6); do
#     if [ -z "$prev" ]; then
#       jid=$(sbatch --parsable kraken/slurm/kraken_run_select.slurm)
#     else
#       jid=$(sbatch --parsable --dependency=afterany:$prev kraken/slurm/kraken_run_select.slurm)
#     fi
#     echo "submitted $jid"; prev=$jid
#   done
#
# afterany (not afterok) is correct: each job is EXPECTED to hit the 24h wall
# and be recorded TIMEOUT. That is not failure. kraken_run_select.py is
# resumable, so the next job picks up where the wall cut the last one off.
#
# --workers 8 is deliberate and must not be raised. --workers 16 was tested and
# is WORSE: NCBI connection contention, no throughput gain.
#
# KRSEL_GZIP_THREADS sizes pigz inside _compress(). Compression is a rare,
# bursty event relative to the network-bound downloads, so oversubscribing
# briefly is harmless and finishing each compression fast matters more.

set -uo pipefail

module load python/3.11.6

export KRSEL_GZIP_THREADS=8

cd /scratch/pawsey1168/llenzo/crypt || exit 1

echo "host=$(hostname) job=$SLURM_JOB_ID cpus=${SLURM_CPUS_PER_TASK:-?} pigz_threads=$KRSEL_GZIP_THREADS"
python kraken/select/kraken_select.py --download --workers 8
