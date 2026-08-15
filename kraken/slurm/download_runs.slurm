#!/bin/bash
#SBATCH --job-name=download_runs
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --array=0-4
#SBATCH --output=kraken/slurm/logs/download_runs_%A_%a.out
#SBATCH --error=kraken/slurm/logs/download_runs_%A_%a.err

# Download complete FASTQ.gz files for all biosample-representative runs that
# reach the Sankey endpoint (~4,200 runs, ~10 TB total).
#
# Runs as a 5-task array so each task lands on a separate compute node with a
# distinct IP — this multiplies ENA's per-IP concurrency limit (empirically ~12
# workers/IP) by the array size without changing worker count per task.
#
# Each task processes a non-overlapping slice of the run list (0-4 of 5).
# Resumable: already-downloaded files are skipped via download_index.txt.
#
# Prerequisites:
#   stat/output/data/runs.tsv must exist (stat/filter_runs.py)
#
# Run from crypt/: sbatch kraken/slurm/download_runs.slurm
# Monitor:         tail -f kraken/slurm/logs/download_runs_<JOBID>_*.out
# After complete:  python kraken/host_split.py  (TBD)

set -euo pipefail

module load python/3.11.6

SCRATCH="/scratch/pawsey1168/llenzo"
READS_DIR="$SCRATCH/kraken/reads_full"

mkdir -p "$READS_DIR"
mkdir -p kraken/slurm/logs

echo "Array task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}  node: $(hostname)"
echo "Reads dir: $READS_DIR"
echo "Started:   $(date)"

python kraken/download.py \
    --runs-tsv  stat/output/data/runs.tsv \
    --biosample-rep \
    --reads-dir "$READS_DIR" \
    --full \
    --workers   12 \
    --array-id    "$SLURM_ARRAY_TASK_ID" \
    --array-count "$SLURM_ARRAY_TASK_COUNT"

echo "Finished: $(date)"
