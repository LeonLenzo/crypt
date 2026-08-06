#!/bin/bash
#SBATCH --job-name=kraken_classify
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-49
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_classify_%A_%a.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_classify_%A_%a.err

# Run from crypt/ on Setonix:
#   sbatch kraken/slurm/kraken_classify.slurm
#
# 50-task array — each task gets a strided 1/50 slice of all runs.
# Shared cache (output/kraken_classify/data/kraken_cache.jsonl) is safe to
# write from all tasks simultaneously via fcntl file locking.
# Safe to resubmit: completed runs are skipped on resume.
#
# Scale: ~593k runs ÷ 50 tasks ≈ 11,860 runs/task.
# At 8 workers × ~90s/run (500k reads, paired): ~37h per task — within 24h limit
# if only one mode is queued, or resubmit once. Increase --array upper bound to
# reduce per-task load (e.g. 0-99 → ~6k runs/task → ~18h).

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
CRYPT="$SCRATCH/crypt"
DB_DIR="$SCRATCH/kraken/db"

module load python/3.11.6
module load spack/0.23.1
spack load kraken2

cd "$CRYPT"
mkdir -p kraken/slurm/logs

echo "Array task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} starting: $(date)"
echo "DB: $DB_DIR"

python kraken/classify.py \
    --mode both \
    --db "$DB_DIR" \
    --workers 8 \
    --kraken-threads 4 \
    --tmp-dir "$SCRATCH/kraken_tmp/${SLURM_ARRAY_TASK_ID}" \
    --array-id "${SLURM_ARRAY_TASK_ID}" \
    --array-count "${SLURM_ARRAY_TASK_COUNT}"

echo "Array task ${SLURM_ARRAY_TASK_ID} complete: $(date)"
