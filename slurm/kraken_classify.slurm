#!/bin/bash
#SBATCH --job-name=kraken_classify
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/logs/kraken_classify_%j.out
#SBATCH --error=slurm/logs/kraken_classify_%j.err

# Run from crypt/: sbatch slurm/kraken_classify.slurm
# Runs both MAL and HAL modes sequentially (shared cache, safe to resume).
# At ~500k reads/run with 8 workers × 4 kraken threads = 32 cores:
#   MAL (~48k runs):  ~6h
#   HAL (~560k runs): ~3 days — resubmit if job times out; cache is resumable.

set -euo pipefail

module load kraken2

DB_DIR="$MYSCRATCH/kraken_db"

mkdir -p slurm/logs

echo "Starting Kraken2 classification: $(date)"
echo "Mode: $1 (pass 'mal', 'hal', or 'both')"
echo "DB: $DB_DIR"

python kraken_classify.py \
    --mode "${1:-both}" \
    --db "$DB_DIR" \
    --workers 8 \
    --kraken-threads 4 \
    --n-reads 500000 \
    --tmp-dir "$MYSCRATCH/kraken_tmp" \
    --reads-dir "$MYSCRATCH/kraken_reads"

# Sync reads to Acacia after classification (update bucket name):
# aws s3 sync "$MYSCRATCH/kraken_reads" s3://pawsey1168-kraken/reads/ --profile acacia

echo "Classification complete: $(date)"
