#!/bin/bash
#SBATCH --job-name=kraken_classify_test
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_classify_test_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_classify_test_%j.err

# Run from crypt/ on Setonix:
#   sbatch kraken/slurm/kraken_classify_test.slurm
# Classifies 50 HC biosample-representative runs to validate the masked DB.
# Prerequisite: scp output/02_filter_runs/data/runs.tsv to Setonix first.

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
CRYPT="$SCRATCH/crypt"
DB_DIR="$SCRATCH/kraken/db"

module load python/3.11.6
module load spack/0.23.1
spack load kraken2

cd "$CRYPT"
mkdir -p kraken/slurm/logs

echo "Starting test classification (50 runs): $(date)"
echo "DB: $DB_DIR"

python kraken/classify.py \
    --runs-tsv output/02_filter_runs/data/runs.tsv \
    --biosample-rep \
    --hc \
    --db "$DB_DIR" \
    --limit 50 \
    --workers 8 \
    --kraken-threads 4 \
    --tmp-dir "$SCRATCH/kraken_tmp"

echo "Test classification complete: $(date)"
