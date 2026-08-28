#!/bin/bash
#SBATCH --job-name=kraken_classify_path
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_classify_pathogens_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_classify_pathogens_%j.err

# Classify 50 HC biosample-rep runs against the pathogen-only DB (db_pathogens).
# FASTQs from the previous masked/unmasked test runs are reused from kraken_tmp.
#
# Run from crypt/ on Setonix:
#   sbatch kraken/slurm/kraken_classify_pathogens_test.slurm

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
CRYPT="$SCRATCH/crypt"
DB_DIR="$SCRATCH/crypt/kraken/output/kraken_db_build/data/db_pathogens"

module load python/3.11.6
module load spack/0.23.1
spack load kraken2

cd "$CRYPT"
mkdir -p kraken/slurm/logs

echo "Starting pathogen-only DB classification (50 runs): $(date)"
echo "DB: $DB_DIR"

python kraken/classify.py \
    --runs-tsv output/02_filter_runs/data/runs.tsv \
    --biosample-rep \
    --hc \
    --db "$DB_DIR" \
    --out-dir output/kraken_classify_pathogens \
    --limit 50 \
    --workers 8 \
    --kraken-threads 4 \
    --tmp-dir "$SCRATCH/kraken_tmp"

echo "Classification complete: $(date)"
