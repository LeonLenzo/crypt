#!/bin/bash
#SBATCH --job-name=kraken_build
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/kraken_build_%j.out
#SBATCH --error=slurm/logs/kraken_build_%j.err

# Run from crypt/: sbatch slurm/kraken_build.slurm
# Prerequisite: run slurm/kraken_download.slurm first (or --download-only pass)

set -euo pipefail

module load kraken2
module load ncbi-datasets-cli   # adjust module name to match Setonix modules

DB_DIR="$MYSCRATCH/kraken_db"
GENOMES_DIR="$MYSCRATCH/kraken_transcriptomes"

mkdir -p slurm/logs

echo "Starting Kraken2 build: $(date)"
echo "DB dir: $DB_DIR"
echo "Genomes dir: $GENOMES_DIR"

python kraken_build.py \
    --build-only \
    --db-dir "$DB_DIR" \
    --genomes-dir "$GENOMES_DIR" \
    --threads 32

echo "Build complete: $(date)"
