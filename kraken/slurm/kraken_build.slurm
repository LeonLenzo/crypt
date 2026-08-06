#!/bin/bash
#SBATCH --job-name=kraken_build
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=230G
#SBATCH --time=12:00:00
#SBATCH --output=kraken/slurm/logs/kraken_build_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_build_%j.err

# Run from crypt/: sbatch kraken/slurm/kraken_build.slurm
# Prerequisite: run slurm/kraken_download.slurm first (or --download-only pass)

set -euo pipefail

module load python/3.11.6
module load spack/0.23.1
spack load kraken2
module load bbmap/38.96--h5c4e2a8_0

SCRATCH="/scratch/pawsey1168/llenzo"
DB_DIR="$SCRATCH/kraken/db"
GENOMES_DIR="$SCRATCH/kraken/cds_from_genomic"
DB_JSON="output/00_build/data/phibase_db.json"

mkdir -p slurm/logs

echo "Starting Kraken2 build: $(date)"
echo "DB dir: $DB_DIR"
echo "Genomes dir: $GENOMES_DIR"

python kraken/build.py \
    --build-only \
    --db-path "$DB_JSON" \
    --db-dir "$DB_DIR" \
    --genomes-dir "$GENOMES_DIR" \
    --threads 32

echo "Build complete: $(date)"
