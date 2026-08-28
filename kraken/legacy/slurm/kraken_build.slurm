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

# Build expanded Kraken2 pathogen DB (fungi + oomycetes, CDS-only, no BBDuk masking).
# Kraken2 LCA handles shared k-mers natively; no pre-masking required.
#
# Prerequisites:
#   1. kraken_download.slurm completed successfully
#   2. FASTAs present in $GENOMES_DIR (one subdir per accession)
#
# Run from crypt/: sbatch kraken/slurm/kraken_build.slurm

set -euo pipefail

module load python/3.11.6
module load spack/0.23.1
spack load kraken2

SCRATCH="/scratch/pawsey1168/llenzo"
DB_DIR="$SCRATCH/crypt/kraken/output/kraken_db_build/data/db_v2"
GENOMES_DIR="$SCRATCH/crypt/kraken/output/kraken_db_search/data/cds_v2"

mkdir -p kraken/slurm/logs
mkdir -p "$DB_DIR"

echo "Starting Kraken2 build: $(date)"
echo "DB dir:      $DB_DIR"
echo "Genomes dir: $GENOMES_DIR"

python kraken/build.py \
    --build-only \
    --db-dir "$DB_DIR" \
    --genomes-dir "$GENOMES_DIR" \
    --threads 32

echo "Build complete: $(date)"
