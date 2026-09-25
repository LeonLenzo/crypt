#!/bin/bash
#SBATCH --job-name=kraken_db_build
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=230G
#SBATCH --time=12:00:00
#SBATCH --output=kraken/slurm/logs/kraken_db_build_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_db_build_%j.err

# Build the Kraken2 pathogen DB (fungi + oomycetes, CDS-only, no BBDuk masking)
# from BUSCO-selected assemblies. Kraken2 LCA handles shared k-mers natively;
# no pre-masking required.
#
# Prerequisites:
#   1. sbatch kraken/slurm/kraken_db_search.slurm  →  ref_candidates.tsv + cds/pathogen/
#   2. sbatch kraken/slurm/kraken_db_hosts.slurm   →  host_candidates.tsv + cds/host/
#   3. sbatch kraken/slurm/kraken_db_busco.slurm   →  busco_scores.tsv (selected=True rows)
# db_v3 adds host plant CDS (--host-table); db_v2 was fungi + oomycetes only.
#
# Run from crypt/: sbatch kraken/slurm/kraken_db_build.slurm

set -euo pipefail

module load python/3.11.6
module load spack/0.23.1
spack load kraken2

SCRATCH="/scratch/pawsey1168/llenzo"
DB_DIR="$SCRATCH/crypt/kraken/build/data/db_v3"
GENOMES_DIR="$SCRATCH/crypt/kraken/search/data/cds/pathogen"

mkdir -p kraken/slurm/logs
mkdir -p "$DB_DIR"

echo "Starting Kraken2 build: $(date)"
echo "DB dir:      $DB_DIR"
echo "Genomes dir: $GENOMES_DIR"

HOST_TABLE="kraken/search/data/host_candidates.tsv"
HOST_GENOMES_DIR="$SCRATCH/crypt/kraken/search/data/cds/host"

python kraken/build/kraken_build.py \
    --db-dir "$DB_DIR" \
    --genomes-dir "$GENOMES_DIR" \
    --host-table "$HOST_TABLE" \
    --host-genomes-dir "$HOST_GENOMES_DIR" \
    --threads 32

echo "Build complete: $(date)"
