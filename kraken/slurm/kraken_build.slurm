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

# Pathogen-only Kraken2 DB build (fungi + oomycetes, no host sequences).
# BBDuk masking: pathogen-vs-pathogen shared k-mers only.
#
# Prerequisites:
#   1. FASTAs on scratch (from kraken_download.slurm or previous build)
#   2. Stale BBDuk intermediates removed (done below — they used old bidirectional masking)
#
# Run from crypt/: sbatch kraken/slurm/kraken_build.slurm

set -euo pipefail

module load python/3.11.6
module load spack/0.23.1
spack load kraken2
module load bbmap/38.96--h5c4e2a8_0

SCRATCH="/scratch/pawsey1168/llenzo"
DB_DIR="$SCRATCH/kraken/db_pathogens"
GENOMES_DIR="$SCRATCH/kraken/cds_from_genomic"
DB_JSON="output/00_build/data/phibase_db.json"

mkdir -p slurm/logs
mkdir -p "$DB_DIR"

echo "Starting pathogen-only Kraken2 build: $(date)"
echo "DB dir:      $DB_DIR"
echo "Genomes dir: $GENOMES_DIR"

# Remove stale BBDuk intermediates — the old masked files used bidirectional
# (pathogen + host) masking; the new build uses pathogen-vs-pathogen only.
echo "Removing stale BBDuk intermediates …"
rm -f "$GENOMES_DIR"/_pathogen_combined.fna \
      "$GENOMES_DIR"/_pathogen_masked.fna \
      "$GENOMES_DIR"/_shared_pathogen_kmers.fna \
      "$GENOMES_DIR"/_host_combined.fna \
      "$GENOMES_DIR"/_host_masked.fna \
      "$GENOMES_DIR"/_host_tmp.fna \
      "$GENOMES_DIR"/_shared_host_kmers.fna
echo "Done."

python kraken/build.py \
    --build-only \
    --db-path "$DB_JSON" \
    --db-dir "$DB_DIR" \
    --genomes-dir "$GENOMES_DIR" \
    --threads 32 \
    --bbduk-mem 200g

echo "Build complete: $(date)"
