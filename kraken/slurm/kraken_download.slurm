#!/bin/bash
#SBATCH --job-name=kraken_download
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=kraken/slurm/logs/kraken_download_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_download_%j.err

# Download transcriptome/genomic FASTAs for all PHI-base euk pathogen seeds.
# Run BEFORE kraken_build.slurm (which does --build-only).
#
# Prerequisites:
#   - kraken/ref_screen.tsv exists (run: python scripts/screen_references.py)
#   - NCBI datasets CLI available (module load or in PATH)
#
# Run from crypt/: sbatch kraken/slurm/kraken_download.slurm

set -euo pipefail

module load ncbi-datasets-cli   # adjust to Setonix module name

GENOMES_DIR="$MYSCRATCH/kraken_transcriptomes"

mkdir -p slurm/logs

echo "Starting FASTA download: $(date)"
echo "Genomes dir: $GENOMES_DIR"
echo "Ref screen: kraken/ref_screen.tsv"

python kraken/build.py \
    --download-only \
    --genomes-dir "$GENOMES_DIR"

echo "Download complete: $(date)"
