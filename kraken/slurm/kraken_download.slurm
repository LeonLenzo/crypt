#!/bin/bash
#SBATCH --job-name=kraken_download
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=kraken/slurm/logs/kraken_download_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_download_%j.err

# Pull pre-downloaded transcriptome FASTAs from Acacia into scratch.
# FASTAs were downloaded locally and uploaded via:
#   python kraken/build.py --download-only --upload-to-acacia
#
# Run BEFORE kraken_build.slurm (which does --build-only).
# Run from crypt/: sbatch kraken/slurm/kraken_download.slurm

set -euo pipefail

GENOMES_DIR="$MYSCRATCH/cds_from_genomic"
BUCKET="pawsey1168-llenzo-kraken-db"
ENDPOINT="https://projects.pawsey.org.au"

mkdir -p "$GENOMES_DIR"
mkdir -p kraken/slurm/logs

echo "Pulling FASTAs from Acacia: $(date)"
echo "Source: s3://${BUCKET}/cds_from_genomic/"
echo "Dest:   $GENOMES_DIR"

aws s3 sync "s3://${BUCKET}/cds_from_genomic/" "$GENOMES_DIR/" \
    --profile acacia \
    --endpoint-url "$ENDPOINT"

echo "Download complete: $(date)"
echo "Contents:"
ls "$GENOMES_DIR" | wc -l
