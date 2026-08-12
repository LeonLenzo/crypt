#!/bin/bash
#SBATCH --job-name=kraken_download
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=kraken/slurm/logs/kraken_download_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_download_%j.err

# Download 605 CDS assemblies from NCBI for the expanded Kraken2 pathogen DB.
# Selected by kraken/screen_refs.py (pan-genome + genus fill-in, CDS-only).
# Downloads directly from NCBI datasets CLI; no Acacia transfer needed.
#
# Run BEFORE kraken_build.slurm (which does --build-only).
# Run from crypt/: sbatch kraken/slurm/kraken_download.slurm

set -euo pipefail

module load python/3.11.6
module load spack/0.23.1
spack load ncbi-datasets-cli

SCRATCH="/scratch/pawsey1168/llenzo"
GENOMES_DIR="$SCRATCH/kraken/cds_v2"

mkdir -p "$GENOMES_DIR"
mkdir -p kraken/slurm/logs

echo "Starting CDS download: $(date)"
echo "Assemblies: $(tail -n +2 kraken/ref_screen.tsv | grep -c 'cds') CDS entries in ref_screen.tsv"
echo "Dest: $GENOMES_DIR"

python kraken/build.py \
    --download-only \
    --genomes-dir "$GENOMES_DIR" \
    --ref-screen kraken/ref_screen.tsv

echo "Download complete: $(date)"
echo "Subdirs written: $(ls "$GENOMES_DIR" | wc -l)"
