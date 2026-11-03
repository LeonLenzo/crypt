#!/bin/bash
#SBATCH --job-name=predict_fetch
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_fetch_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_fetch_%j.err

# Genomic FASTA for all 119 benchmark assemblies, ~7 GB. `copy` partition: this moves bytes
# and does no computation. db_v4 holds CDS only, so there is nothing to reconcile against
# Acacia here — these accessions were never downloaded in genomic form.
#
#   sbatch 02_classify/slurm/predict_fetch.slurm

set -uo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load cray-python/3.12.12

DATASETS_BIN="$SCRATCH/bin/datasets"
if [ ! -x "$DATASETS_BIN" ]; then
    mkdir -p "$SCRATCH/bin"
    curl -sSL -o "$DATASETS_BIN" \
        https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$DATASETS_BIN"
fi
export PATH="$SCRATCH/bin:$PATH"
echo "datasets: $(datasets version)"

echo "START $(date)"
python 02_classify/01_refs/predict/fetch_genomes.py \
    --genomes-dir "$SCRATCH/crypt/02_classify/data/bronze/genomes" \
    --workers 8
echo "END $(date)"
