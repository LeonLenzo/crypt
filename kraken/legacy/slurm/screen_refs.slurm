#!/bin/bash
#SBATCH --job-name=screen_refs
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=kraken/slurm/logs/screen_refs_%j.out
#SBATCH --error=kraken/slurm/logs/screen_refs_%j.err

# Query NCBI for all candidate assemblies (no N cap).
# Produces kraken/ref_candidates.tsv for busco_screen.slurm.
# Requires NCBI datasets CLI (bootstrapped below).
#
# Run from crypt/: sbatch kraken/slurm/screen_refs.slurm
# Then:           sbatch kraken/slurm/busco_screen.slurm

set -euo pipefail

module load python/3.11.6

SCRATCH="/scratch/pawsey1168/llenzo"
DATASETS_BIN="$SCRATCH/bin/datasets"

if [ ! -x "$DATASETS_BIN" ]; then
    echo "Downloading NCBI datasets CLI …"
    mkdir -p "$SCRATCH/bin"
    curl -fsSL -o "$DATASETS_BIN" \
        https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$DATASETS_BIN"
fi
export PATH="$SCRATCH/bin:$PATH"
echo "datasets version: $(datasets version)"

mkdir -p kraken/slurm/logs

echo "Starting assembly screen: $(date)"

python kraken/screen_refs.py --workers 8

echo "Screen complete: $(date)"
echo "Candidates: $(wc -l < kraken/ref_candidates.tsv) rows (incl. header)"
echo "Next: sbatch kraken/slurm/busco_screen.slurm"
