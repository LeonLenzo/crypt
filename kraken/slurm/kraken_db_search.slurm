#!/bin/bash
#SBATCH --job-name=kraken_db_search
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=kraken/slurm/logs/kraken_db_search_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_db_search_%j.err

# Select candidate assemblies (query NCBI, no N cap) AND download their CDS FASTA.
# This is the ONE place in the kraken_db_* pipeline that downloads — kraken_db_busco.py
# and kraken_db_build.py both read from data/cds/pathogen/ produced here.
# Requires NCBI datasets CLI (bootstrapped below).
#
# Run from crypt/: sbatch kraken/slurm/kraken_db_search.slurm
# Then:           sbatch kraken/slurm/kraken_db_busco.slurm

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

echo "Starting assembly search: $(date)"

python kraken/kraken_db_search.py --download --workers 16

echo "Search complete: $(date)"
echo "Candidates: $(wc -l < kraken/output/kraken_db_search/data/ref_candidates.tsv) rows (incl. header)"
echo "Next: sbatch kraken/slurm/kraken_db_busco.slurm"
