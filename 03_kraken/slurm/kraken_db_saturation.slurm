#!/bin/bash
#SBATCH --job-name=kraken_db_saturation
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=02:00:00
#SBATCH --output=03_kraken/slurm/logs/kraken_db_saturation_%j.out
#SBATCH --error=03_kraken/slurm/logs/kraken_db_saturation_%j.err

# Measure how many assemblies per species the Kraken2 build actually needs, by
# FracMinHash k-mer accumulation over the CDS already on disk. Never downloads.
#
# Two stages, both resumable. Sketching is the expensive one (one .npy per
# accession, cached); the curves rebuild from cached sketches in seconds, so
# re-running with different --min-assemblies or --permutations is cheap.
#
# Prerequisites:
#   1. sbatch 03_kraken/slurm/kraken_db_search.slurm  →  cds/pathogen/
#   2. sbatch 03_kraken/slurm/kraken_db_busco.slurm   →  busco_scores.tsv
#
# Run from crypt/: sbatch 03_kraken/slurm/kraken_db_saturation.slurm

set -euo pipefail

module load cray-python/3.11.7

# numpy must not thread inside each worker process — the parallelism is the pool
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

SCRATCH="/scratch/pawsey1168/llenzo"
GENOMES_DIR="$SCRATCH/crypt/03_kraken/01_search/data/cds/pathogen"

mkdir -p 03_kraken/slurm/logs

echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-none}"
echo "genomes dir : $GENOMES_DIR"
echo "python      : $(python3 -c 'import sys,numpy;print(sys.version.split()[0], "numpy", numpy.__version__)')"
echo

python3 -u 03_kraken/utilities/saturation/kraken_saturation.py \
    --genomes-dir "$GENOMES_DIR" \
    --stage all \
    --workers 64 \
    --k 31 \
    --scaled 200 \
    --min-assemblies 4 \
    --permutations 25

echo
echo "done: $(date)"
