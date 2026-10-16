#!/bin/bash
#SBATCH --job-name=predict_compare
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_compare_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_compare_%j.err

# Score the predicted CDS against the k-mer space a classifier actually uses.
#
# --mem=16G deliberately small: on Setonix the memory request sets the CPU bill (~1.9 GB per
# core), and this is single-threaded numpy sketching, so asking for more would buy nothing and
# cost cores. See memory setonix-compute-accounting.
#
# Arm A only on this pass. The Arm B novelty measure needs sketches of db_v4's OWN CDS, and
# the cache currently holds 1,989 sketches from the db_v3 era — a SUBSET of db_v4, which would
# understate db_v4's k-mer space and so overstate novelty. Sketching all 4,271 is a separate
# job; running it before db_v4's quality is known would be spending compute on the wrong
# question first.
#
#   sbatch 02_classify/slurm/predict_compare.slurm

set -uo pipefail
SCRATCH="/scratch/pawsey1168/llenzo"
cd "$SCRATCH/crypt"
mkdir -p 02_classify/slurm/logs
module load cray-python/3.12.12
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "START $(date)"
python 02_classify/01_refs/predict/compare.py \
    --tool vipsania_pretrained \
    --predicted-dir "$SCRATCH/crypt/02_classify/data/bronze/predicted/vipsania_pretrained" \
    --cds-dir "$SCRATCH/crypt/02_classify/data/bronze/cds/pathogen" \
    --out "$SCRATCH/crypt/02_classify/data/gold/predict_benchmark.tsv"
echo "END $(date)"
