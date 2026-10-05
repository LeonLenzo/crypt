#!/bin/bash
#SBATCH --job-name=kallisto_nb
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=58G
#SBATCH --time=02:00:00
#SBATCH --output=04_kallisto/slurm/logs/kallisto_nb_%j.out
#SBATCH --error=04_kallisto/slurm/logs/kallisto_nb_%j.err

# Measure k-mer containment between every flagged species and every other species in the
# reference set, to decide index membership in 02_build. Reuses the saturation module's cached
# FracMinHash sketches, so it downloads and sketches nothing and runs in minutes.
#
# Needs numpy, hence cray-python and NOT a login node.
#
# Run from crypt/: sbatch 04_kallisto/slurm/kallisto_neighbours.slurm

set -uo pipefail

module load cray-python/3.11.7

# numpy must not thread inside the intersections; the work is already one pass per pair
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 04_kallisto/slurm/logs

echo "host : $(hostname)"
echo "job  : ${SLURM_JOB_ID:-none}"
echo

python3 04_kallisto/02_build/kallisto_neighbours.py
echo
echo "exit $?"
