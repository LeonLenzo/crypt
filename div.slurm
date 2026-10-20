#!/bin/bash
#SBATCH --job-name=rust_div
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/04_kallisto/slurm/logs/rust_div_%j.out
set -uo pipefail
module load cray-python/3.11.7
export OMP_NUM_THREADS=1
cd /scratch/pawsey1168/llenzo/crypt || exit 1
python3 divergence.py
