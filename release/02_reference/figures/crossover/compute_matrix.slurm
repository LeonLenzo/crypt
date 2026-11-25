#!/bin/bash
#SBATCH --job-name=wp_matrix
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/_analysis/wpjob/wp_matrix_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/_analysis/wpjob/wp_matrix_%j.out
set -uo pipefail
cd /scratch/pawsey1168/llenzo/crypt/_analysis/wpjob || exit 1
module load py-numpy/2.1.2 2>/dev/null
python3 wp_matrix.py wheatrow_meta.tsv wheatrow_matrix.tsv
echo "RC=$?"
