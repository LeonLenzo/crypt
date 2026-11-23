#!/bin/bash
#SBATCH --job-name=congener
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --output=release/03_classify/scripts/slurm/logs/congener_%j.out
#SBATCH --error=release/03_classify/scripts/slurm/logs/congener_%j.err
#
# Pairwise minimizer sharing within every genus the cohort calls.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/03_classify/scripts/slurm/congener_sharing.slurm
#
# Sized from the crossover run it reuses (wp_matrix 53387299: 254 taxa, 64 cores, 200G, under
# 2 h). This reads more assemblies, 1,342 against that run's set, but computes HALF the pairs
# because only within-genus comparisons are made, 15,526 against 32,131. 6 h is headroom on
# the assembly-reading phase, which is the part that scales with the bigger taxon set.
#
# 200G because phase 2 forks the whole minimizer table to every worker; the sets total a few
# billion 64-bit minimizers across 396 taxa. The memory request is the CPU bill on Setonix,
# so this is sized to that fork and not above it.

set -uo pipefail
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p release/03_classify/scripts/slurm/logs release/03_classify/data/raw
module load py-numpy/2.1.2 2>/dev/null

echo "host    $(hostname)"
echo "started $(date -Is)"
echo "commit  $(git rev-parse --short HEAD)"

python3 release/03_classify/scripts/congener_sharing.py \
        release/03_classify/data/inputs/congener_meta.tsv \
        release/03_classify/data/raw/congener_sharing.tsv
rc=$?
echo "finished $(date -Is)  rc=$rc"
exit $rc
