#!/bin/bash --login
#SBATCH --job-name=align_redundancy
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
# Measure within-species transcript redundancy over every pick in a plan.
#
#   sbatch research/03_align/01_index/slurm/redundancy.slurm pleosporaceae
#
# Single-threaded on purpose: it is one pass of hashing per assembly and the whole plan is
# ~62 of them, so the 16G is for the inverted index of the largest and nothing else. The
# memory request is the CPU bill on Setonix, so it is not rounded up for comfort.
set -euo pipefail

LABEL="${1:-pleosporaceae}"
REPO="/scratch/pawsey1168/llenzo/crypt"

cd "$REPO"
echo "host $(hostname), label ${LABEL}"
time python3 research/03_align/01_index/redundancy.py --label "${LABEL}"

# The head-to-head the flag exists for: A. alternata's two Scaffold annotations, one contig
# apart in contiguity and 11,000 genes apart in content.
echo
echo "=== A. alternata, the two candidates ==="
time python3 research/03_align/01_index/redundancy.py --cds GCA_016097525.1 GCA_001642055.1
echo "done"
