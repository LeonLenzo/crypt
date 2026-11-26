#!/bin/bash --login
#SBATCH --job-name=align_index
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
# Build the bowtie2 index over the tagged reference.
#
#   sbatch research/03_align/01_index/slurm/build_index.slurm pleosporaceae
#
# --export=NONE is Setonix's default, so this script sets up its own environment and never
# inherits one (crypt-two-machines-one-tree). The memory request is the CPU bill on Setonix,
# so it is sized to bowtie2-build over ~1.1 Gbp and no further; a wider scope needs more.
#
# THE INDEX LIVES ON SCRATCH, under the 21-day purge. It is a derived artefact and rebuilding
# it is one job, so it is not kept alive; the reference FASTA beside it is what makes the
# rebuild reproducible, and the plan in the repo is what makes the FASTA reproducible.
set -euo pipefail

LABEL="${1:-pleosporaceae}"
REPO="/scratch/pawsey1168/llenzo/crypt"
WORK="${ALIGN_WORK:-/scratch/pawsey1168/llenzo/_align}/${LABEL}"

module load bowtie2/2.5.4--he96a11b_7

cd "$WORK"
test -s reference.fa || { echo "no reference.fa in $WORK; run build_reference.py first" >&2; exit 1; }

echo "host $(hostname), ${SLURM_CPUS_PER_TASK} threads, label ${LABEL}"
echo "reference: $(grep -c '^>' reference.fa) records, $(du -h reference.fa | cut -f1)"

time bowtie2-build --threads "${SLURM_CPUS_PER_TASK}" reference.fa index

ls -la index*.bt2* | tee "${WORK}/index_files.txt"
echo "done"
