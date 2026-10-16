#!/bin/bash
#SBATCH --job-name=saturation
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=230G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/saturation_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/saturation_%j.err

# The saturation QC record for db_v4: how much k-mer space each species' assemblies actually
# add, and what depth costs in congener misassignment. A RECORD, not a gate — the no-cap
# decision is already made and this is the evidence behind it. See the script docstring.
#
# Prerequisite: 02_classify/02_build/fetch_cds.py, so every manifest row has a FASTA.
#
# k=31 and scaled=200 ARE NOT FREE TO CHANGE. Two reasons:
#   * The sketch cache is parameter-specific. 1,989 sketches computed for db_v3 at these exact
#     values are hardlinked into data/bronze/sketches, so only the ~2,500 newly downloaded
#     accessions need sketching. A different k or scaled silently invalidates all of them —
#     silently because the filename is the accession only, with no parameters in it.
#   * db_v3's published curves (pangenome_ratio, n_for_95pct, depth_bias) were measured at
#     these values. Changing them makes db_v4's numbers incomparable to the ones already
#     quoted in CLAUDE.md and in the 02_build docstrings.
#
#   sbatch 02_classify/slurm/saturation.slurm

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load cray-python/3.11.7

# numpy must not thread inside each worker process — the parallelism IS the process pool, and
# letting BLAS spawn 64 threads per worker oversubscribes the node by 64x.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

GENOMES_DIR="$SCRATCH/crypt/02_classify/data/bronze/cds/pathogen"

echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-none}"
echo "genomes dir : $GENOMES_DIR"
echo "fasta on disk: $(find "$GENOMES_DIR" -name '*.fna' -size +100k | wc -l)"
echo "sketches cached: $(ls "$SCRATCH/crypt/02_classify/data/bronze/sketches" 2>/dev/null | wc -l)"
echo "python      : $(python3 -c 'import sys,numpy;print(sys.version.split()[0], "numpy", numpy.__version__)')"
echo

python3 -u 02_classify/utilities/saturation/saturation.py \
    --genomes-dir "$GENOMES_DIR" \
    --stage all \
    --workers 64 \
    --k 31 \
    --scaled 200 \
    --min-assemblies 4 \
    --permutations 25

echo
echo "done: $(date)"
