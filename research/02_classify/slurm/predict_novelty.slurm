#!/bin/bash
#SBATCH --job-name=predict_novelty
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_novelty_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_novelty_%j.err

# ARM B: how much k-mer space would predicting the unannotated species actually ADD to db_v4?
#
# NO REBUILD. The question can be answered from sketches alone, because FracMinHash sketches
# of different sequences sample the SAME hash subspace: set operations on the sketches
# estimate set operations on the full k-mer sets. So the added space is
#   |predicted \ db_v4| x scaled
# computed over cached sketches, with no kraken2-build anywhere in it. Rebuilding to find out
# would cost ~6 h and ~5,000 SU; this costs minutes.
#
# Two stages:
#   1. sketch  — db_v4's own CDS. 1,749 of the 4,271 are already cached from the db_v3 era
#                (the cache is keyed by accession, so every reused accession carries over);
#                ~2,522 remain at roughly 3.4 s each.
#   2. compare — sketch each predicted CDS and measure containment against the db_v4 union.
#
# --mem=64G is sized for stage 2, not stage 1: unioning 4,271 sketches means concatenating
# ~427 M uint64 hashes (~3.4 GB) and np.unique sorts a copy on top of that. Sketching itself
# is small. Note the memory request also sets the CPU bill on Setonix, so this is deliberately
# no larger than it needs to be.
#
#   sbatch 02_classify/slurm/predict_novelty.slurm

set -uo pipefail
SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load cray-python/3.12.12
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

CDS="$SCRATCH/crypt/02_classify/data/bronze/cds/pathogen"
SKETCH="$SCRATCH/crypt/02_classify/data/bronze/sketches"
PRED="$SCRATCH/crypt/02_classify/data/bronze/predicted/vipsania_pretrained"

echo "START $(date)"
echo "sketches cached before: $(ls "$SKETCH"/*.npy 2>/dev/null | wc -l)"

echo
echo "──── stage 1: sketch db_v4's CDS ────"
# k=31 scaled=200 MUST match the cached sketches and the Arm A containments; the cache
# filename is the accession alone and carries no parameters, so a different value silently
# mixes incompatible sketches.
python3 -u 02_classify/utilities/saturation/saturation.py \
    --genomes-dir "$CDS" --stage sketch --workers 32 --k 31 --scaled 200 2>&1 | tail -20

echo "sketches cached after: $(ls "$SKETCH"/*.npy 2>/dev/null | wc -l)"

echo
echo "──── stage 2: novelty of the predicted CDS against db_v4 ────"
python3 -u 02_classify/01_refs/predict/compare.py \
    --tool vipsania_pretrained \
    --predicted-dir "$PRED" \
    --cds-dir "$CDS" \
    --db-sketches "$SKETCH" \
    --out "$SCRATCH/crypt/02_classify/data/gold/predict_benchmark.tsv" 2>&1 | tail -30

echo "END $(date)"
