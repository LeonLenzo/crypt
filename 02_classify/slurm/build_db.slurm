#!/bin/bash
#SBATCH --job-name=build_db_v4
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/build_db_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/build_db_%j.err
#SBATCH --dependency=singleton

# Build db_v4 from the CDS that fetch_cds.py placed on scratch. Fetches nothing.
#
# SIZING, measured not estimated. db_v3: 2,081 FASTAs / ~25.5 GB of CDS -> 4,800,752,701
# distinct minimizers -> hash.k2d 27.4 GB at table capacity 6,856,360,594 (4 bytes per cell,
# ~70% load). db_v4 holds roughly 90 GB of CDS, so the LINEAR upper bound is ~17 B minimizers
# and ~97 GB of hash. Real dedup will land it below that, which is why the number is not
# estimated here: kraken2-inspect reports it and build_db.py prints it.
#
# 230G is therefore headroom over the 97 GB upper bound, not a guess. If --build is OOM-killed,
# move to a highmem node rather than adding --max-db-size: that flag makes Kraken2 SUBSAMPLE
# the minimizer set to fit, silently discarding k-mer space and invalidating every coverage
# denominator in 05_detect.
#
# Time: db_v3's library+build pass took ~4.3 h. db_v4 is ~3x the sequence, so 24 h is set.
#
#   sbatch 02_classify/slurm/build_db.slurm

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load python/3.11.6

# Pin kraken2 by PATH, never `spack load`. Two reasons, both bitten before: `spack load`
# resolves against ~/.pawsey_project rather than the job's --account, and in a
# non-interactive shell `spack` is not even on PATH after `module load spack`, so the load
# fails SILENTLY and the build dies later on a missing binary.
# Two architectures are installed; compute nodes are zen3, login nodes zen2.
K2_ZEN3="/software/projects/pawsey1168/llenzo/setonix/2025.08/software/linux-sles15-zen3/gcc-14.2.0/kraken2-2.1.2-qk2lek64msmwg57n2fowuf5b7etidqaj/bin"
K2_ZEN2="/software/projects/pawsey1168/llenzo/setonix/2025.08/software/linux-sles15-zen2/gcc-14.2.0/kraken2-2.1.2-errlhphjipv3gwjwefhruqqogybieteb/bin"
if "$K2_ZEN3/kraken2-build" --version >/dev/null 2>&1; then
    export PATH="$K2_ZEN3:$PATH"
elif "$K2_ZEN2/kraken2-build" --version >/dev/null 2>&1; then
    echo "note: zen3 kraken2 would not run here, falling back to zen2"
    export PATH="$K2_ZEN2:$PATH"
else
    echo "ABORT: no runnable kraken2-build found at either pinned path."
    exit 1
fi
echo "kraken2-build: $(which kraken2-build)  $(kraken2-build --version 2>&1 | head -1)"

DB_DIR="$SCRATCH/crypt/02_classify/data/silver/db/db_v4"
CDS_DIR="$SCRATCH/crypt/02_classify/data/bronze/cds/pathogen"
mkdir -p "$DB_DIR"

echo "START $(date)"
echo "node: $(hostname)   mem: $(free -g | awk '/^Mem:/{print $2"G total"}')"
# --threads 128, not 32. The 2026-10-16 first attempt asked for --cpus-per-task=32 and used
# --threads 32, but requesting 230G of a 245G node reserves the WHOLE node, so SLURM allocated
# cpu=256 and billed for 256 while kraken2-build used 32. The node is 128 physical cores
# (64 per socket, 2 sockets, SMT 2); 128 matches physical cores, and going to 256 would only
# add SMT siblings, which do not help a compute-bound minimizer pass.
#
# Extra flags pass straight through, so a build can be resumed without redoing the library:
#   sbatch 02_classify/slurm/build_db.slurm --build-only
python 02_classify/02_build/build_db.py \
    --cds-dir "$CDS_DIR" \
    --db-dir "$DB_DIR" \
    --threads 128 "$@"
echo "END $(date)"
