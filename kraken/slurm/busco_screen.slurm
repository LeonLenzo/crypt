#!/bin/bash
#SBATCH --job-name=busco_screen
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=24:00:00
#SBATCH --output=kraken/slurm/logs/busco_screen_%j.out
#SBATCH --error=kraken/slurm/logs/busco_screen_%j.err

# Download CDS FASTA + run BUSCO transcriptome mode for all candidates.
# Workers × cpus-per-busco = total CPUs (16 × 8 = 128).
#
# Prerequisites:
#   1. python kraken/screen_refs.py  →  kraken/ref_candidates.tsv
#   2. BUSCO lineage DBs pre-downloaded to $BUSCO_DB_PATH (see below)
#   3. git pull on Setonix to get latest ref_candidates.tsv
#
# Download BUSCO lineage DBs once (if not already present):
#   busco --download fungi_odb10 stramenopiles_odb10 basidiomycota_odb10 \
#         --download_path $BUSCO_DB_PATH --offline
#
# Run from crypt/: sbatch kraken/slurm/busco_screen.slurm

set -euo pipefail

module load python/3.11.6

# BUSCO: load via conda env or module — adjust for your Setonix setup.
# Option A (conda):
#   source /software/projects/pawsey1168/miniconda3/etc/profile.d/conda.sh
#   conda activate busco
# Option B (module, if available):
#   module load busco/5.7.0

SCRATCH="/scratch/pawsey1168/llenzo"
GENOMES_DIR="$SCRATCH/kraken/cds_v2"
BUSCO_OUT="$SCRATCH/kraken/busco"
BUSCO_DB_PATH="$SCRATCH/busco_dbs"
DATASETS_BIN="$SCRATCH/bin/datasets"

# Bootstrap NCBI datasets CLI if not already present
if [ ! -x "$DATASETS_BIN" ]; then
    echo "Downloading NCBI datasets CLI …"
    mkdir -p "$SCRATCH/bin"
    curl -fsSL -o "$DATASETS_BIN" \
        https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$DATASETS_BIN"
fi
export PATH="$SCRATCH/bin:$PATH"
echo "datasets version: $(datasets version)"

mkdir -p kraken/slurm/logs
mkdir -p "$GENOMES_DIR"
mkdir -p "$BUSCO_OUT"
mkdir -p "$BUSCO_DB_PATH"

echo "Starting BUSCO screen: $(date)"
echo "Genomes dir:  $GENOMES_DIR"
echo "BUSCO out:    $BUSCO_OUT"
echo "BUSCO DBs:    $BUSCO_DB_PATH"

python kraken/busco_screen.py \
    --genomes-dir    "$GENOMES_DIR" \
    --busco-out      "$BUSCO_OUT" \
    --busco-db-path  "$BUSCO_DB_PATH" \
    --workers        16 \
    --cpus-per-busco 8

echo "BUSCO screen complete: $(date)"
echo ""
echo "Next steps:"
echo "  1. Copy kraken/busco_scores.tsv back to local (if run on Setonix)"
echo "  2. python kraken/filter_refs.py"
echo "  3. sbatch kraken/slurm/kraken_build.slurm"
