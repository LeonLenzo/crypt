#!/bin/bash
#SBATCH --job-name=fetch_cds
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/fetch_cds_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/fetch_cds_%j.err

# Put a CDS FASTA on disk for every accession in ref_candidates.tsv (4,286 as of 2026-10-16).
# Reconciles against Acacia FIRST: 1,753 of them were downloaded for db_v3 and are in the
# project bucket, and pulling those is free where re-pulling from NCBI is not.
#
# `copy` partition, not `work`: this job moves ~50 GB and does no computation. Running it on
# `work` wastes a compute allocation on I/O wait, and running it on a login node is forbidden.
#
# Run from the repo root on Setonix:
#   sbatch 02_classify/slurm/fetch_cds.slurm

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load python/3.11.6
module load aws-cli/2.13.0

# The aws-cli module is a Singularity wrapper that overrides --home to
# $MYSOFTWARE/setonix/<ver>/.shpc_home, so it reads that tree's .aws and NOT ~/.aws.
# Credentials in the real home are invisible to it; the symptom is "Unable to locate
# credentials". Fail here rather than 4,286 accessions later.
if ! aws s3 ls --profile pawsey1168 --endpoint-url https://projects.pawsey.org.au \
        s3://pawsey1168-kraken-db/ >/dev/null 2>&1; then
    echo "ABORT: cannot read the project bucket as profile pawsey1168."
    echo "  Credentials must live in \$MYSOFTWARE/setonix/*/.shpc_home/.aws, not ~/.aws."
    exit 1
fi

# NCBI datasets CLI, bootstrapped to scratch (no module provides it).
DATASETS_BIN="$SCRATCH/bin/datasets"
if [ ! -x "$DATASETS_BIN" ]; then
    echo "fetching NCBI datasets CLI"
    mkdir -p "$SCRATCH/bin"
    curl -sSL -o "$DATASETS_BIN" \
        https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$DATASETS_BIN"
fi
export PATH="$SCRATCH/bin:$PATH"
echo "datasets version: $(datasets version)"

CDS_DIR="$SCRATCH/crypt/02_classify/data/bronze/cds/pathogen"
mkdir -p "$CDS_DIR"

echo "START $(date)"
# --workers 8, not 16: measured. Past ~8, NCBI throttles and the retries cost more than the
# added concurrency returns.
python 02_classify/02_build/fetch_cds.py \
    --cds-dir "$CDS_DIR" \
    --profile pawsey1168 \
    --workers 8
echo "END $(date)"
