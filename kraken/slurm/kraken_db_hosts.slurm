#!/bin/bash
#SBATCH --job-name=kraken_db_hosts
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=kraken/slurm/logs/kraken_db_hosts_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_db_hosts_%j.err

# Select and download host plant CDS for db_v3. Reads the resolved host taxids from
# run_list.tsv, picks assemblies from Ensembl Plants (cereal pangenomes) and NCBI
# (the tail), and downloads CDS into cds/host/. See kraken_db_hosts.py.
#
# Prerequisite: kraken_run_select.py has produced run_list.tsv with host_taxid.
# Run from crypt/: sbatch kraken/slurm/kraken_db_hosts.slurm

set -euo pipefail
module load python/3.11.6
SCRATCH="/scratch/pawsey1168/llenzo"
DATASETS_BIN="$SCRATCH/bin/datasets"
if [ ! -x "$DATASETS_BIN" ]; then
    mkdir -p "$SCRATCH/bin"
    curl -fsSL -o "$DATASETS_BIN" \
        https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$DATASETS_BIN"
fi
export PATH="$SCRATCH/bin:$PATH"
mkdir -p kraken/slurm/logs

echo "Starting host selection + download: $(date)"
# Downloads the committed host_candidates.tsv selection (no --reselect), so the
# cluster fetches exactly what was reviewed locally, not a fresh Ensembl/NCBI query.
python kraken/db/kraken_db_hosts.py --download --workers 16
echo "Host download complete: $(date)"
echo "Next: sbatch kraken/slurm/kraken_db_build.slurm"
