#!/bin/bash
#SBATCH --job-name=fetch_pf
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-2
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/f125_%A_%a.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/f125_%A_%a.err

# ONE GENOME PER ARRAY TASK. The serial fetcher did ~3 in a few minutes; the limit is NCBI
# per-connection latency, not the data movers, so independent downloads parallelise cleanly.
# 20 concurrent is well under any NCBI contention threshold (the datasets CLI is one stream
# each). Into the SAME store the pilot used.
set -uo pipefail
R=/scratch/pawsey1168/llenzo/crypt
G=$R/02_classify/data/bronze/genomes
DATASETS=/software/projects/pawsey1168/llenzo/bin/datasets
mkdir -p "$G"
row=$(tail -n +2 "$R/research/predict/data/pfermentans_alt.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -z "$row" ] && exit 0
ACC=$(echo "$row" | cut -f1)
if [ -n "$(find "$G/$ACC" -name "*_genomic.fna" 2>/dev/null)" ]; then echo "skip $ACC (present)"; exit 0; fi
rm -rf "$G/$ACC"   # empty dir from a cancelled fetch would block the download
echo "fetch $ACC"
if "$DATASETS" download genome accession "$ACC" --include genome --filename "/tmp/$ACC.zip" >/dev/null 2>&1 \
   && unzip -oq "/tmp/$ACC.zip" -d "$G/$ACC.tmp" && mv "$G/$ACC.tmp" "$G/$ACC"; then
  echo "ok $ACC $(find "$G/$ACC" -name "*_genomic.fna" | wc -l) fasta"
else
  echo "FAILED $ACC"; rm -rf "$G/$ACC.tmp"; rm -f "/tmp/$ACC.zip"; exit 1
fi
rm -f "/tmp/$ACC.zip"
