#!/bin/bash
#SBATCH --job-name=predict_fetch
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/fetch_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/fetch_%j.err

# Fetch the pilot genomes with the NCBI datasets CLI, into the SAME tree the 184 benchmark
# genomes already live in. That tree is at the PRE-research/-move path; keeping one genome
# store matters more than the path being tidy, and moving 184 genomes is a separate job.
set -uo pipefail
R=/scratch/pawsey1168/llenzo/crypt
G=$R/02_classify/data/bronze/genomes
DATASETS=/software/projects/pawsey1168/llenzo/bin/datasets
mkdir -p "$G" "$R/research/predict/logs"
cd "$G" || exit 1
echo "START $(date -Is)"
ok=0; skip=0; fail=0
while IFS=$'\t' read -r acc taxon mb; do
  [ "$acc" = "accession" ] && continue
  if [ -n "$(find "$G/$acc" -name "*_genomic.fna" 2>/dev/null)" ]; then echo "  skip  $acc (fasta present)"; skip=$((skip+1)); continue; fi
  rm -rf "$G/$acc"   # empty dir from a cancelled fetch would otherwise block the download
  echo "  fetch $acc  $taxon  ${mb} Mb"
  if "$DATASETS" download genome accession "$acc" --include genome \
       --filename "/tmp/$acc.zip" >/dev/null 2>&1 \
     && unzip -oq "/tmp/$acc.zip" -d "$G/$acc.tmp" \
     && mv "$G/$acc.tmp" "$G/$acc"; then
    n=$(find "$G/$acc" -name "*_genomic.fna" | wc -l)
    echo "        ok, $n fasta"; ok=$((ok+1))
  else
    echo "        FAILED"; rm -rf "$G/$acc.tmp"; fail=$((fail+1))
  fi
  rm -f "/tmp/$acc.zip"
done < "$R/research/predict/data/pilot10.tsv"
echo "fetched=$ok skipped=$skip failed=$fail"
echo "END $(date -Is)"
[ "$fail" -eq 0 ]
