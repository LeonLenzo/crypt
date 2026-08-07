#!/bin/bash
#SBATCH --job-name=kraken_unmask_debug
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0:30:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_unmask_debug_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/kraken_unmask_debug_%j.err

set -euo pipefail

module load spack/0.23.1
spack load kraken2

DB="/scratch/pawsey1168/llenzo/kraken/db_unmasked/db"
TMP="/scratch/pawsey1168/llenzo/kraken_tmp"
REPORT="/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/debug_report.txt"

echo "kraken2 version: $(kraken2 --version 2>&1 | head -1)"
echo "DB: $DB"
echo "hash.k2d size: $(du -sh $DB/hash.k2d | cut -f1)"
echo ""

echo "--- Running kraken2 on SRR13349969 (paired, existing temp FASTQs) ---"
kraken2 \
    --db "$DB" \
    --memory-mapping \
    --report "$REPORT" \
    --output /dev/null \
    --threads 4 \
    --paired \
    "$TMP/SRR13349969/r1.fastq" \
    "$TMP/SRR13349969/r2.fastq"

echo ""
echo "--- Report (top 10 lines) ---"
head -10 "$REPORT"
