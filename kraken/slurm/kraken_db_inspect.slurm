#!/bin/bash
#SBATCH --job-name=kraken_db_inspect
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=kraken/slurm/logs/kraken_db_inspect_%j.out
#SBATCH --error=kraken/slurm/logs/kraken_db_inspect_%j.err

# Read the built index back out: how many minimizers sit at each taxon, and so at
# which rank each clade's content actually resolves.
#
# This is the one thing the database can answer about its own contents and the
# sketches cannot. Kraken2 stores every distinct minimizer once, against the LCA of
# everything containing it, so a minimizer shared by two Puccinia species is filed
# under Puccinia and the pairwise information is gone. What survives is the marginal
# profile: minimizers exclusive to a species versus pushed up to genus and above.
# A species holding many species-level minimizers with little at its genus node is an
# isolated node — confidently callable with no competitor to pull the call up, which
# is the Phakopsora failure mode. Pairwise similarity needs the FracMinHash sketches
# from kraken_saturation.py --stage sketch instead.
#
# Run from crypt/: sbatch kraken/slurm/kraken_db_inspect.slurm

set -euo pipefail
module load kraken2/2.1.2-qk2lek6
cd /scratch/pawsey1168/llenzo/crypt || exit 1

DB=kraken/build/data/db_v3
OUT=kraken/build/data/db_v3_inspect.tsv

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} db=$DB"
kraken2-inspect --db "$DB" --threads 16 --skip-counts | tee kraken/build/data/db_v3_summary.txt
echo "=== full per-taxon minimizer counts $(date) ==="
kraken2-inspect --db "$DB" --threads 16 --report-zero-counts > "$OUT"
echo "rows: $(wc -l < "$OUT")"
echo "done: $(date)"
