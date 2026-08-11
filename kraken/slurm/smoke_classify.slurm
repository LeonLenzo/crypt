#!/bin/bash
#SBATCH --account=pawsey1168
#SBATCH --job-name=krak_smoke
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --partition=work
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/output/classify/logs/smoke_%j.log

module load python/3.11.6
module load spack/0.23.1
spack load kraken2

cd /scratch/pawsey1168/llenzo/crypt

# 5-run smoke test: 2 workers × 4 threads, limit 5 runs.
# Verify DB path, gzipped reads, and output format before the full control job.
python -u kraken/classify.py \
    --run-list kraken/control/output/data/run_ids.txt \
    --reads-dir /scratch/pawsey1168/llenzo/kraken_reads/control \
    --db /scratch/pawsey1168/llenzo/kraken/db_pathogens \
    --workers 2 \
    --kraken-threads 4 \
    --limit 5
