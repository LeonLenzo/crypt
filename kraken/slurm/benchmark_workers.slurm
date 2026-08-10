#!/bin/bash
#SBATCH --account=pawsey1168
#SBATCH --job-name=bench_workers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=work
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/output/benchmark_workers_%j.log

module load python/3.11.6

cd /scratch/pawsey1168/llenzo/crypt

mkdir -p kraken/output

python kraken/benchmark_download.py \
    --runs 15 \
    --methods ena-https \
    --workers 8 16 24 32 48 \
    --tmp-dir /scratch/pawsey1168/llenzo/kraken_tmp/bench \
    --out kraken/output/benchmark_workers.tsv
