#!/bin/bash
#SBATCH --account=pawsey1168
#SBATCH --job-name=bench_dl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=work
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/output/benchmark_setonix_%j.log

module load python/3.11.6

cd /scratch/pawsey1168/llenzo/crypt

python kraken/benchmark_download.py \
    --runs 5 \
    --methods ena-ftp ena-https \
    --workers 1 8 \
    --tmp-dir /scratch/pawsey1168/llenzo/kraken_tmp/bench \
    --out kraken/output/benchmark_setonix.tsv
