#!/bin/bash
#SBATCH --account=pawsey1168
#SBATCH --job-name=dl_control
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --partition=work
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/output/download/logs/dl_control_%j.log

module load python/3.11.6

cd /scratch/pawsey1168/llenzo/crypt

python kraken/download.py \
    --run-list kraken/control/output/data/run_ids.txt \
    --reads-dir /scratch/pawsey1168/llenzo/kraken_reads/control \
    --workers 16 \
    --n-reads 500000
