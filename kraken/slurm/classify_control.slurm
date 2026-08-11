#!/bin/bash
#SBATCH --account=pawsey1168
#SBATCH --job-name=krak_control
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=work
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/output/classify/logs/krak_control_%j.log
#SBATCH --dependency=singleton

module load python/3.11.6
module load spack/0.23.1
spack load kraken2
module load gcc-native/12
module load kraken2/2.1.2-qk2lek6

cd /scratch/pawsey1168/llenzo/crypt

# Reads pre-downloaded by download_control.slurm — classify from local gzipped files.
# 32 CPUs: 8 workers × 4 kraken2 threads each.
python kraken/classify.py \
    --run-list kraken/control/output/data/run_ids.txt \
    --reads-dir /scratch/pawsey1168/llenzo/kraken_reads/control \
    --db /scratch/pawsey1168/llenzo/kraken/db_pathogens \
    --workers 8 \
    --kraken-threads 4
