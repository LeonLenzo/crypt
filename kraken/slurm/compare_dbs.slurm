#!/bin/bash
#SBATCH --job-name=cmp_dbv2v3
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=230G
#SBATCH --time=6:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/cmp_dbv2v3_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/cmp_dbv2v3_%j.err
set -euo pipefail
module load spack/0.23.1 && spack load kraken2
BASE=/scratch/pawsey1168/llenzo/crypt/kraken
READS=$BASE/select/data/reads
OUT=$BASE/assign/data/compare_v2v3; mkdir -p $OUT
DBV2=$BASE/build/data/db_v2
DBV3=$BASE/build/data/db_v3
for RUN in ERR3170248 ERR3170234 ERR3170229; do
  for TAG in v2 v3; do
    DB=$([ $TAG = v2 ] && echo $DBV2 || echo $DBV3)
    echo "=== $RUN $TAG $(date) ==="
    kraken2 --db $DB --threads 32 --confidence 0.15 --minimum-hit-groups 3 \
      --paired $READS/${RUN}_1.fastq.gz $READS/${RUN}_2.fastq.gz \
      --report $OUT/${RUN}.${TAG}.report --output /dev/null 2>$OUT/${RUN}.${TAG}.log
  done
done
echo "=== ALL DONE $(date) ==="
