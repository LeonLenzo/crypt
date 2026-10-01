#!/bin/bash
#SBATCH --job-name=smoke_ram
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/smoke_ram_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/smoke_ram_%j.err
# db_v3 smoke test on the 5 ground-truth runs, one config each so both start with a
# cold page cache. 3 of these runs have db_v3 numbers from compare_dbs (49970298) to
# check against: ERR3170248 -> 20.0% unclassified, 47.52% P.striiformis, 2.98% Triticum.
set -uo pipefail
module load python/3.11.6
module load kraken2/2.1.2-qk2lek6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
echo "host=$(hostname) job=$SLURM_JOB_ID config=ram"
free -g | head -2
START=$SECONDS
python kraken/assign/kraken_assign.py \
  --run-list    kraken/assign/data/ground_truth_runs.txt \
  --reads-dir   kraken/select/data/reads \
  --db          kraken/build/data/db_v3 \
  --out-dir     kraken/assign/smoke_ram \
  --reports-dir kraken/assign/smoke_ram/reports \
  --workers 4 --kraken-threads 8
echo "=== CONFIG ram WALL-CLOCK: $((SECONDS-START))s ==="
