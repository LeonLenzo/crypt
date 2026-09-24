#!/bin/bash
#SBATCH --job-name=krassign
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-15
#SBATCH --output=kraken/slurm/logs/krassign_%A_%a.out
#SBATCH --error=kraken/slurm/logs/krassign_%A_%a.err

# Submodule 2 step 3/3: Kraken2 classification of the downloaded cohort.
#
# No new script was needed. kraken/run/kraken_run_assign.py already handles pre-downloaded
# reads via --reads-dir, finds {run}_1.fastq.gz/_2.fastq.gz, skips empty archives,
# and is resumable through kraken_cache.jsonl + kraken_cache_index.txt.
#
# The split/host-removal step is deliberately skipped. Host-read removal was
# measured to change the reported percentage, not what Kraken2 finds, because
# db_v2 is pathogen-only (ratio 1.00x across 5 ground-truth runs).
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch kraken/slurm/kraken_run_assign.slurm
#
# 16 array tasks, strided (todo[array_id::array_count]) so each gets a mixed
# spread of run sizes rather than a contiguous block. ~200 runs per task.
#
# --mem=200G with --workers 8: each kraken2 process loads db_v2 (20 GB) into RAM
# unless --memory-mapping is used. Node has 245 GB. If memory becomes the limit,
# add --memory-mapping to share one on-disk copy instead, at some speed cost.
#
# The kraken2 timeout inside classify.py scales with input size. Do not replace
# it with a flat value: that is precisely what silently truncated 224 files
# during the download step, where large runs blew the limit and small ones did
# not, so the failure looked sporadic.

set -uo pipefail

module load python/3.11.6
module load kraken2/2.1.2-qk2lek6

cd /scratch/pawsey1168/llenzo/crypt || exit 1

DB=kraken/db/build/data/db_v2
READS=kraken/run/kraken_run_select/data/reads
RUNS=kraken/run/kraken_run_select/data/run_list.tsv

echo "host=$(hostname) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID cpus=${SLURM_CPUS_PER_TASK:-?}"
kraken2 --version | head -1

python kraken/run/kraken_run_assign/kraken_run_assign.py \
  --runs-tsv  "$RUNS" \
  --reads-dir "$READS" \
  --db        "$DB" \
  --reports-dir kraken/run/kraken_run_assign/data/reports \
  --workers 8 \
  --kraken-threads 4 \
  --array-id    "$SLURM_ARRAY_TASK_ID" \
  --array-count "$SLURM_ARRAY_TASK_COUNT"
