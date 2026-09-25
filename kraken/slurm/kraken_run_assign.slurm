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
# No new script was needed. kraken/assign.py already handles pre-downloaded
# reads via --reads-dir, finds {run}_1.fastq.gz/_2.fastq.gz, skips empty archives,
# and is resumable through kraken_cache.jsonl + kraken_cache_index.txt.
#
# The split/host-removal step is deliberately skipped. Host-read removal was
# measured to change the reported percentage, not what Kraken2 finds, because
# db_v2 was pathogen-only (ratio 1.00x across 5 ground-truth runs). db_v3 holds host
# CDS, but the same finding stands: see crypt-host-removal-not-needed.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch kraken/slurm/kraken_run_assign.slurm
#
# 16 array tasks, strided (todo[array_id::array_count]) so each gets a mixed
# spread of run sizes rather than a contiguous block. ~200 runs per task.
#
# --mem=200G with --workers 4: each kraken2 process loads db_v3's 27.4 GB hash.k2d
# into RAM, so 4 x 27.4 = 110 GB on a 245 GB node. db_v2 was 6.3 GB and allowed 8.
#
# Do NOT use --memory-mapping here, despite what the db_v2-era note said. Benchmarked
# 2026-09-25 on the 5 ground-truth runs: the RAM config finished in 3m42s, while the
# mmap config wrote zero reports in 13 minutes with kraken2 pinned at 1.6% CPU. The
# index lives on Lustre and Kraken2's hash lookups are random access, so every page
# fault becomes a filesystem RPC. It does not amortise; it is I/O bound throughout.
#
# The kraken2 timeout inside classify.py scales with input size. Do not replace
# it with a flat value: that is precisely what silently truncated 224 files
# during the download step, where large runs blew the limit and small ones did
# not, so the failure looked sporadic.

set -uo pipefail

module load python/3.11.6
module load kraken2/2.1.2-qk2lek6

cd /scratch/pawsey1168/llenzo/crypt || exit 1

DB=kraken/build/data/db_v3
READS=kraken/select/data/reads
RUNS=kraken/select/data/run_list.tsv

# --tmp-dir must be set explicitly. kraken_assign.py defaults it to $MYSCRATCH, and
# MYSCRATCH follows PAWSEY_PROJECT, which is still the deprecated fl3 on this account:
# /scratch/fl3/llenzo. The 2026-09-25 run lost 443 of 3,225 runs to
# "[Errno 2] No such file or directory: /scratch/fl3/llenzo/kraken_tmp/<run>", the
# per-run dir disappearing after it was created, consistent with fl3 being decommissioned
# underneath the job. Never let this default through on Setonix.
TMP=/scratch/pawsey1168/llenzo/kraken_tmp

echo "host=$(hostname) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID cpus=${SLURM_CPUS_PER_TASK:-?}"
kraken2 --version | head -1

python kraken/assign/kraken_assign.py \
  --runs-tsv  "$RUNS" \
  --reads-dir "$READS" \
  --db        "$DB" \
  --reports-dir kraken/assign/data/reports \
  --tmp-dir     "$TMP" \
  --workers 4 \
  --kraken-threads 8 \
  --array-id    "$SLURM_ARRAY_TASK_ID" \
  --array-count "$SLURM_ARRAY_TASK_COUNT"
