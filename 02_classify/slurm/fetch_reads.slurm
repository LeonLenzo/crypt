#!/bin/bash
#SBATCH --job-name=fetch
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=02_classify/slurm/logs/fetch_%j.out
#SBATCH --error=02_classify/slurm/logs/fetch_%j.err

# 03_cohort: pull a BioProject's reads from ENA.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch --export=ALL,PROJECT=PRJEB83049,PILOT=1 02_classify/slurm/fetch_reads.slurm
#
# PARTITION=copy, not work. This is pure network and filesystem I/O: ENA serves ready-made
# fastq.gz, so nothing is decompressed or recompressed and no core is doing arithmetic. The
# `copy` partition exists for exactly this and keeps a 24h data transfer off the compute nodes.
#
# --mem=8G because the memory request sets the CPU bill and this needs almost none: eight curl
# workers streaming to disk. The OLD route (prefetch + fasterq-dump + pigz) cost ~10,000 SU for
# 14 TB precisely because decompress-then-recompress is CPU-bound.
#
# --workers 24. The "8 not 16" courtesy limit belongs to NCBI, where 16 prefetch streams caused
# connection contention with no throughput gain, and it does NOT apply here: against EBI every
# worker gets its own ~0.7 MB/s because the limit is latency, not contention, so throughput
# scales with streams. Measured on job 51341245: 24 workers held 17.5 MB/s aggregate at
# 0.73 MB/s per stream, exactly linear. The default was left at 8 long after fetch_reads.py's
# DEFAULT_WORKERS moved to 24, so submitting without WORKERS= silently cost two thirds of the
# rate. Raise it only with the retry policy in mind: that job ended with EBI refusing
# connections after ~13 h, and whether sustained concurrency provoked that is still unknown.

set -uo pipefail

module load python/3.11.6

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 02_classify/slurm/logs

: "${PROJECT:?set PROJECT, e.g. --export=ALL,PROJECT=PRJEB83049}"
ARGS="--project $PROJECT --workers ${WORKERS:-24}"
[ -n "${PILOT:-}" ] && ARGS="$ARGS --pilot-only"
[ -n "${LIMIT:-}" ] && ARGS="$ARGS --limit $LIMIT"

echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
python3 02_classify/03_cohort/fetch_reads.py $ARGS
rc=$?
echo "exit     $rc after ${SECONDS}s"
echo "finished $(date -Is)"
exit $rc
