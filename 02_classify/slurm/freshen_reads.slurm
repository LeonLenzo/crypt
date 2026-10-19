#!/bin/bash
#SBATCH --job-name=freshen_reads
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/freshen_reads_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/freshen_reads_%j.err

# Keep the SRA reads out of /scratch's 21-day purge by touching their atime.
#
# URGENT when submitted 2026-10-19: 3,200 of 3,225 runs had 5 days left. /scratch is Lustre
# mounted WITHOUT noatime, so the purge clock runs on atime and reading 4 KB refreshes a file
# as well as reading all of it — a full pass over 6,029 fastq is ~24 MB of I/O, not 14 TB.
#
# --all, not --strata, and the reason is specific. Of the 3,225 runs on disk, only 1,184 are
# in the 7,688-run literature cohort; the other 2,041 are the old STAT-gated set. Those 2,041
# are ALSO the exact runs db_v3 classified, which makes them the comparison set for the
# db_v3-vs-db_v4 question CLAUDE.md still lists as outstanding. Scratch has no quota and
# 3.7 PB free, so holding them costs nothing but tidiness.
#
# This is NOT a licence to park 14 TB indefinitely, which is what freshen.py's own docstring
# warns against. Once the db_v3/db_v4 comparison is done, the 2,041 non-cohort runs should be
# allowed to expire — they are public SRA accessions, re-pullable from run_list.tsv.
#
# `copy` partition: this is pure I/O. Not the login node, per the standing rule.
#
#   sbatch 02_classify/slurm/freshen_reads.slurm

set -uo pipefail
SCRATCH="/scratch/pawsey1168/llenzo"
cd "$SCRATCH/crypt"
mkdir -p 02_classify/slurm/logs
module load cray-python/3.12.12

echo "START $(date)"
echo "── before ──"
python3 03_kraken/03_select/freshen.py --report 2>&1 | head -8

echo
echo "── freshening every run on disk ──"
# --older-than 0 so nothing is skipped: the default of 10 would leave the 25 runs that were
# touched 13 days ago, and they expire only 3 days after the rest.
python3 -u 03_kraken/03_select/freshen.py --all --older-than 0 2>&1 | tail -15

echo
echo "── after ──"
python3 03_kraken/03_select/freshen.py --report 2>&1 | head -8
echo "END $(date)"
