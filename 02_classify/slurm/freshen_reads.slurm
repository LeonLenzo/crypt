#!/bin/bash
#SBATCH --job-name=freshen_reads
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/freshen_reads_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/freshen_reads_%j.err

# Keep the COHORT reads out of /scratch's 21-day purge by touching their atime.
#
# /scratch is Lustre mounted WITHOUT noatime, so the purge clock runs on atime and reading
# 4 KB refreshes a file as well as reading all of it: a full pass over the 2,442 fastq in the
# cohort tree is ~10 MB of I/O, not 4 TB.
#
# ── WHAT CHANGED, 2026-10-20 ─────────────────────────────────────────────────────────────────
# This job used to call `03_kraken/03_select/freshen.py`, whose READS constant was the FLAT
# pre-migration directory. `migrate_reads.py` moved the cohort into reads/<BioProject>/ and
# preserved atime on purpose, so the purge clock carried across unbroken -- and the old script
# kept reporting a healthy tree while freshening nothing the cohort still used. The 4.1 TiB
# that matters survived only because classify jobs happened to read it.
#
# It now calls `02_classify/03_cohort/freshen.py`, which takes its root from `_layout.READS`.
# The path cannot drift from the one the rest of the module writes to.
#
# ── SCOPE: the cohort tree only ──────────────────────────────────────────────────────────────
# The flat tree `03_kraken/03_select/data/reads` still holds ~9.1 TB: the 2,041 non-cohort
# STAT-gated runs. They are NOT freshened and are meant to EXPIRE. They were held open for the
# db_v3-vs-db_v4 comparison, which is closed -- db_v4 is validated against the off-host null
# and is the database in use. Every one of those runs is a public SRA accession, re-pullable
# from run_list.tsv.
#
# --older-than 0 so nothing is skipped; the default of 10 leaves recently touched files to
# expire a few days after the rest, which is the staggered-expiry trap.
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
python3 02_classify/03_cohort/freshen.py --report 2>&1 | head -12

echo
echo "── freshening every run in the cohort tree ──"
python3 -u 02_classify/03_cohort/freshen.py --all --older-than 0 2>&1 | tail -15

echo
echo "── after ──"
python3 02_classify/03_cohort/freshen.py --report 2>&1 | head -12
echo "END $(date)"
