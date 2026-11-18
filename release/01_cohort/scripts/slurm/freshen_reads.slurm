#!/bin/bash
#SBATCH --job-name=freshen_reads
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/release/01_cohort/scripts/slurm/logs/freshen_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/release/01_cohort/scripts/slurm/logs/freshen_%j.err
#
# Keep the cohort reads out of /scratch's 21-day purge by touching their atime. Run it by hand
# when something has just landed; the weekly schedule is freshen_reads.scron.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/01_cohort/scripts/slurm/freshen_reads.slurm
#
# PARTITION=copy. Reading one 4 KB block per file is ~12 MB of I/O over the whole cohort and
# no core does arithmetic; --mem=4G because the memory request sets the CPU bill.
#
# WHAT KILLED THE LAST ONE. research/02_classify/slurm/freshen_reads.slurm called a script
# whose reads root came from `_layout.READS`, pointing into the research tree. The cohort had
# moved to release/01_cohort/raw/reads, so it renewed 8.24 TB meant to expire and protected
# none of the 4.90 TB that matters, while reporting a healthy tree and exiting 0. It last ran
# 2026-10-19; the purge took 991 of 1,560 runs on 2026-11-13. freshen_reads.py now derives its
# root from its own location and fails on low cohort coverage, so the same drift cannot pass
# silently. That old pair is superseded; do not resurrect it.
exec /scratch/pawsey1168/llenzo/crypt/release/01_cohort/scripts/slurm/freshen_reads.sh
