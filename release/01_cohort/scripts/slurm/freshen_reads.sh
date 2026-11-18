#!/bin/bash
# The body of the cohort reads keep-alive. ONE copy, because it has two callers:
#
#   freshen_reads.slurm   manual `sbatch`, #SBATCH supplies the resources
#   freshen_reads.scron   the weekly scrontab entry, #SCRON supplies the same resources
#
# scrontab runs a plain command rather than a script with #SBATCH directives in it, so the two
# callers cannot be one file. They can share this one, and must: a keep-alive that behaves
# differently depending on who started it is the next version of the bug this replaces.
#
# Setonix submits with --export=NONE, so this loads its own module and assumes no environment.
set -uo pipefail
module load python/3.11.6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p release/01_cohort/scripts/slurm/logs
tools/repo_sync.sh require-current || exit 1
echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
# --older-than 0 so nothing is skipped. The old script's default of 10 left recently touched
# files to expire a few days after the rest: the tree stops having one deadline and starts
# having a dozen.
python3 -u release/01_cohort/scripts/freshen_reads.py --all --older-than 0
rc=$?
echo "exit     $rc after ${SECONDS}s"
echo "finished $(date -Is)"
exit $rc
