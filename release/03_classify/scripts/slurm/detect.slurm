#!/bin/bash
#SBATCH --job-name=detect
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=release/03_classify/scripts/slurm/logs/detect_%j.out
#SBATCH --error=release/03_classify/scripts/slurm/logs/detect_%j.err
#
# Re-decide the calls for BOTH arms and rebuild the supplementary tables.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/03_classify/scripts/slurm/detect.slurm
#
# This exists because detect.py reads the 1,560 reports per arm on scratch, so it can only run
# on the machine holding them, and a 3,120-file pass is not login-node work.
#
# ARGUMENTS ARE FLAGS, NOT ENVIRONMENT VARIABLES: Setonix defaults sbatch to --export=NONE.
# This script takes none; it always does both arms, because the two tables are only comparable
# when they were produced by the same code, and the thing being fixed is in that code.
#
# --mem=16G and 8 cores. The memory request is the CPU bill (about 1.9 GB per physical core),
# so 16G bills roughly 9 cores. detect.py is a single-threaded pass over small TSVs plus the
# inspect table; the hash is not loaded and 100G would bill the node for nothing.
#
# BOTH ARMS, THEN THE TABLES, IN ORDER. make_supp_tables.py reads detections_db2.tsv and
# compare_arms.py reads both, so a partial run leaves S5 describing one arm and the comparison
# describing another. -e is not set so that a failure is reported with its step named rather
# than vanishing.

set -uo pipefail
cd /scratch/pawsey1168/llenzo/crypt || exit 1
tools/repo_sync.sh require-current || exit 1
mkdir -p release/03_classify/scripts/slurm/logs
module load python/3.11.6

echo "host     $(hostname)"
echo "started  $(date -Is)"
echo "commit   $(git rev-parse --short HEAD)"
echo

rc=0
for arm in db1p db2; do
    echo "=== detect.py --arm $arm ==="
    python3 release/03_classify/scripts/detect.py --arm "$arm" || { rc=1; echo "FAILED: detect $arm"; }
    echo
done

if [ $rc -eq 0 ]; then
    echo "=== make_supp_tables.py ==="
    python3 release/03_classify/scripts/make_supp_tables.py || { rc=1; echo "FAILED: make_supp_tables"; }
    echo
    echo "=== compare_arms.py ==="
    python3 release/03_classify/scripts/compare_arms.py || { rc=1; echo "FAILED: compare_arms"; }
else
    echo "SKIPPING the tables: a detect arm failed, and a table built on one arm would be worse than none."
fi

echo
echo "finished $(date -Is)  rc=$rc"
exit $rc
