#!/bin/bash
#SBATCH --job-name=dedup_pred
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/dedup_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/dedup_%j.err

# Repair the 36 predictions that two concurrent array jobs wrote twice. See the docstring
# in dedup_predictions.py for the cause and for why truncation is lossless.
#
# Batched rather than run on the login node: it reads ~4.5 GB across 405 files, which is
# exactly the kind of I/O the login nodes are not for.
#
# cray-python, NOT python/3.11.6. The login nodes still default to 3.6.15 and this script
# uses PEP 604 unions in its annotations.

set -uo pipefail
module load cray-python/3.12.12
P=/scratch/pawsey1168/llenzo/crypt/research/predict
cd "$P" || exit 1

MODE="${MODE:---check}"
echo "START $(date -Is)  mode=$MODE"
python3 dedup_predictions.py "$MODE"
rc=$?
echo "END $(date -Is)  rc=$rc"
exit $rc
