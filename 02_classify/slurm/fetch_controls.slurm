#!/bin/bash
#SBATCH --job-name=fetchctl
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=02_classify/slurm/logs/fetchctl_%j.out
#SBATCH --error=02_classify/slurm/logs/fetchctl_%j.err

# utilities/controls: pull a negative-control study's reads from ENA by accession.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch --dependency=afterany:<cohort fetch> \
#          --export=ALL,ACC="PRJEB29046 PRJNA613154" 02_classify/slurm/fetch_controls.slurm
#
# ONE ENA STREAM AT A TIME. The dependency above is not optional housekeeping: job 51341245
# ended with EBI REFUSING CONNECTIONS after ~13 h of 24 sustained streams, and 623 of 820 runs
# were marked failed in a flat 13 s as a result. Whether sustained concurrency provoked that was
# never established, so a control batch is queued BEHIND the cohort fetch rather than beside it.
# Chain on `afterany`, not `afterok`: a cohort fetch that hits its 24 h wall exits non-zero
# having done its work, and `afterok` would silently cancel the controls for it.
#
# --workers 6 (the script's default), well under fetch_reads' 24, for the same reason. A control
# batch is tens of GB; there is nothing to gain by crowding the archive for it.
#
# PARTITION=copy and --mem=4G: pure network I/O, as fetch_reads.slurm explains at length.
# Nothing is decompressed here, so no core does arithmetic and the memory request -- which sets
# the CPU bill -- is held at the floor.

set -uo pipefail

module load python/3.11.6

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 02_classify/slurm/logs

: "${ACC:?set ACC, e.g. --export=ALL,ACC=\"PRJEB29046 PRJNA613154\"}"
echo "host       $(hostname)"
echo "accessions $ACC"
echo "started    $(date -Is)"

SECONDS=0
# Metadata first and in full, so the log records what the archive said BEFORE any bytes moved.
# That capture is evidence, not noise: PRJEB29046's ENA title for ERR2822761 contradicts the
# paper, and a control's role is the whole reason it is being fetched.
python3 02_classify/utilities/controls/fetch_controls.py --meta-only --accession $ACC

python3 02_classify/utilities/controls/fetch_controls.py \
    --accession $ACC --workers "${WORKERS:-6}"
rc=$?

echo "exit       $rc after ${SECONDS}s"
echo "finished   $(date -Is)"
exit $rc
