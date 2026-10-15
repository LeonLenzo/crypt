#!/bin/bash
#SBATCH --job-name=kraken_annotate
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=03_kraken/utilities/annotate/logs/annotate_%j.out
#SBATCH --error=03_kraken/utilities/annotate/logs/annotate_%j.err

# Predict CDS for an unannotated assembly by projecting annotated congeners' proteins onto it
# with miniprot, so it can enter the CDS-based index as a competitor. See kraken_annotate.py.
#
# First target: Puccinia hordei (barley leaf rust) GCA_062009444.1, guided by the three
# annotated wheat rusts whose CDS are on disk (P. triticina, P. striiformis, P. graminis).
# This is the 04_kallisto blind-spot test: does barley leaf rust surface, currently misattributed
# to P. triticina.
#
# Builds miniprot on first run (into $MYSOFTWARE/.../bin), which needs a compiler, so it happens
# in the job and never on a login node.
#
# Run from crypt/:
#   sbatch 03_kraken/slurm/kraken_annotate.slurm
#   TARGET=GCA_xxx GUIDES=GCA_a,GCA_b sbatch 03_kraken/slurm/kraken_annotate.slurm

set -uo pipefail
module load cray-python/3.11.7
module load gcc/12.2.0 2>/dev/null || module load gcc 2>/dev/null || true
export PATH="/software/projects/pawsey1168/llenzo/bin:$PATH"   # datasets CLI + miniprot
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 03_kraken/utilities/annotate/logs

TARGET="${TARGET:-GCA_062009444.1}"
GUIDES="${GUIDES:-GCA_977970215.1,GCA_025169535.1,GCA_008522505.1}"

echo "host   : $(hostname)"
echo "job    : ${SLURM_JOB_ID:-none}"
echo "target : $TARGET"
echo "guides : $GUIDES"
echo

python3 03_kraken/utilities/annotate/kraken_annotate.py \
    --target "$TARGET" --guides "$GUIDES" --threads "${SLURM_CPUS_PER_TASK:-8}"
echo "exit $?"
