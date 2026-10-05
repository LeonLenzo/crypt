#!/bin/bash
#SBATCH --job-name=kallisto_quant
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=58G
#SBATCH --time=12:00:00
#SBATCH --array=0-18
#SBATCH --output=04_kallisto/slurm/logs/kallisto_quant_%A_%a.out
#SBATCH --error=04_kallisto/slurm/logs/kallisto_quant_%A_%a.err

# Pseudoalign the stratified set: 174 runs, 0.49 TB of reads, against the per-host indices.
#
# Array index -> host, so each task handles every selected run for one host and loads that
# host's index once instead of once per run. Wheat carries the most runs, which is why the
# walltime is 12h rather than sized on the average.
#
# Two kallisto passes per run (quant for abundances, bus for equivalence classes). The second
# is required: `kallisto quant` in 0.50.1 has no EC output, so unique-versus-shared reads are
# otherwise unavailable. See the kallisto_quant.py docstring.
#
# Resume is integrity-checked, not existence-checked, so a re-run after a timeout is cheap and
# safe. Reads are under a 21-day atime purge; quant refreshes what it touches, but run
# 03_kraken/03_select/freshen.py --report first if it has been a while.
#
# Run from crypt/:
#   sbatch 04_kallisto/slurm/kallisto_quant.slurm
#   sbatch --array=0 04_kallisto/slurm/kallisto_quant.slurm    # one host
#
# Prerequisites: 01_select/data/strata.tsv and all 19 indices in 02_build/data/idx/.

set -uo pipefail

module load cray-python/3.11.7
module load spack/default
spack load kallisto@0.50.1

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 04_kallisto/slurm/logs

HOST=$(python3 - "$SLURM_ARRAY_TASK_ID" <<'PY'
import csv, sys
hosts = sorted({r["host"] for r in
                csv.DictReader(open("04_kallisto/01_select/data/strata.tsv"), delimiter="\t")})
i = int(sys.argv[1])
print(hosts[i] if i < len(hosts) else "", end="")
PY
)
[ -z "$HOST" ] && { echo "no host for array index ${SLURM_ARRAY_TASK_ID}"; exit 0; }

echo "host      : $(hostname)"
echo "job       : ${SLURM_ARRAY_JOB_ID:-none}_${SLURM_ARRAY_TASK_ID:-none}"
echo "target    : $HOST"
echo "kallisto  : $(kallisto version 2>&1 | head -1)"
echo "cpus      : ${SLURM_CPUS_PER_TASK:-?}"
echo

python3 04_kallisto/03_quant/kallisto_quant.py \
    --host "$HOST" --threads "${SLURM_CPUS_PER_TASK:-8}"
rc=$?
echo
echo "exit $rc"
exit $rc
