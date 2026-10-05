#!/bin/bash
#SBATCH --job-name=kallisto_idx
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=58G
#SBATCH --time=06:00:00
#SBATCH --array=0-18
#SBATCH --output=04_kallisto/slurm/logs/kallisto_idx_%A_%a.out
#SBATCH --error=04_kallisto/slurm/logs/kallisto_idx_%A_%a.err

# One kallisto index per host, for 04_kallisto's competitive read assignment.
#
# Array index -> host, from 04_kallisto/01_select/data/strata.tsv (19 hosts, every host with a
# flagged detection). Each index holds, as tagged CDS: the host, every species Kraken2 saw in
# that host at ANY level, and the measured >=1% k-mer containment partners of the flagged ones.
# Rationale in 04_kallisto/README.md and the kallisto_build.py docstring.
#
# Sizes from --dry-run: 0.26 GB (Eleusine coracana) to 3.38 GB (Triticum aestivum) of CDS.
# Measured peak RSS is ~10x the CDS size (13.7 GB for a 1.3 GB host, 2.9 GB for Eleusine), so
# wheat needs ~34 GB. 58G is the 32-core quantum on a 230G node and fits 4 tasks per node;
# an earlier 115G request was 3x oversized and halved throughput for nothing.
#
# `spack load` prepends its own python to PATH, so cray-python is loaded FIRST and the scripts
# are stdlib-only by design. If a step here ever needs numpy, load spack first or use srun.
#
# Run from crypt/:
#   sbatch 04_kallisto/slurm/kallisto_build_index.slurm              # all 19
#   sbatch --array=2 04_kallisto/slurm/kallisto_build_index.slurm    # one host
#
# Prerequisites:
#   python 04_kallisto/01_select/kallisto_select.py            -> 01_select/data/strata.tsv
#   sbatch 04_kallisto/slurm/kallisto_neighbours.slurm         -> 02_build/data/neighbours.tsv

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

python3 04_kallisto/02_build/kallisto_build.py \
    --host "$HOST" --threads "${SLURM_CPUS_PER_TASK:-8}"
rc=$?
echo
echo "exit $rc"
exit $rc
