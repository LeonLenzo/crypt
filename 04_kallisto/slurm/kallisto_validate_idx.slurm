#!/bin/bash
#SBATCH --job-name=kallisto_val
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=58G
#SBATCH --time=01:00:00
#SBATCH --output=04_kallisto/slurm/logs/kallisto_val_%j.out
#SBATCH --error=04_kallisto/slurm/logs/kallisto_val_%j.err

# Confirm every per-host index is complete and readable. `kallisto inspect` parses the index
# header, so a half-written index fails it where a file-existence check would pass -- the
# distinction that matters after a timed-out or cancelled build array.
#
# This belongs in a job, not on a login node: inspect loads multi-GB indices off Lustre.
#
# Run from crypt/: sbatch 04_kallisto/slurm/kallisto_validate_idx.slurm

set -uo pipefail
module load spack/default
spack load kallisto@0.50.1
cd /scratch/pawsey1168/llenzo/crypt || exit 1

D=04_kallisto/02_build/data/idx
ok=0; bad=0
for f in "$D"/*.idx; do
  gb=$(stat -c%s "$f" | awk '{printf "%.2f", $1/1e9}')
  if kallisto inspect "$f" >/dev/null 2>&1; then
    ok=$((ok+1)); printf "  valid    %-30s %8s GB\n" "$(basename "$f")" "$gb"
  else
    bad=$((bad+1)); printf "  INVALID  %-30s %8s GB\n" "$(basename "$f")" "$gb"
  fi
done
echo "valid=$ok invalid=$bad"
[ "$bad" -eq 0 ] || exit 1
