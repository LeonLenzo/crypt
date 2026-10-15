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

# Check against the hosts 01_select actually needs, not against whatever files happen to be in
# the directory. Globbing *.idx alone passes when an index is MISSING, which is the failure that
# matters after an OOM kill: the broken file gets deleted and the glob stops seeing it.
#
# VALIDATE_HOSTS overrides the list with space-separated slugs, for checking one rebuilt index
# instead of re-inspecting all 19. inspect costs ~1.5 min per index, so a full pass is ~30 min
# and is only worth it after a whole-array build.
#   sbatch --export=ALL,VALIDATE_HOSTS=zea_mays 04_kallisto/slurm/kallisto_validate_idx.slurm
module load cray-python/3.11.7
if [ -n "${VALIDATE_HOSTS:-}" ]; then
  read -r -a HOSTS <<< "$VALIDATE_HOSTS"
  echo "checking only: ${HOSTS[*]}"
else
  mapfile -t HOSTS < <(python3 - <<'PY'
import csv
for h in sorted({r["host"] for r in csv.DictReader(
        open("04_kallisto/01_select/data/strata.tsv"), delimiter="\t")}):
    print(h.lower().replace(" x ", "_x_").replace(" ", "_"))
PY
)
fi

ok=0; bad=0; missing=0
for h in "${HOSTS[@]}"; do
  f="$D/$h.idx"
  if [ ! -f "$f" ]; then
    missing=$((missing+1)); printf "  MISSING  %-30s\n" "$h.idx"; continue
  fi
  gb=$(stat -c%s "$f" | awk '{printf "%.2f", $1/1e9}')
  if kallisto inspect "$f" >/dev/null 2>&1; then
    ok=$((ok+1)); printf "  valid    %-30s %8s GB\n" "$h.idx" "$gb"
  else
    bad=$((bad+1)); printf "  INVALID  %-30s %8s GB\n" "$h.idx" "$gb"
  fi
done
echo "expected=${#HOSTS[@]} valid=$ok invalid=$bad missing=$missing"
[ "$bad" -eq 0 ] && [ "$missing" -eq 0 ] || exit 1
