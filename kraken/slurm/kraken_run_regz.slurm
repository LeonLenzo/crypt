#!/bin/bash
#SBATCH --job-name=krsel_regz
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --array=0-15
#SBATCH --output=kraken/slurm/logs/krsel_regz_%A_%a.out
#SBATCH --error=kraken/slurm/logs/krsel_regz_%A_%a.err

# Repair pass for reads left behind by the krsel_full chain (2026-09-14).
#
# Two classes of damage, both caused by compression not completing:
#
#   pair   a .gz AND its raw .fastq both exist. gzip only unlinks the source on
#          success, so a surviving raw means the .gz is truncated. 228 of these;
#          3 of 3 sampled failed `gzip -t`. Cause was kraken_run_select.py's
#          flat `timeout=1800` on gzip (four runs of 174-230 GB each produced a
#          .gz of almost exactly 4.1 GB) plus SIGKILL at the 24h job wall.
#   orphan a raw .fastq with no .gz at all. 131 of these. fasterq-dump finished
#          but compression never produced anything.
#
# Both are fully recoverable because the raw reads survived in every case.
#
# This job verifies before it destroys anything: for a pair it tests the .gz
# first, and if the .gz turns out to be valid it keeps it and drops the
# redundant raw instead of recompressing. Every file is handled independently
# and the script is idempotent, so a timed-out array task can simply be rerun.
#
# Run from crypt/: sbatch kraken/slurm/kraken_run_regz.slurm
#
# Array of 8 tasks, each taking every 8th file, so one 230 GB run cannot block
# the rest and a single task dying loses only its own slice.

set -uo pipefail

READS=/scratch/pawsey1168/llenzo/crypt/kraken/output/run/select/data/reads
THREADS=${SLURM_CPUS_PER_TASK:-8}
TASK=${SLURM_ARRAY_TASK_ID:-0}
NTASK=${SLURM_ARRAY_TASK_COUNT:-1}

cd "$READS" || { echo "FATAL: cannot cd to $READS"; exit 1; }
command -v pigz >/dev/null || { echo "FATAL: pigz not found"; exit 1; }

# Build the worklist fresh each run rather than trusting a stale list, so the
# job stays correct if it is rerun after a partial pass.
mapfile -t ALL < <(ls *.fastq 2>/dev/null | grep -v '\.gz$' | sort)
echo "task $TASK/$NTASK  candidate raw files: ${#ALL[@]}  threads=$THREADS"

recovered=0; kept=0; failed=0; bytes_before=0; bytes_after=0

for i in "${!ALL[@]}"; do
    (( i % NTASK != TASK )) && continue
    raw="${ALL[$i]}"
    gz="$raw.gz"
    [[ -f "$raw" ]] || continue                 # handled by a concurrent rerun
    rawsz=$(stat -c %s "$raw")

    if [[ -f "$gz" ]]; then
        # Test before destroying. A valid .gz here means the raw is genuinely
        # redundant, which is the cheap outcome: drop the raw, keep the .gz.
        if pigz -t "$gz" 2>/dev/null; then
            rm -f "$raw"
            kept=$((kept+1)); bytes_before=$((bytes_before+rawsz))
            echo "  [keep-gz ] $gz valid, removed redundant raw ($(numfmt --to=iec $rawsz))"
            continue
        fi
        echo "  [truncated] $gz failed gzip -t, recompressing from raw"
        rm -f "$gz"
    else
        echo "  [orphan   ] $raw has no .gz, compressing"
    fi

    t0=$SECONDS
    # -k keeps the raw until the verify below passes, so a failure here is
    # recoverable rather than destructive.
    if pigz -k -f -p "$THREADS" "$raw" 2>/dev/null && [[ -f "$gz" ]] \
       && pigz -t "$gz" 2>/dev/null; then
        gzsz=$(stat -c %s "$gz")
        rm -f "$raw"
        recovered=$((recovered+1))
        bytes_before=$((bytes_before+rawsz)); bytes_after=$((bytes_after+gzsz))
        echo "  [ok       ] $gz  $(numfmt --to=iec $rawsz) -> $(numfmt --to=iec $gzsz)  in $((SECONDS-t0))s"
    else
        rm -f "$gz"
        failed=$((failed+1))
        echo "  [FAIL     ] $raw  compression or verify failed after $((SECONDS-t0))s, raw kept"
    fi
done

echo
echo "task $TASK done: recovered=$recovered kept-existing-gz=$kept failed=$failed"
echo "  raw removed:   $(numfmt --to=iec ${bytes_before:-0})"
echo "  gz written:    $(numfmt --to=iec ${bytes_after:-0})"
echo "  net reclaimed: $(numfmt --to=iec $(( bytes_before - bytes_after )))"
