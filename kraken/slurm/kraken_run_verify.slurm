#!/bin/bash
#SBATCH --job-name=krsel_verify
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-15
#SBATCH --output=kraken/slurm/logs/krsel_verify_%A_%a.out
#SBATCH --error=kraken/slurm/logs/krsel_verify_%A_%a.err

# READ ONLY. Deletes nothing, writes nothing into the reads directory.
#
# Checks every downloaded read file on two independent axes:
#
#   integrity    is the gzip stream complete? (a truncated member fails to
#                decompress). This catches the damage from the old flat 1800s
#                gzip timeout described in kraken_run_select.py's _compress
#                docstring.
#   completeness does the file contain as many reads as SRA says the run has?
#                Integrity alone cannot answer this: gzip will happily produce
#                a perfectly valid stream wrapping a SHORT read set, so a
#                `gzip -t` pass is not evidence that the download finished.
#
# Expected read counts come from NCBI STAT's own tax_totals.total, already held
# in stat/stat_fetch/data/stat_cache.jsonl for 607k runs (0.3% of which
# STAT subsampled, where analysed < total; total is still the true spot count).
# Extracted to spots.tsv and shipped alongside this job.
#
# Both axes are measured in ONE decompression pass: `pigz -dc | wc -l` validates
# the stream and counts records together, so completeness costs nothing extra
# over the integrity check.
#
# Paired runs: each mate file should hold `spots` records, so expected records
# per file == spots regardless of mate.
#
# Output: kraken/run/kraken_run_select/data/verify_%a.tsv
#   file  kind(gz|raw)  bytes  twin(yes|no)  records  expected  ratio  integrity  verdict
#
#   verdict  COMPLETE    integrity ok and records within 1% of expected
#            SHORT       integrity ok but records well below expected
#            OVER        more records than expected (investigate, not damage)
#            CORRUPT     stream failed to decompress
#            NO_EXPECTED no spot count available for this run
#
# Run from crypt/: sbatch kraken/slurm/kraken_run_verify.slurm
# Collate:         cat kraken/run/kraken_run_select/data/verify_*.tsv > verify.tsv

set -uo pipefail

READS=/scratch/pawsey1168/llenzo/crypt/kraken/run/kraken_run_select/data/reads
OUT=/scratch/pawsey1168/llenzo/crypt/kraken/run/kraken_run_select/data
SPOTS=/scratch/pawsey1168/llenzo/crypt/kraken/run/kraken_run_select/data/spots.tsv
THREADS=${SLURM_CPUS_PER_TASK:-8}
TASK=${SLURM_ARRAY_TASK_ID:-0}
NTASK=${SLURM_ARRAY_TASK_COUNT:-1}

cd "$READS" || { echo "FATAL: cannot cd to $READS"; exit 1; }
command -v pigz >/dev/null || { echo "FATAL: pigz not found"; exit 1; }
[[ -f "$SPOTS" ]] || { echo "FATAL: $SPOTS missing"; exit 1; }

TSV="$OUT/verify_${TASK}.tsv"
: > "$TSV"

# Every read file, gz and raw alike. Raw files are the ones about to be
# recompressed, so their completeness matters most.
mapfile -t ALL < <( { ls *.fastq.gz 2>/dev/null; ls *.fastq 2>/dev/null | grep -v '\.gz$'; } | sort -u )
echo "task $TASK/$NTASK  files: ${#ALL[@]}  threads=$THREADS"

# spots.tsv covers 607k runs; a per-file awk scan over that would be ~2.5e9
# comparisons across this task's share. Narrow it to the runs actually on disk
# once, then hold it in an associative array for O(1) lookup in the loop.
MYRUNS=$(mktemp); SMALL=$(mktemp)
printf '%s\n' "${ALL[@]}" | sed -E 's/\.fastq(\.gz)?$//; s/_[12]$//' | sort -u > "$MYRUNS"
awk -F'\t' 'NR==FNR{want[$1];next} ($1 in want){print $1"\t"$2}' "$MYRUNS" "$SPOTS" > "$SMALL"
declare -A EXP
while IFS=$'\t' read -r _r _t; do [[ -n "$_r" ]] && EXP["$_r"]="$_t"; done < "$SMALL"
echo "  distinct runs on disk: $(wc -l < "$MYRUNS")   with spot counts: ${#EXP[@]}"
rm -f "$MYRUNS" "$SMALL"

for i in "${!ALL[@]}"; do
    (( i % NTASK != TASK )) && continue
    f="${ALL[$i]}"
    [[ -e "$f" ]] || continue

    if [[ "$f" == *.gz ]]; then kind=gz; twinf="${f%.gz}"; else kind=raw; twinf="$f.gz"; fi
    [[ -e "$twinf" ]] && twin=yes || twin=no
    bytes=$(stat -c %s "$f" 2>/dev/null || echo 0)

    # Run accession: strip the mate suffix and the extension(s).
    run=$(echo "$f" | sed -E 's/\.fastq(\.gz)?$//; s/_[12]$//')
    exp=${EXP["$run"]:-0}

    # One pass: decompression validates the stream, wc counts the records.
    if [[ "$kind" == gz ]]; then
        lines=$(pigz -dc "$f" 2>/dev/null | wc -l)
        st=${PIPESTATUS[0]}
    else
        lines=$(wc -l < "$f" 2>/dev/null)
        st=0
    fi

    if (( st != 0 )); then
        integrity=FAIL; recs=0; ratio=0; verdict=CORRUPT
    else
        integrity=OK
        recs=$(( lines / 4 ))
        if (( exp == 0 )); then
            ratio=0; verdict=NO_EXPECTED
        else
            ratio=$(awk -v a="$recs" -v b="$exp" 'BEGIN{printf "%.4f", a/b}')
            verdict=$(awk -v r="$ratio" 'BEGIN{ if (r>=0.99 && r<=1.01) print "COMPLETE"; else if (r>1.01) print "OVER"; else print "SHORT" }')
        fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$f" "$kind" "$bytes" "$twin" "$recs" "$exp" "$ratio" "$integrity" "$verdict" >> "$TSV"
done

echo "task $TASK done -> $TSV"
awk -F'\t' '{c[$9]++} END{for (k in c) printf "  %-12s %d\n", k, c[k]}' "$TSV"
