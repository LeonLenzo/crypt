#!/bin/bash
#SBATCH --job-name=busco135
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-2
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/busco_%A_%a.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/busco_%A_%a.err

# BUSCO on the PREDICTED PROTEINS. This is the check predict_counts cannot do: the worklist
# is by definition the set with no annotation, so there is no truth to compare against, and
# containment against a reference is unavailable. A real fungal gene set recovers most of
# its lineage s BUSCO complement, and that needs no reference at all.
#
# -m proteins, not transcriptome. The old kraken_db_busco ran transcriptome mode on
# downloaded CDS; here vipsania already emitted proteins.fa, so score the gene models
# directly rather than re-deriving ORFs from the CDS.
#
# WHAT WOULD FALSIFY THE PREDICTIONS. 100% of these CDS start ATG, end on a stop codon and
# are length %3==0, but an ab-initio finder ALWAYS emits syntactically valid ORFs, so that
# rules out malformed output and nothing more. If completeness comes back low while the
# gene COUNT looked plausible, the models are ORFs in the right organism but in the wrong
# places, and predicted CDS would be competitor-only in a much stronger sense than policy
# already assumes.
#
# LINEAGE DBs GO ON /software, NOT /scratch. The previous copies were wiped by the 21-day
# purge, which is the documented failure where a purged BUSCO asset silently marks every
# assembly busco_error in ~90 s (2026-09-24). Same reason the container already lives there.

# WHY THIS EXISTS. The first full pass (job 53232778) returned NO_SUMMARY for 36 of 135,
# every one of them with "ERROR: Duplicate of sequence >...". That was not a BUSCO fault
# and not a prediction fault: two prediction arrays ran the same list concurrently and
# vipsania appended, so those 36 held two complete copies of every record. They are
# repaired by dedup_predictions.py (record-level, identity verified). This re-scores them.
#
#   LIST=data/redo36_busco.tsv sbatch research/predict/busco_redo.slurm
#

set -uo pipefail
module load singularity/3.11.4-nompi   # NOT on PATH by default; task 1 died with
                                        # "singularity: command not found" without it
R=/scratch/pawsey1168/llenzo/crypt
P=$R/research/predict
SIF=/software/projects/pawsey1168/llenzo/containers/busco_v5.7.0.sif
DB=$MYSOFTWARE/busco_downloads
OUT=$P/data/busco
mkdir -p "$DB/lineages" "$OUT" "$P/logs"

row=$(tail -n +2 "$P/${LIST:-data/all135_busco.tsv}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -z "$row" ] && { echo "no row ${SLURM_ARRAY_TASK_ID}"; exit 1; }
ACC=$(echo "$row" | cut -f1); TAXON=$(echo "$row" | cut -f2); LIN=$(echo "$row" | cut -f4)
FAA=$P/data/predicted/$ACC/proteins.fa
echo "START $(date -Is)  $ACC  $TAXON  lineage=$LIN"
[ -s "$FAA" ] || { echo "FAIL no proteins.fa at $FAA"; exit 1; }
echo "  proteins: $(grep -c "^>" "$FAA") sequences"

# Lineage download is serialised with a lock: five array tasks sharing $DB would otherwise
# race on the same tarball and leave a half-extracted directory that looks present.
exec 9>"$DB/.${LIN}.lock"; flock 9
if [ ! -s "$DB/lineages/$LIN/dataset.cfg" ] || [ -z "$(ls "$DB/lineages/$LIN/prfl" 2>/dev/null)" ]; then
  echo "  downloading $LIN"
  rm -rf "$DB/lineages/$LIN"
  curl --retry 3 --retry-delay 10 -fsSL -o "$DB/lineages/$LIN.tgz" \
       "https://busco-data.ezlab.org/v5/data/lineages/${LIN}.2024-01-08.tar.gz" \
    && tar -xzf "$DB/lineages/$LIN.tgz" -C "$DB/lineages" && rm -f "$DB/lineages/$LIN.tgz"
fi
flock -u 9; exec 9>&-
[ -s "$DB/lineages/$LIN/dataset.cfg" ] || { echo "FAIL lineage $LIN unavailable"; exit 1; }

cd "$OUT" || exit 1
singularity exec --bind /scratch:/scratch --bind /software:/software \
  --env HOME="$OUT" "$SIF" \
  busco -m proteins -i "$FAA" -l "$LIN" -o "$ACC" --download_path "$DB" \
        --offline -c "${SLURM_CPUS_PER_TASK:-8}" -f 2>&1 | tail -12

S=$(ls "$OUT/$ACC"/short_summary.specific.*.txt 2>/dev/null | head -1)
if [ -s "$S" ]; then
  line=$(grep -m1 "C:" "$S" | tr -d " \t")
  printf "BUSCO\t%s\t%s\t%s\t%s\n" "$ACC" "$TAXON" "$LIN" "$line"
else
  printf "BUSCO\t%s\t%s\t%s\tNO_SUMMARY\n" "$ACC" "$TAXON" "$LIN"
fi
echo "END $(date -Is)"
