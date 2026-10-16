#!/bin/bash
#SBATCH --job-name=vipsania
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --array=1-119%40
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/vipsania_%A_%a.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/vipsania_%A_%a.err

# One genome per array task. Measured on the probe (job 51033066): a 37 Mbp chromosome-level
# fungal genome took 2,447 s (40.8 min) and peaked at 67 GB RSS, so --mem=120G is headroom
# over the measured figure rather than a guess — the probe itself was OOM-killed at 64G when
# it went on to the finetuning arm.
#
# NO --finetune, deliberately. It costs a training pass per genome and more memory than the
# annotation, for Fungi locus F1 0.660 against 0.645 pretrained. On our actual metric the
# pretrained model already reached recall 0.936 / precision 0.973 against the true CDS of
# P. nodorum, so there is nothing for finetuning to buy that justifies 119 training passes.
#
# -B 8 is the batch size the probe validated at this memory. -p is NOT passed: the HMM's
# parallel factor must divide the context length and vipsania selects the largest valid
# divisor itself (~400 for T=200000), so passing it can only break or reduce it.
#
# %40 throttles to 40 concurrent tasks; `work` had ~44 idle nodes and two tasks fit a node at
# 120G, so this fills what is free without monopolising the partition.
#
#   sbatch 02_classify/slurm/predict_vipsania.slurm

set -uo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load cray-python/3.12.12
# shellcheck disable=SC1091
source "$MYSOFTWARE/venv/vipsania/bin/activate"

# $HOME is at 100% of its 1 GiB quota, so every cache must be redirected or the model load
# fails with "Cannot create the model directory". Model weights sit on /software and survive
# the scratch purge.
export VIPSANIA_CACHE="$MYSOFTWARE/cache/vipsania"
export KERAS_HOME="$MYSOFTWARE/cache/keras"
export XDG_CACHE_HOME="$MYSOFTWARE/cache/xdg"
export HF_HOME="$MYSOFTWARE/cache/hf"
export WANDB_MODE=offline
export WANDB_DIR="$SCRATCH/cache/wandb"
export TF_CPP_MIN_LOG_LEVEL=2

MODEL=fh1kg88z
BATCH="${BATCH:-8}"
MANIFEST="$REPO/02_classify/data/silver/predict_genomes.tsv"
OUTROOT="$SCRATCH/crypt/02_classify/data/bronze/predicted/vipsania_pretrained"

# Row $SLURM_ARRAY_TASK_ID of the manifest, skipping the header. Keyed off the FETCH manifest,
# not the sample, so a task can never be handed an accession whose genome is not on disk.
LINE=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$MANIFEST")
if [ -z "$LINE" ]; then echo "no manifest row $SLURM_ARRAY_TASK_ID; nothing to do"; exit 0; fi
ACC=$(echo "$LINE" | cut -f2)
REL=$(echo "$LINE" | cut -f9)
NAME=$(echo "$LINE" | cut -f4)
FASTA="$SCRATCH/crypt/02_classify/data/bronze/genomes/$REL"

if [ ! -s "$FASTA" ]; then echo "ABORT: no genome at $FASTA"; exit 1; fi
OUT="$OUTROOT/$ACC"
mkdir -p "$OUT"

# Resumable: a completed task is skipped, so the array can be resubmitted after a timeout
# without redoing finished genomes.
if [ -s "$OUT/coding.fa" ]; then
    echo "already done: $ACC ($NAME)  $(grep -c '^>' "$OUT/coding.fa") sequences"
    exit 0
fi

echo "task $SLURM_ARRAY_TASK_ID  $ACC  $NAME"
echo "  genome: $FASTA ($(du -h "$FASTA" | cut -f1))"
t0=$SECONDS
vipsania annotate "$MODEL" "$FASTA" \
    -o "$OUT/genes.gff3" \
    --coding "$OUT/coding.fa" \
    --protein "$OUT/proteins.fa" \
    -B "$BATCH" 2>&1 | tail -12
rc=${PIPESTATUS[0]}
echo "  WALL $((SECONDS - t0)) s  (exit: $rc)"
if [ -s "$OUT/coding.fa" ]; then
    echo "  predicted: $(grep -c '^>' "$OUT/coding.fa") sequences"
else
    echo "  FAILED: no coding.fa written"; exit 1
fi
