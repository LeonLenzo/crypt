#!/bin/bash
#SBATCH --job-name=vipsania
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --array=1-120%40
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/vipsania_%A_%a.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/vipsania_%A_%a.err

# One genome per array task. Measured on the probe (job 51033066): a 37 Mbp chromosome-level
# fungal genome took 2,447 s (40.8 min) and peaked at 67 GB RSS, so --mem=120G is headroom
# over the measured figure rather than a guess — the probe itself was OOM-killed at 64G when
# it went on to the finetuning arm.
#
# MEMORY 48G, AND THAT IS A BILLING DECISION AS MUCH AS A CORRECTNESS ONE.
# Setonix couples memory to cores at ~1.9 GB per core, so the memory request — not
# --cpus-per-task — decides how many CPUs are ALLOCATED AND BILLED. The first attempt asked
# for 32 cores and 120G and was billed 134 cores per task, which over 2,138 genomes would be
# ~196,000 core-hours, 49% of the project's 400,000 SU allocation, for a benchmark.
#
# Measured instead of assumed. The probe's 67 GB MaxRSS was the FINETUNING arm, which is the
# part that OOM-killed it; annotation alone sits at a flat 30 GB. Sampled across five live
# tasks at 6-10 min elapsed: 29.4, 29.7, 30.0, 30.1, 30.2 GB. It is the model plus batch
# working set and does not grow with genome size, because vipsania works in windows of the
# context length.
#
# 48G therefore carries ~60% headroom over the real figure and bills ~50 cores, against the
# ~48 that vipsania actually drives (observed 4,827% CPU). Going lower buys nothing: below
# CPU saturation the runtime stretches and core-hours stay flat. This cuts the benchmark from
# ~11,000 to ~4,100 core-hours and the full 2,138-genome run from 49% of the allocation to
# ~18%.
#
# WALLTIME 10 h, not 4. Runtime scales with genome SIZE and the sample is not uniform: median
# 38 MB but the rusts and ectomycorrhizal fungi are genuinely large — Cenococcum geophilum
# 180 MB, Puccinia hordei 163 MB, P. triticina 122 MB. At the probe's rate a 180 MB genome is
# ~3.3 h, which a 4 h limit would clear only barely. MEMORY does not need the same treatment:
# vipsania works in windows of the context length, so peak RSS tracks batch x T, not total
# genome size.
#
# Plant hosts are excluded at the SAMPLE stage, not here. An unfiltered Arm A draw had pulled
# Triticum aestivum GCA_040256815.2 — 14.69 GB, ~270 h of inference, and scored against the
# Fungi model it would have been meaningless anyway.
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
