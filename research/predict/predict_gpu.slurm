#!/bin/bash
#SBATCH --job-name=vipsania_gpu
#SBATCH --account=fl3-gpu
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --array=1-10%4
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/vip_%A_%a.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/vip_%A_%a.err

# vipsania gene prediction, one genome per array task, on Setonix GPUs.
#
#   LIST=data/pilot10.tsv sbatch research/predict/predict_gpu.slurm
#
# GPU, NOT CPU, and the trial settled it (jobs 51036963 / 51037118, 16 Oct 2026). On the same
# genome, GPU took 448 s against the CPU route s 2,447 s and produced byte-comparable output:
# 13,634 CDS both ways. Billing moves from pawsey1168 (CPU, 400,000 SU, 38% used) to fl3-gpu
# (87,500 SU, essentially untouched), and the per-genome cost drops from ~51 SU to ~8.
#
# DO NOT PASS --rocm. That was the ONLY reason the first trial died. It injects the host s
# ROCm libraries into /.singularity.d/libs, shadowing the ROCm 6.3.3 the Pawsey container
# already bundles, and the injected rocBLAS then cannot find its Tensile kernels:
#   rocBLAS error: Cannot read /.singularity.d/libs/rocblas/library/TensileLibrary.dat
#                  for GPU arch : gfx90a
# It aborted a plain tf.matmul, so it was never a vipsania problem. The container s own ROCm
# stack is complete and needs only device access, which singularity gives by default.
#
# NO --cpus-per-task. Setonix refuses it on a GPU allocation (each GPU allocates 8 cores).
#
# PARTITION=gpu, not gpu-dev. gpu-dev allows MaxSubmitJobs=4 / MaxJobs=1 for this account,
# so a 10-task array is rejected with AssocMaxSubmitJobLimit. gpu allows 1024 / 64 and its
# 1-day limit also covers the largest genomes.
#
# WALLTIME 4 h. Runtime scales with genome size: 36 Mb took 448 s on GPU, so the largest in
# the worklist (Puccinia hordei, 143 Mb) should be ~30 min. 4 h is ample headroom.

set -uo pipefail
R=/scratch/pawsey1168/llenzo/crypt
LIST="${LIST:-data/pilot10.tsv}"
G=$R/02_classify/data/bronze/genomes
OUT=$R/research/predict/data/predicted
MODEL=fh1kg88z
cd "$R/research/predict" || exit 1
mkdir -p logs "$OUT"

# row N of the list, skipping the header
row=$(tail -n +2 "$LIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -z "$row" ] && { echo "no row ${SLURM_ARRAY_TASK_ID} in $LIST"; exit 1; }
ACC=$(echo "$row" | cut -f1); TAXON=$(echo "$row" | cut -f2); MB=$(echo "$row" | cut -f3)
echo "START $(date -Is)  task=${SLURM_ARRAY_TASK_ID}  $ACC  $TAXON  ${MB} Mb"

FASTA=$(find "$G/$ACC" -name "*_genomic.fna" 2>/dev/null | head -1)
[ -s "$FASTA" ] || { echo "FAIL no genome fasta under $G/$ACC"; exit 1; }
echo "  fasta $(basename "$FASTA")  $(du -h "$FASTA" | cut -f1)"

D="$OUT/$ACC"; mkdir -p "$D"
if [ -s "$D/coding.fa" ]; then echo "  already predicted, skipping"; exit 0; fi

# pawseyenv FIRST. The October probe worked without it, but the site default has moved on
# and the bare module now fails with "exist but cannot be loaded as requested", which
# surfaces as exit 127 (tensorflow-exec: command not found) at 0 s per task.
module load pawseyenv/2025.08
module load tensorflow/2.17-rocm6.3.3
V=$MYSOFTWARE/venv/vipsania-rocm
# $HOME is at 100% quota, so every cache must live on /software
export VIPSANIA_CACHE=$MYSOFTWARE/cache/vipsania
export KERAS_HOME=$MYSOFTWARE/cache/keras
export XDG_CACHE_HOME=$MYSOFTWARE/cache/xdg
export HF_HOME=$MYSOFTWARE/cache/hf
export WANDB_MODE=offline
export WANDB_DIR=/scratch/pawsey1168/llenzo/cache/wandb
export TF_CPP_MIN_LOG_LEVEL=2

t0=$SECONDS
tensorflow-exec "$V/bin/vipsania" annotate "$MODEL" "$FASTA" \
    -o "$D/genes.gff3" --coding "$D/coding.fa" --protein "$D/proteins.fa" -B 8 2>&1 | tail -8
rc=${PIPESTATUS[0]}
wall=$((SECONDS - t0))

if [ -s "$D/coding.fa" ]; then
  n=$(grep -c "^>" "$D/coding.fa"); bp=$(grep -v "^>" "$D/coding.fa" | tr -d "\n" | wc -c)
  printf "RESULT\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n" "$ACC" "$TAXON" "$MB" "$wall" "$rc" "$n" "$bp"
  echo "  coding.fa $n CDS, $bp bp in ${wall}s (exit $rc)"
else
  printf "RESULT\t%s\t%s\t%s\t%d\t%d\t0\t0\n" "$ACC" "$TAXON" "$MB" "$wall" "$rc"
  echo "  FAIL no coding.fa after ${wall}s (exit $rc)"
fi
echo "END $(date -Is)"
[ -s "$D/coding.fa" ]
