#!/bin/bash
#SBATCH --job-name=vipsania_gpu
#SBATCH --account=fl3-gpu
#SBATCH --partition=gpu-dev
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/vipsania_gpu_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/vipsania_gpu_%j.err

# Can vipsania use Setonix's AMD GPUs, and is it cheaper than CPU?
#
# WHY THIS MATTERS FOR THE BUDGET, not just the clock. GPU work bills to fl3-gpu (87,500 SU,
# untouched) and NOT to pawsey1168, which has no gpu association at all. The CPU route puts
# the full 2,138-genome run at ~78,000 core-hours, ~20% of the pawsey1168 allocation. If this
# works, that cost moves to a different and idle budget.
#
# WHY IT MIGHT NOT WORK. vipsania depends on hidten[tensorflow], which pins
# tensorflow[and-cuda], and Setonix has no CUDA — the GPUs are AMD MI250X. The pip wheel
# installed for the CPU route reports is_cuda_build=True and sees no device here. The way
# round it: Pawsey's ROCm container (quay.io/pawsey/tensorflow:2.17-rocm6.3.3) reports
# is_rocm_build=True with TF 2.17.0, which satisfies vipsania's >=2.17,<2.20, and a venv built
# with --system-site-packages against it lets vipsania use ROCm TF instead of its own.
# Installed with --no-deps so pip cannot drag CUDA TF back in.
#
# Residual risk is hidten's custom HMM and parallel scan: if those lower to ops ROCm TF does
# not implement, this fails regardless of versions. The existence of a --nojit flag upstream
# suggests JIT trouble was already anticipated, so that is the first fallback to try.
#
# NO --cpus-per-task. Setonix refuses it on a GPU allocation:
#   lua: cli_filter: cannot explicitly request CPU resources for GPU allocation;
#        each allocated GPU allocates 8 cores
# so the core count follows the GPU count and must not be stated.
#
# gpu-dev, not gpu: gpu had 351 jobs PENDING while gpu-dev had 3 running and none queued, and
# 4 h is ample for one genome.
#
#   sbatch 02_classify/slurm/predict_probe_gpu.slurm

set -uo pipefail          # no -e: a probe must report its failure

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

say () { echo; echo "──── $* ────"; }

module load tensorflow/2.17-rocm6.3.3
V="$MYSOFTWARE/venv/vipsania-rocm"

# Singularity needs the ROCm devices bound in or the container sees no GPU. The wrapper passes
# $SINGULARITY_COMMAND_OPTS through to `singularity exec`.
export SINGULARITY_COMMAND_OPTS="--rocm"

# $HOME is at 100% quota, so every cache must live elsewhere; weights go on /software.
export VIPSANIA_CACHE="$MYSOFTWARE/cache/vipsania"
export KERAS_HOME="$MYSOFTWARE/cache/keras"
export XDG_CACHE_HOME="$MYSOFTWARE/cache/xdg"
export HF_HOME="$MYSOFTWARE/cache/hf"
export WANDB_MODE=offline
export WANDB_DIR="$SCRATCH/cache/wandb"
export TF_CPP_MIN_LOG_LEVEL=2

say "node and GPU"
echo "host: $(hostname)  cpus=${SLURM_CPUS_PER_TASK:-?}  gpus=${SLURM_GPUS_PER_NODE:-?}"
rocm-smi --showproductname 2>/dev/null | head -8 || echo "(rocm-smi not on PATH outside the container)"

say "does TensorFlow SEE the GPU?"
cat > /tmp/gputest.$$.py <<'PY'
import tensorflow as tf
bi = tf.sysconfig.get_build_info()
print("  tf", tf.__version__, " rocm_build:", bi.get("is_rocm_build"), " cuda_build:", bi.get("is_cuda_build"))
devs = tf.config.list_physical_devices()
print("  devices:", [(d.device_type, d.name) for d in devs])
gpus = tf.config.list_physical_devices("GPU")
print("  GPU count:", len(gpus))
if gpus:
    with tf.device("/GPU:0"):
        a = tf.random.normal([2048, 2048]); b = tf.random.normal([2048, 2048])
        c = tf.matmul(a, b)
    print("  matmul on GPU OK, result sum finite:", bool(tf.math.is_finite(tf.reduce_sum(c))))
else:
    print("  NO GPU VISIBLE — the rest of this probe measures nothing new")
PY
tensorflow-exec "$V/bin/python" /tmp/gputest.$$.py 2>&1 | grep -vE "^(W0|I0|E0|2026-)" | tail -12
rm -f /tmp/gputest.$$.py

say "TIMED annotate on the same genome the CPU probe used"
# P. nodorum GCA_016801405.1, 37.0 Mbp, 23 contigs. CPU baseline: 2,447 s (40.8 min) at 30 GB.
MODEL=fh1kg88z
FASTA="$SCRATCH/crypt/02_classify/data/bronze/genomes/GCA_016801405.1/ncbi_dataset/data/GCA_016801405.1/GCA_016801405.1_ASM1680140v1_genomic.fna"
OUT="$SCRATCH/crypt/02_classify/data/bronze/predicted/probe_gpu"
mkdir -p "$OUT"
if [ ! -s "$FASTA" ]; then echo "PROBE FAIL: genome missing at $FASTA"; exit 1; fi
echo "target: $(basename "$FASTA") ($(du -h "$FASTA" | cut -f1))"

# -B 8 matches the CPU run so the timings are comparable. On a GPU vipsania would normally
# auto-size this from nvidia-smi, which does not exist here either.
t0=$SECONDS
tensorflow-exec "$V/bin/vipsania" annotate "$MODEL" "$FASTA" \
    -o "$OUT/genes.gff3" --coding "$OUT/coding.fa" --protein "$OUT/proteins.fa" \
    -B 8 2>&1 | tail -20
rc=${PIPESTATUS[0]}
echo "WALL_GPU $((SECONDS - t0)) s   (exit: $rc)   CPU baseline was 2447 s"

say "result"
if [ -s "$OUT/coding.fa" ]; then
    echo "  coding.fa: $(du -h "$OUT/coding.fa" | cut -f1), $(grep -c '^>' "$OUT/coding.fa") sequences"
    echo "  CPU run gave 13,634 sequences / 19.29 Mb — these should match closely"
else
    echo "  coding.fa ABSENT: GPU route did not produce output"
fi
echo "END $(date)"
