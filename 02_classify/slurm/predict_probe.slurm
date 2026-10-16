#!/bin/bash
#SBATCH --job-name=predict_probe
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_probe_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/predict_probe_%j.err

# STEP 0 of the predict evaluation: can vipsania even run here, and how fast on CPU?
# Everything downstream is moot if this fails, so it runs before the 7 GB genomic fetch.
#
# THE RISK, stated so the failure is recognisable: vipsania documents "TensorFlow with CUDA
# support", and Setonix has no CUDA. The GPU partition is AMD MI250X behind a ROCm 6.3.3
# TensorFlow container (quay.io/pawsey/tensorflow:2.17-rocm6.3.3). A CUDA-pinned wheel will
# not run on that hardware at all.
#
# WHY CPU AND WHY `work`: the gpu partition had 374 jobs PENDING against 8 idle nodes, while
# work had 44 idle. vipsania's docs say CPU annotation is "possible but slow". For 119 fungal
# genomes of ~35-60 Mbp, a tolerable CPU time beats queueing behind 374 jobs, and it removes
# the ROCm question entirely. So measure CPU first; only fight for a GPU if CPU is hopeless.
# That is what this job answers, and the number it prints is the one that decides.
#
#   sbatch 02_classify/slurm/predict_probe.slurm

set -uo pipefail          # NOT -e: a probe must report a failure, not vanish on it

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

VENV="$SCRATCH/venv/vipsania"
GDIR="$SCRATCH/crypt/02_classify/data/bronze/genomes"
OUT="$SCRATCH/crypt/02_classify/data/bronze/predicted/probe"
mkdir -p "$OUT"

say () { echo; echo "──── $* ────"; }

say "environment"
module load cray-python/3.12.12
echo "python: $(python3 --version 2>&1)"
echo "host:   $(hostname)   cpus=${SLURM_CPUS_PER_TASK:-?}"

# Every cache these libraries use defaults somewhere under $HOME, which a Setonix compute node
# cannot write. The first probe died on exactly that:
#   OSError: Cannot create the model directory /home/llenzo/.cache/vipsania/models
# so all of them are redirected to scratch. wandb comes in as a vipsania dependency and will
# block on a login prompt if it is left enabled, which on a batch node means a silent hang
# until the walltime expires rather than an error.
export VIPSANIA_CACHE="$SCRATCH/cache/vipsania"
export KERAS_HOME="$SCRATCH/cache/keras"
export XDG_CACHE_HOME="$SCRATCH/cache/xdg"
export HF_HOME="$SCRATCH/cache/hf"
export WANDB_MODE=offline
export WANDB_DIR="$SCRATCH/cache/wandb"
export TF_CPP_MIN_LOG_LEVEL=2
mkdir -p "$VIPSANIA_CACHE" "$KERAS_HOME" "$XDG_CACHE_HOME" "$HF_HOME" "$WANDB_DIR"
echo "VIPSANIA_CACHE=$VIPSANIA_CACHE"

say "venv + install"
if [ ! -x "$VENV/bin/python" ]; then
    python3 -m venv "$VENV" || { echo "PROBE FAIL: venv creation"; exit 1; }
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install --quiet --upgrade pip
echo "installing vipsania (this is the step that can fail on a CUDA pin)"
python -m pip install vipsania 2>&1 | tail -15
echo "pip exit: $?"

say "what actually got installed"
python - <<'PY'
import importlib, sys
for mod in ("vipsania", "tensorflow", "keras"):
    try:
        m = importlib.import_module(mod)
        print(f"  {mod:12s} {getattr(m, '__version__', '?')}")
    except Exception as exc:
        print(f"  {mod:12s} IMPORT FAILED: {type(exc).__name__}: {exc}")
try:
    import tensorflow as tf
    print("  TF devices:", [d.device_type for d in tf.config.list_physical_devices()])
    print("  built with cuda:", getattr(tf.test, "is_built_with_cuda", lambda: "?")())
except Exception as exc:
    print("  TF probe failed:", type(exc).__name__, exc)
PY

say "cli surface"
vipsania --help 2>&1 | head -25
echo "(--help exit: $?)"

say "pre-download the Fungi model, so it is not inside the timing"
# Fungi model id fh1kg88z, from the README clade table. `Fungi` also resolves, but the id is
# pinned here so a README edit cannot silently change which model we benchmarked.
MODEL=fh1kg88z
vipsania download "$MODEL" 2>&1 | tail -5
echo "(download exit: $?)"

say "fetch ONE genome to predict on"
python 02_classify/01_refs/predict/fetch_genomes.py \
    --genomes-dir "$GDIR" --arm annotated --limit 1 --workers 1 2>&1 | tail -8

FASTA=$(find "$GDIR" -name '*.fna' -size +1M | head -1)
if [ -z "$FASTA" ]; then echo "PROBE FAIL: no genome on disk to test"; exit 1; fi
echo "target: $FASTA  ($(du -h "$FASTA" | cut -f1))"
echo "contigs: $(grep -c '^>' "$FASTA")   bases: $(grep -v '^>' "$FASTA" | tr -d '\n' | wc -c)"

# TWO timed runs, because --finetune buys very little and may cost a great deal.
# The README reports Fungi locus F1 0.645 pretrained against 0.660 finetuned: +0.015 for a
# training pass over the target genome. On CPU that trade could easily be hours for nothing,
# so the pretrained run is timed FIRST and is the one that decides feasibility.
say "TIMED run 1 of 2: PRETRAINED only (no --finetune)"
# /usr/bin/time does not exist on a Setonix compute node (the first probe exited 127 on it).
# $SECONDS is a shell built-in and needs nothing installed. Peak memory comes from
# `sacct -j <id> --format=MaxRSS` afterwards.
t0=$SECONDS
vipsania annotate "$MODEL" "$FASTA" \
    -o "$OUT/pretrained.gff3" \
    --coding "$OUT/pretrained_coding.fa" \
    --protein "$OUT/pretrained_proteins.fa" 2>&1 | tail -25
rc=${PIPESTATUS[0]}
echo "WALL_PRETRAINED $((SECONDS - t0)) s   (exit: $rc)"

say "TIMED run 2 of 2: WITH --finetune"
t1=$SECONDS
vipsania annotate "$MODEL" "$FASTA" \
    -o "$OUT/finetuned.gff3" \
    --coding "$OUT/finetuned_coding.fa" \
    --protein "$OUT/finetuned_proteins.fa" --finetune 2>&1 | tail -25
rc=${PIPESTATUS[0]}
echo "WALL_FINETUNED $((SECONDS - t1)) s   (exit: $rc)"

say "result"
for f in "$OUT"/*_coding.fa "$OUT"/*.gff3; do
    if [ -s "$f" ]; then
        n=""; case "$f" in *.fa) n=", sequences: $(grep -c '^>' "$f")";; esac
        echo "  $(basename "$f"): $(du -h "$f" | cut -f1)$n"
    else
        echo "  $(basename "$f"): ABSENT OR EMPTY"
    fi
done
echo
echo "Compare against the TRUTH we already hold for this accession:"
echo "  real CDS lives under data/bronze/cds/pathogen/<accession>/ and its gene count is in"
echo "  ref_candidates.tsv (protein_coding_genes). compare.py does this properly."
echo "END $(date)"
