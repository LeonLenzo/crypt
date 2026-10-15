#!/bin/bash
#SBATCH --job-name=phordei_test
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=115G
#SBATCH --time=03:00:00
#SBATCH --output=04_kallisto/slurm/logs/phordei_test_%j.out
#SBATCH --error=04_kallisto/slurm/logs/phordei_test_%j.err

# Blind-spot experiment: does barley leaf rust (Puccinia hordei) surface when its gene-predicted
# CDS is added to the Hordeum index? Currently P. hordei is absent from db_v3 and its reads are
# misattributed to P. triticina (04_confirm). See crypt-align-crosstalk-measurement.
#
# Builds an augmented Hordeum index (name hordeum_phordei) with P. hordei added as a competitor,
# re-quants every barley run in the stratified set against it, and 04_confirm reads the result.
# Does NOT touch the validated hordeum_vulgare index or the baseline quant; compare the
# phordei quant's P. triticina/striiformis/hordei unique reads against the baseline calls.tsv.

set -uo pipefail
module load cray-python/3.11.7
module load spack/default
spack load kallisto@0.50.1
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 04_kallisto/slurm/logs

echo "=== 1. build augmented Hordeum index (P. hordei added) ==="
python3 04_kallisto/02_build/kallisto_build.py \
    --host "Hordeum vulgare" --name hordeum_phordei --force \
    --add "GCA_062009444.1:Puccinia hordei" --threads "${SLURM_CPUS_PER_TASK:-8}"

echo
echo "=== 2. re-quant the barley strata runs against it ==="
IDX=04_kallisto/02_build/data/idx/hordeum_phordei.idx
READS=03_kraken/03_select/data
OUT=04_kallisto/03_quant/data_phordei
python3 - "$IDX" "$READS" "$OUT" "${SLURM_CPUS_PER_TASK:-8}" <<'PY'
import csv, subprocess, sys
from pathlib import Path
idx, reads_root, out_root, threads = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
runs = {}
for r in csv.DictReader(open("04_kallisto/01_select/data/strata.tsv"), delimiter="\t"):
    if r["host"] == "Hordeum vulgare" and r["run"] not in runs:
        runs[r["run"]] = (r["layout"], [reads_root / p for p in r["read_files"].split(",") if p])
print(f"{len(runs)} barley runs to re-quant", flush=True)
for run, (layout, reads) in runs.items():
    od = out_root / run; od.mkdir(parents=True, exist_ok=True)
    cmd = ["kallisto", "quant", "-i", idx, "-o", str(od), "-t", threads]
    if layout != "paired":
        cmd += ["--single", "-l", "200", "-s", "20"]
    cmd += [str(p) for p in reads]
    subprocess.run(cmd, check=True)
    bd = od / "bus"; bd.mkdir(exist_ok=True)
    cmd = ["kallisto", "bus", "-i", idx, "-o", str(bd), "-x", "Bulk", "-t", threads, "--num"]
    if layout == "paired":
        cmd.append("--paired")
    cmd += [str(p) for p in reads]
    subprocess.run(cmd, check=True)
    print(f"  done {run}", flush=True)
PY
echo "exit $?"
