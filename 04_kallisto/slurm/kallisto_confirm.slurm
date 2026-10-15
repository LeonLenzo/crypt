#!/bin/bash
#SBATCH --job-name=kallisto_confirm
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=04_kallisto/slurm/logs/kallisto_confirm_%j.out
#SBATCH --error=04_kallisto/slurm/logs/kallisto_confirm_%j.err

# Turn the quant output into the three per-detection facts -> 04_confirm/data/calls.tsv.
#
# Single-threaded: it re-reads each run's matrix.ec and output.bus (matrix.ec is per-run, see
# kallisto_confirm.py) and parses the BUS binary directly, ~5-35s per run, so ~40 min for 174
# runs. transcripts.txt (the million-line file) is cached per host. 2h walltime is generous
# headroom; a 30-min limit timed out at run 136/174.
#
# Progress prints per run to 04_confirm/logs/, NOT to the SLURM .out (step_log redirects stdout).
#
# Run from crypt/: sbatch 04_kallisto/slurm/kallisto_confirm.slurm

set -uo pipefail
module load cray-python/3.11.7
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 04_kallisto/slurm/logs

echo "host : $(hostname)"
echo "job  : ${SLURM_JOB_ID:-none}"
echo

python3 04_kallisto/04_confirm/kallisto_confirm.py --all
echo "exit $?"
