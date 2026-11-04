#!/bin/bash
#SBATCH --job-name=fetch_reads
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=release/01_cohort/slurm/logs/fetch_%j.out
#SBATCH --error=release/01_cohort/slurm/logs/fetch_%j.err
#
# Fetch the reads ENA still owes the cohort, into release/01_cohort/raw/reads/<BioProject>/.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/01_cohort/slurm/fetch_reads.slurm
#
# PARTITION=copy, not work. This is pure network and filesystem I/O: ENA serves ready-made
# fastq.gz, so nothing is decompressed or recompressed and no core does arithmetic. The copy
# partition exists for exactly this and keeps a long transfer off the compute nodes; --mem=8G
# because the memory request sets the CPU bill and this needs almost none.
#
# WORKERS=24. Against EBI each stream gets its own ~0.7 MB/s (the limit is latency, not
# contention), so throughput scales with streams; 24 held ~17.5 MB/s. The gap is ~1 TB, so
# ~15 h, inside the 24 h wall. If EBI starts refusing connections late in a long run (it has
# before), the job can end short of finished: just resubmit. A completed file is never
# refetched (fetch_reads.py checks the disk, not the manifest), so the next job resumes.
set -uo pipefail
module load python/3.11.6
cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p release/01_cohort/slurm/logs
echo "host     $(hostname)"
echo "started  $(date -Is)"
SECONDS=0
python3 release/01_cohort/fetch_reads.py \
    --reads-dir release/01_cohort/raw/reads --fetch --workers "${WORKERS:-24}"
rc=$?
echo "exit     $rc after ${SECONDS}s"
echo "finished $(date -Is)"
exit $rc
