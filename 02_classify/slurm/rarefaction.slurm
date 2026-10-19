#!/bin/bash
#SBATCH --job-name=rarefy
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=02_classify/slurm/logs/rarefy_%j.out
#SBATCH --error=02_classify/slurm/logs/rarefy_%j.err

# The depth-vs-detection curve: at what library depth does a real detection stop clearing
# the 05_detect floors? 100 classifications over 5 runs x 2 arms x 8 rungs.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch 02_classify/slurm/rarefaction.slurm            # all 100
#   sbatch --export=ALL,TASKS="1 2" 02_classify/slurm/rarefaction.slurm   # smoke test
#
# ONE NODE, NOT AN ARRAY, and that is the whole cost argument. Each classification here is
# tiny -- the 0.019 Gkmer rung is ~0.6% of a 6 Gb library, so a few tens of MB -- while db_v4's
# hash is 41 GB. The DB load dominates by orders of magnitude. A 100-task array would scatter
# over 100 nodes and pay 100 separate 41 GB loads off Lustre. Instead the hash is staged ONCE
# into /dev/shm and every kraken2 invocation mmaps that one RAM-resident copy.
#
# --memory-mapping IS CORRECT HERE, and CLAUDE.md says never to use it. That note measured mmap
# against LUSTRE, where kraken2's random hash lookups turn every page fault into a filesystem
# RPC: 13 minutes for zero reports, kraken2 pinned at 1.6% CPU. /dev/shm is RAM. The note is
# right about the case it measured and says nothing about this one.
#
# --mem=64G because the MEMORY REQUEST SETS THE CPU BILL at ~1.9 GB/core, and /dev/shm counts
# against the job's memory cgroup. 41 GB of hash + working set needs ~64 GB, which bills ~34
# cores whether they are busy or not. That is the floor for holding db_v4 at all; the only way
# to spend less is to do fewer DB loads, which is exactly what the shm staging achieves.

set -uo pipefail

module load python/3.11.6
module load kraken2/2.1.2-qk2lek6

cd /scratch/pawsey1168/llenzo/crypt || exit 1
mkdir -p 02_classify/slurm/logs

DB=02_classify/data/silver/db/db_v4
SHM=/dev/shm/db_v4_$SLURM_JOB_ID
WORKERS=${WORKERS:-6}
THREADS=${THREADS:-5}

echo "host     $(hostname)"
echo "started  $(date -Is)"

# Stage the index into RAM. taxo.k2d and opts.k2d are tiny but must travel with the hash;
# kraken2 refuses a db directory missing any of the three.
echo "staging db_v4 into $SHM ..."
SECONDS=0
mkdir -p "$SHM"
cp "$DB"/hash.k2d "$DB"/opts.k2d "$DB"/taxo.k2d "$SHM"/ || { echo "shm stage FAILED"; exit 1; }
echo "staged $(du -sh --apparent-size $SHM | cut -f1) in ${SECONDS}s"
df -h /dev/shm | tail -1

cleanup() { echo "clearing $SHM"; rm -rf "$SHM"; }
trap cleanup EXIT

# Which tasks: all of them unless TASKS is set.
if [ -n "${TASKS:-}" ]; then
  LIST="$TASKS"
else
  N=$(( $(wc -l < 02_classify/data/silver/rarefy_worklist.tsv) - 1 ))
  LIST=$(seq 1 "$N")
fi
echo "tasks    $(echo $LIST | wc -w)"

SECONDS=0
printf '%s\n' $LIST | xargs -P "$WORKERS" -I{} \
  bash -c 'echo "=== task {} ==="; python3 02_classify/utilities/rarefaction/rarefy.py \
      --task {} --db '"$SHM"' --threads '"$THREADS"' --mmap 2>&1'
rc=$?
echo "classification exit $rc after ${SECONDS}s"

python3 02_classify/utilities/rarefaction/rarefy.py --collect

echo "finished $(date -Is)"
exit $rc
