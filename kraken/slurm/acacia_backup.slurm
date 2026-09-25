#!/bin/bash
#SBATCH --job-name=acacia_backup
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/acacia_backup_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/acacia_backup_%j.err

# 2026-09-25. Replaces backup_dbv3.slurm, which looped 6x on a hard failure and wrote
# 14 MB of "argument of type NoneType is not iterable". That message is aws-cli 2.13.0
# crashing in s3errormsg.py while formatting a Ceph <Code>QuotaExceeded</Code> whose
# <Message> is empty. It was never SSL. The cause was the bucket sitting on the
# personal llenzo allocation (100 GiB, full) instead of the project allocation
# (pawsey1168, 1.1 TB, empty). Destination is now the project bucket.
#
# Order is deliberate: the files needed to CLASSIFY with db_v3 go first, then the
# files needed to EXTEND it, then the bulk that is regenerable.

module load aws-cli/2.13.0
P="--profile pawsey1168"          # endpoint_url now lives in the profile config
B=s3://pawsey1168-kraken-db
SRC=/scratch/pawsey1168/llenzo/crypt/kraken
V3=$SRC/build/data/db_v3

fail=0

# Fail fast rather than retrying a refusal: a 6-byte canary proves we can write.
echo probe > /tmp/acacia_probe.$$
if ! aws s3 cp /tmp/acacia_probe.$$ $B/.probe $P --only-show-errors 2>/dev/null; then
  echo "ABORT: cannot write to $B. Check the quota panel; aws-cli hides QuotaExceeded."
  rm -f /tmp/acacia_probe.$$; exit 1
fi
aws s3 rm $B/.probe $P --only-show-errors; rm -f /tmp/acacia_probe.$$

step () {  # step <label> <aws args...>
  local label=$1; shift
  echo "=== $label START $(date) ==="
  if "$@"; then echo "=== $label OK $(date) ==="
  else echo "=== $label FAILED $(date) ==="; fail=1; fi
}

# 1. Enough to classify with db_v3.
for f in opts.k2d taxo.k2d seqid2taxid.map hash.k2d; do
  step "db_v3/$f" aws s3 cp $V3/$f $B/db/build/db_v3/$f $P --only-show-errors
done
# 2. Enough to extend it.
step "db_v3/taxonomy" aws s3 sync $V3/taxonomy/ $B/db/build/db_v3/taxonomy/ $P --only-show-errors
# 3. The NCBI downloads: expensive to re-pull, and the seed for any db_v4.
step "cds/pathogen" aws s3 sync $SRC/search/data/cds/pathogen/ $B/db/search/cds/pathogen/ $P --only-show-errors
step "cds/host"     aws s3 sync $SRC/search/data/cds/host/     $B/db/search/cds/host/     $P --only-show-errors
# 4. Build intermediates, regenerable from the cds above.
step "db_v3/library" aws s3 sync $V3/library/ $B/db/build/db_v3/library/ $P --only-show-errors
# 5. Small and cheap.
step "busco"  aws s3 sync $SRC/utilities/busco/data/ $B/db/busco/ $P --only-show-errors
step "db_v2"  aws s3 sync $SRC/build/data/db_v2/ $B/db/build/db_v2/ $P --only-show-errors

echo "=== USAGE $(date) ==="
aws s3 ls $B --recursive --summarize $P 2>/dev/null | tail -2
[ $fail -eq 0 ] && echo "=== ALL STEPS OK ===" || echo "=== ONE OR MORE STEPS FAILED, see above ==="
exit $fail
