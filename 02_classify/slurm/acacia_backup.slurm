#!/bin/bash
#SBATCH --job-name=acacia_backup_v4
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/acacia_backup_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/02_classify/slurm/logs/acacia_backup_%j.err

# Back up db_v4's CDS tree to Acacia's PROJECT bucket. ~2,518 of the 4,271 accessions were
# downloaded fresh from NCBI on 2026-10-16 and exist nowhere else; /scratch runs a 21-day
# purge, and purge removes FILES and leaves DIRECTORIES, which is how the previous tree came
# to show 1,997 accession directories holding only 1,162 FASTAs.
#
# `aws s3 sync` skips what is already there, so the ~1,753 reused accessions cost nothing.
#
# Destination is the PROJECT bucket on the 1.1 TB allocation, NOT
# `pawsey1168-llenzo-kraken-db`, which is the personal 100 GiB allocation and was full at
# 115.6 GB. Getting that wrong produces `argument of type 'NoneType' is not iterable` —
# aws-cli 2.13 crashing while formatting a Ceph QuotaExceeded whose <Message> is empty. It is
# never SSL and never the key.
#
# NOTE: the FASTAs in this tree are UNTAGGED (pristine NCBI headers). build_db.py writes its
# kraken:taxid| tags to a staging copy and deletes it, deliberately leaving bronze immutable.
# The 2026-09-25 backup did not do this, so the older objects in this prefix ARE tagged with
# db_v3 taxids — see the build_db.py docstring.
#
#   sbatch 02_classify/slurm/acacia_backup.slurm

set -euo pipefail

SCRATCH="/scratch/pawsey1168/llenzo"
REPO="$SCRATCH/crypt"
mkdir -p "$REPO/02_classify/slurm/logs"
cd "$REPO"

module load aws-cli/2.13.0

P="--profile pawsey1168"
B="s3://pawsey1168-kraken-db"
SRC="$SCRATCH/crypt/02_classify/data/bronze/cds/pathogen"

# A 6-byte canary proves we can write, so a refusal fails here rather than 2,518 objects in.
echo probe > /tmp/acacia_probe.$$
if ! aws s3 cp /tmp/acacia_probe.$$ $B/.probe $P --only-show-errors 2>/dev/null; then
    echo "ABORT: cannot write to $B. Check the PROJECT allocation quota in the portal's"
    echo "       S3 Bucket Viewer; aws-cli hides QuotaExceeded behind a NoneType crash."
    rm -f /tmp/acacia_probe.$$; exit 1
fi
aws s3 rm $B/.probe $P --only-show-errors; rm -f /tmp/acacia_probe.$$

echo "START $(date)"
echo "source : $SRC"
echo "local  : $(find "$SRC" -name '*.fna' -size +100k | wc -l) fasta, $(du -sh "$SRC" | cut -f1)"
echo "remote before: $(aws s3 ls $B/db/search/cds/pathogen/ --recursive $P 2>/dev/null | grep -c 'cds_from_genomic.fna' || true) cds_from_genomic.fna"

aws s3 sync "$SRC/" "$B/db/search/cds/pathogen/" $P --only-show-errors

echo "remote after : $(aws s3 ls $B/db/search/cds/pathogen/ --recursive $P 2>/dev/null | grep -c 'cds_from_genomic.fna' || true) cds_from_genomic.fna"
echo "END $(date)"
