#!/bin/bash
#SBATCH --job-name=backup_db
#SBATCH --account=pawsey1168
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/release/02_reference/scripts/slurm/logs/backup_db_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/release/02_reference/scripts/slurm/logs/backup_db_%j.err
#
# Back the built Kraken2 databases up to the Acacia TEAM bucket.
#
#   cd /scratch/pawsey1168/llenzo/crypt
#   sbatch release/02_reference/scripts/slurm/backup_db_acacia.slurm
#
# WHAT GOES UP, AND WHAT DOES NOT.
#
#   db1p  EVERYTHING, library included. It is the live database (7h21m of build in job
#         53334256) and the tree db2 will extend, so its library has to survive for the
#         extension to be possible without a full rebuild.
#   db0   k2d files + seqid2taxid.map + taxonomy/ only.
#   db1   same.
#
# db0's and db1's libraries (69 GB and 101 GB) are deliberately LEFT OFF. They regenerate
# from the CDS tree, which is already in this bucket under db/search, and they overlap db1p's
# 105 GB library almost entirely -- three near-identical copies of the same tagged FASTAs is
# not a backup, it is 275 GB of the same thing. The k2d files are what cost compute, and
# taxonomy/ is what makes a database EXTENSIBLE rather than merely usable: that was the
# lesson of db_v2, whose library/ and taxonomy/ were purged from scratch while only the k2d
# files survived.
#
# NEVER `aws s3 cp --recursive` ON THIS ENDPOINT. Under the Setonix aws-cli/2.13.0 Singularity
# wrapper it recreates the directory tree, copies NO files, and can still exit 0: a db0 CDS
# pull once reported 2,864 successes and put zero .fna on disk. `sync` is not affected and is
# what db_v4's 53 GB CDS backup used (job 51017506, 8 minutes). Top-level files go up as
# individual `cp` calls, directories via `sync`, and EVERYTHING is verified afterwards by
# comparing remote object counts and byte totals against local. A backup nobody checked is
# the same bug as a keep-alive nobody checked.
#
# --endpoint-url is mandatory even though it is also set in the profile: omitting it sends the
# request to Amazon, which answers InvalidAccessKeyId and reads exactly like rotated keys.
set -uo pipefail
module load aws-cli/2.13.0

PROFILE=pawsey1168
BUCKET=s3://pawsey1168-kraken-db
ENDPOINT=https://projects.pawsey.org.au
AWS="aws --profile $PROFILE --endpoint-url $ENDPOINT"
CRYPT=/scratch/pawsey1168/llenzo/crypt
cd "$CRYPT" || exit 1
mkdir -p release/02_reference/scripts/slurm/logs

echo "host     $(hostname)"
echo "started  $(date -Is)"
echo "aws      $(which aws)"
SECONDS=0
rc=0

# The credentials live inside the wrapper's HOME, which tracks the software stack version, so
# a stack upgrade strands them and the symptom is "profile could not be found" for a profile
# you can cat. Fail loudly here rather than halfway through a 280 GB upload.
if ! $AWS s3 ls "$BUCKET/" >/dev/null 2>&1; then
    echo "ABORT: cannot list $BUCKET. Check that \$MYSOFTWARE/setonix/<stack>/.shpc_home/.aws"
    echo "       matches the stack version in $(which aws)."
    exit 1
fi

put_file() {   # put_file <local> <remote-key>
    local src="$1" key="$2"
    [ -f "$src" ] || { echo "  skip (absent) $src"; return 0; }
    local sz; sz=$(stat -c %s "$src")
    printf '  cp   %-22s %8.2f GB -> %s\n' "$(basename "$src")" "$(echo "$sz" | awk '{print $1/1e9}')" "$key"
    $AWS s3 cp "$src" "$BUCKET/$key" --only-show-errors || return 1
}

put_dir() {    # put_dir <local-dir> <remote-prefix>
    local src="$1" key="$2"
    [ -d "$src" ] || { echo "  skip (absent) $src"; return 0; }
    printf '  sync %-22s %8s    -> %s\n' "$(basename "$src")/" "$(du -sh "$src" | cut -f1)" "$key"
    $AWS s3 sync "$src" "$BUCKET/$key" --only-show-errors || return 1
}

backup_db() {  # backup_db <local-db-dir> <name> <with-library: yes|no>
    local P="$1" name="$2" lib="$3"
    echo
    echo "── $name ── $P"
    for f in hash.k2d opts.k2d taxo.k2d seqid2taxid.map excluded_no_cds.tsv; do
        put_file "$P/$f" "db/build/$name/$f" || rc=1
    done
    put_dir "$P/taxonomy" "db/build/$name/taxonomy" || rc=1
    if [ "$lib" = yes ]; then
        put_dir "$P/library" "db/build/$name/library" || rc=1
    else
        echo "  library SKIPPED on purpose (regenerates from db/search CDS)"
    fi
}

SILVER="$CRYPT/research/02_classify/data/silver/db"
backup_db "$CRYPT/release/02_reference/data/db/db1p" db1p yes
backup_db "$SILVER/db0" db0 no
backup_db "$SILVER/db1" db1 no

# ── VERIFY ───────────────────────────────────────────────────────────────────────────────────
# Compare what is on Acacia against what is on disk, per database. Counts AND bytes: the
# cp --recursive bug produced the right tree with no files in it, so a directory existing
# proves nothing.
echo
echo "── verification ──"
verify() {     # verify <local-dir> <name> <with-library>
    local P="$1" name="$2" lib="$3"
    local want_n=0 want_b=0
    for f in hash.k2d opts.k2d taxo.k2d seqid2taxid.map excluded_no_cds.tsv; do
        [ -f "$P/$f" ] && { want_n=$((want_n+1)); want_b=$((want_b+$(stat -c %s "$P/$f"))); }
    done
    for d in taxonomy $([ "$lib" = yes ] && echo library); do
        [ -d "$P/$d" ] || continue
        want_n=$((want_n + $(find "$P/$d" -type f | wc -l)))
        want_b=$((want_b + $(find "$P/$d" -type f -printf '%s\n' | awk '{s+=$1} END {print s+0}')))
    done
    local out got_n got_b
    out=$($AWS s3 ls "$BUCKET/db/build/$name/" --recursive --summarize 2>/dev/null | tail -2)
    got_n=$(echo "$out" | awk '/Total Objects/ {print $3}')
    got_b=$(echo "$out" | awk '/Total Size/ {print $3}')
    printf '  %-6s local %5d files %10.2f GB   acacia %5d objects %10.2f GB   %s\n' \
        "$name" "$want_n" "$(echo "$want_b" | awk '{print $1/1e9}')" \
        "${got_n:-0}" "$(echo "${got_b:-0}" | awk '{print $1/1e9}')" \
        "$([ "${got_n:-0}" -eq "$want_n" ] && [ "${got_b:-0}" -eq "$want_b" ] && echo OK || echo MISMATCH)"
    [ "${got_n:-0}" -eq "$want_n" ] && [ "${got_b:-0}" -eq "$want_b" ]
}
verify "$CRYPT/release/02_reference/data/db/db1p" db1p yes || rc=1
verify "$SILVER/db0" db0 no || rc=1
verify "$SILVER/db1" db1 no || rc=1

echo
echo "── bucket total after ──"
$AWS s3 ls "$BUCKET" --recursive --summarize 2>/dev/null | tail -2
echo
echo "exit     $rc after ${SECONDS}s"
echo "finished $(date -Is)"
exit $rc
