#!/bin/bash
#SBATCH --job-name=busco_screen
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=24:00:00
#SBATCH --output=kraken/slurm/logs/busco_screen_%j.out
#SBATCH --error=kraken/slurm/logs/busco_screen_%j.err

# Download CDS FASTA + run BUSCO transcriptome mode for all candidates.
# Workers × cpus-per-busco = total CPUs (16 × 8 = 128).
#
# Prerequisites:
#   1. python kraken/screen_refs.py  →  kraken/ref_candidates.tsv
#   2. git pull on Setonix to get latest ref_candidates.tsv
#
# Lineage DBs are downloaded automatically inside this script to
# $SCRATCH/busco_downloads/lineages/ (BUSCO's natural default with HOME=$SCRATCH).
#
# Run from crypt/: sbatch kraken/slurm/busco_screen.slurm

set -euo pipefail

module load python/3.11.6
module load singularity/3.11.4-nompi

SCRATCH="/scratch/pawsey1168/llenzo"
GENOMES_DIR="$SCRATCH/kraken/cds_v2"
BUSCO_OUT="$SCRATCH/kraken/busco"
BUSCO_DB_PATH="$SCRATCH/busco_downloads"
BUSCO_SIF="$SCRATCH/containers/busco_v5.7.0.sif"
DATASETS_BIN="$SCRATCH/bin/datasets"

# Bootstrap NCBI datasets CLI if not already present
if [ ! -x "$DATASETS_BIN" ]; then
    echo "Downloading NCBI datasets CLI …"
    mkdir -p "$SCRATCH/bin"
    curl -fsSL -o "$DATASETS_BIN" \
        https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
    chmod +x "$DATASETS_BIN"
fi

# Busco wrapper: write to .tmp then mv to avoid "text file busy" if a prior
# job's process still holds the file open.
cat > "$SCRATCH/bin/busco.tmp" << EOF
#!/bin/bash
exec singularity exec --bind $SCRATCH:$SCRATCH --env HOME=$SCRATCH $BUSCO_SIF busco "\$@"
EOF
chmod +x "$SCRATCH/bin/busco.tmp"
mv "$SCRATCH/bin/busco.tmp" "$SCRATCH/bin/busco"

export PATH="$SCRATCH/bin:$PATH"
echo "datasets version: $(datasets version)"
echo "busco version:    $(busco --version)"

# Download lineage DBs — check prfl/ count as health signal (empty prfl = inode-limit
# corruption; re-download in that case even if the directory exists).
# Downloads via curl directly to scratch (no HOME involved; busco --download sends
# to HOME which has a 1 GB quota and silently fails for large databases).
mkdir -p "$BUSCO_DB_PATH/lineages"
for lineage in ascomycota_odb10 basidiomycota_odb10 stramenopiles_odb10 fungi_odb10; do
    prfl_dir="$BUSCO_DB_PATH/lineages/$lineage/prfl"
    prfl_count=$(ls "$prfl_dir" 2>/dev/null | wc -l)
    if [ ! -d "$BUSCO_DB_PATH/lineages/$lineage" ] || [ "$prfl_count" -eq 0 ]; then
        echo "Downloading BUSCO lineage: $lineage (prfl_count=$prfl_count) …"
        rm -rf "$BUSCO_DB_PATH/lineages/$lineage"
        TGZ="$BUSCO_DB_PATH/lineages/${lineage}.tar.gz"
        curl --retry 3 --retry-delay 10 -fsSL \
            -o "$TGZ" \
            "https://busco-data.ezlab.org/v5/data/lineages/${lineage}.2024-01-08.tar.gz"
        tar -xzf "$TGZ" -C "$BUSCO_DB_PATH/lineages"
        rm -f "$TGZ"
        echo "  → prfl after download: $(ls $prfl_dir 2>/dev/null | wc -l) files"
    else
        echo "Lineage OK: $lineage (prfl=$prfl_count)"
    fi
done

mkdir -p kraken/slurm/logs
mkdir -p "$GENOMES_DIR"
mkdir -p "$BUSCO_OUT"

echo "Starting BUSCO screen: $(date)"
echo "Genomes dir:  $GENOMES_DIR"
echo "BUSCO out:    $BUSCO_OUT"
echo "BUSCO DBs:    $BUSCO_DB_PATH"

python kraken/busco_screen.py \
    --genomes-dir    "$GENOMES_DIR" \
    --busco-out      "$BUSCO_OUT" \
    --busco-db-path  "$BUSCO_DB_PATH" \
    --workers        16 \
    --cpus-per-busco 8

echo "BUSCO screen complete: $(date)"
echo ""
echo "Next steps:"
echo "  1. Copy kraken/busco_scores.tsv back to local (if run on Setonix)"
echo "  2. python kraken/filter_refs.py"
echo "  3. sbatch kraken/slurm/kraken_build.slurm"
