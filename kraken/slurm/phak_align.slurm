#!/bin/bash
#SBATCH --job-name=phak_aln
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/phak_aln_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/phak_aln_%j.err
set -euo pipefail
module load bwa/0.7.17 samtools/1.15
CDS=/scratch/pawsey1168/llenzo/crypt/kraken/search/data/cds/pathogen
W=/scratch/pawsey1168/llenzo/crypt/kraken/assign/data/phak_cf; cd $W

echo "=== combined tagged reference $(date) ==="
: > combined.fa
for a in GCA_025201825.1 GCA_943846425.1; do
  sed "s/^>/>PHAK|/" $CDS/$a/ncbi_dataset/data/$a/cds_from_genomic.fna >> combined.fa; done
for a in GCF_021901695.1 GCA_002920065.1 GCA_021901655.1; do
  sed "s/^>/>PSTR|/" $CDS/$a/ncbi_dataset/data/$a/cds_from_genomic.fna >> combined.fa; done
grep -c "^>PHAK" combined.fa; grep -c "^>PSTR" combined.fa
bwa index combined.fa 2>/dev/null

echo "=== map 248,504 Phakopsora-called read pairs $(date) ==="
bwa mem -t 16 combined.fa phak_1.fq phak_2.fq 2>/dev/null | samtools view -b -o aln.bam -
echo "--- flagstat ---"; samtools flagstat aln.bam
echo "--- primary reads: reference species x confidence ---"
samtools view -F 0x904 aln.bam | awk "{split(\$3,a,\"|\"); print a[1], (\$5>=20?\"unique(MAPQ>=20)\":\"ambiguous(MAPQ<20)\")}" | sort | uniq -c
echo "=== DONE $(date) ==="
