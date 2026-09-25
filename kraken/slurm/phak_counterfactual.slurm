#!/bin/bash
#SBATCH --job-name=phak_cf
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=230G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/phak_cf_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/kraken/slurm/logs/phak_cf_%j.err
set -euo pipefail
# bwa and samtools are NOT spack packages here: the 2026-09-25 run died at step 4 with
# "Error: Spec 'bwa' matches no installed packages" after 6 minutes of work. Plain
# modules, matching phak_align.slurm which is why that one worked.
module load kraken2/2.1.2-qk2lek6
module load bwa/0.7.17 samtools/1.15
BASE=/scratch/pawsey1168/llenzo/crypt/kraken
CDS=$BASE/search/data/cds/pathogen
READS=$BASE/select/data/reads
DBV3=$BASE/build/data/db_v3
W=$BASE/assign/data/phak_cf; mkdir -p $W; cd $W
RUN=ERR3170248

echo "=== 1. per-read kraken2 (db_v3) $(date) ==="
kraken2 --db $DBV3 --threads 32 --confidence 0.15 --minimum-hit-groups 3 \
  --paired $READS/${RUN}_1.fastq.gz $READS/${RUN}_2.fastq.gz \
  --output $W/kout.txt --report /dev/null >/dev/null 2>&1

echo "=== 2. extract Phakopsora(170000) read-pair IDs ==="
awk -F"\t" "\$3==170000 {print \$2}" $W/kout.txt > $W/phak_ids.txt
wc -l $W/phak_ids.txt

echo "=== 3. pull those read pairs ==="
python3 - << "PY"
ids=set(l.strip() for l in open("phak_ids.txt"))
import gzip
for m in (1,2):
    fin=gzip.open(f"/scratch/pawsey1168/llenzo/crypt/kraken/select/data/reads/ERR3170248_{m}.fastq.gz","rt")
    fout=open(f"phak_{m}.fq","w")
    while True:
        h=fin.readline()
        if not h: break
        s=fin.readline(); p=fin.readline(); q=fin.readline()
        rid=h[1:].split()[0]
        if rid in ids: fout.write(h+s+p+q)
    fout.close()
PY
echo "  extracted: $(($(wc -l < phak_1.fq)/4)) pairs"

echo "=== 4. combined tagged reference (Phakopsora + P.striiformis CDS) ==="
: > combined.fa
for a in GCA_025201825.1 GCA_943846425.1; do
  sed "s/^>/>PHAK|/" $CDS/$a/ncbi_dataset/data/$a/cds_from_genomic.fna >> combined.fa; done
for a in GCF_021901695.1 GCA_002920065.1 GCA_021901655.1; do
  sed "s/^>/>PSTR|/" $CDS/$a/ncbi_dataset/data/$a/cds_from_genomic.fna >> combined.fa; done
bwa index combined.fa 2>/dev/null

echo "=== 5. map + tally $(date) ==="
bwa mem -t 32 combined.fa phak_1.fq phak_2.fq 2>/dev/null | samtools view -b -o aln.bam -
echo "--- flagstat ---"; samtools flagstat aln.bam
echo "--- primary reads by reference species x MAPQ ---"
samtools view -F 0x904 aln.bam | awk "{split(\$3,a,\"|\"); print a[1], (\$5>=20?\"MAPQ>=20(unique)\":\"MAPQ<20(ambiguous)\")}" | sort | uniq -c
echo "=== DONE $(date) ==="
