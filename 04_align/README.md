# align: confirm co-infections by competitive read assignment

Kraken2 + the `filter` scores tell us a detection has a real accumulation shape and that the
organism's DNA is present (>=1% genome coverage). They cannot tell us **which genes** the
reads come from, or whether an apparent co-infection is two organisms or one organism whose
reads are shared with a congener in the database. This module answers both by aligning reads
competitively to gene sets, using kallisto.

## Why kallisto

kallisto pseudo-aligns reads to a transcriptome and resolves multi-mapping reads with an EM.
That EM is the point: when a read is compatible with transcripts from several species (a
conserved or shared gene), it is apportioned by the overall likelihood rather than counted
for every species. So:

- an organism that is genuinely present keeps many **uniquely-assigned** reads spread across
  **many genes**;
- an organism that only ever appeared because it shares genome with a truly-present congener
  **collapses** toward zero when both are in the index together.

This is the competitive test the k-mer counts cannot do. (Salmon would do the same with an
added decoy-aware mode, host genome as decoy; kept in reserve if host read bleed is a
problem. minimap2/bwa to genomes give locus-level evidence but slower and with cruder
multi-mapping handling.)

## Design

Reference CDS already exist per species in `03_kraken/01_search/data/cds/pathogen/{accession}/`
(downloaded for the db build; on Setonix scratch). Host CDS via the same host table the db
used. Nothing is downloaded here.

1. **Per-host index.** For each host, build one kallisto index from the CDS of every
   candidate pathogen flagged in that host (score >= 0.7 in `filter/data/_classified.tsv`),
   plus their **congeners** (so within-genus reads compete), plus the **host** CDS (so host
   reads are absorbed, not forced onto a pathogen). Transcripts are tagged with their species
   so abundances aggregate up.
2. **Quant each sample** (the same fastqs used for assign, on scratch) against its host
   index. kallisto quant, paired where available.
3. **Per detection, three read-level facts:**
   - *uniquely-assigned reads*: reads the EM gives to this species alone (from the
     equivalence classes), not shared. Real infection >> 0; cross-map ~ 0.
   - *gene breadth*: number of distinct genes with reads. Real infection covers many genes;
     cross-map hits a few conserved ones.
   - *competitive survival*: est_counts for the species with congeners in the index vs alone.
     A collapse when competitors are present marks cross-mapping.
4. **Confirm or reject** each flagged co-infection from the network on those three facts.

## What it feeds back

This is the ground truth the `filter` thresholds (`score`, `coverage`) were never able to
calibrate against (see `03_kraken/05_filter/README.md` and memory `crypt-kmer-fraction-criterion`).
Once alignment labels a stratified set of detections real/artefact, the score weights and the
coverage floor can be set objectively (ROC / cost-weighted), rather than at the current
provisional operating point.

## First step (before running everything)

Do not quant the whole cohort. Take a **stratified sample** of flagged detections: strong
(high score + high coverage), borderline (just over the 1% floor), and suspected cross-map
(congeneric pairs, e.g. the wheat Fusarium cluster; off-host cases). Quant those, look at the
three facts, and confirm the method separates them before scaling. That set doubles as the
filter calibration set.

## Layout (planned)

```
align.py     build per-host kallisto index from CDS, quant samples, aggregate to species+gene
data/        indices + abundances (gitignored, on Setonix scratch / Acacia)
logs/
```

## Status

Sketch. `align.py` is a skeleton: the plan is fixed, the Setonix specifics (kallisto spack
module, read paths on pawsey1168, SLURM array like `03_kraken/04_assign`) are not wired yet.
