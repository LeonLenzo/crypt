# kallisto: confirm co-infections by competitive read assignment

Kraken2 + the `05_filter` scores tell us a detection has a real accumulation shape and that the
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
(downloaded for the db build; on Setonix scratch). Host CDS sit under `cds/host/`, resolved
through `01_search/data/host_candidates.tsv`: Ensembl-sourced hosts in species-named
directories, NCBI-sourced hosts in accession-named ones. Nothing is downloaded here.

1. **Per-host index.** One kallisto index per host, holding the tagged CDS of:
   - the **host**, so host reads are absorbed rather than forced onto a pathogen;
   - **every species Kraken2 saw in that host at any level**, not only the flagged ones. A
     species that is genuinely present but missing from the index has its reads pushed onto
     whatever relative *is* present, inflating exactly what this module is meant to test. One
     assembly each makes this cheap insurance;
   - the **measured k-mer neighbours** of the flagged species (>=1% containment,
     `02_build/data/neighbours.tsv`).

   Every transcript header is rewritten to `{species_slug}|{accession}|{transcript_id}` so
   kallisto's per-target abundances aggregate to species while keeping a gene-level id for
   breadth.

2. **Quant each sample** (the same fastqs used for assign, on scratch) against its host
   index. kallisto quant, paired where available.
3. **Per detection, three read-level facts:**
   - *uniquely-assigned reads*: reads the EM gives to this species alone (from the
     equivalence classes), not shared. Real infection >> 0; cross-map ~ 0.
   - *gene breadth*: number of distinct genes with reads. Real infection covers many genes;
     cross-map hits a few conserved ones.
   - *competitive survival*: est_counts for the species with neighbours in the index vs
     alone. A collapse when competitors are present marks cross-mapping.
4. **Confirm or reject** each flagged co-infection from the network on those three facts.

## Why neighbours are measured, not inferred from taxonomy

The original plan expanded each candidate to its **congeners**. That is the wrong unit.
Phylogenetic closeness is an amino-acid-level statement, while Kraken2 and kallisto both match
exact 31-mers, and at congeneric divergence synonymous substitutions break a 31-mer every few
bases. So containment, not rank, predicts cross-assignment.

`02_build/data/neighbours.tsv` holds it, computed from the saturation module's cached FracMinHash
sketches (`03_kraken/utilities/saturation/data/sketches/`, k=31, scaled=1000): 398 pairs at
>=1% containment across the 45 flagged species. What it showed:

- *Ascochyta rabiei* in *Parastagonospora nodorum* is **0.13%**, and its sharing with
  *Alternaria alternata* (0.16%) is higher, so the 470 low-level *Ascochyta* detections on
  wheat are generic Dothideomycete background, not *P. nodorum* bleed. Real congener
  cross-talk for scale: *F. graminearum* / *F. pseudograminearum* **27.5%**.
- Cross-talk above 1% is almost entirely within-genus (81 same-genus pairs, 1 cross-family).
  Family-level expansion buys nothing.
- **11 of 45 flagged species have no neighbour above 1%**, including all three *Puccinia*,
  *Zymoseptoria tritici* and *P. nodorum*. The competitive-survival arm cannot run for them;
  they rest on unique reads and gene breadth alone.
- *P. striiformis* / *P. triticina* is **0.08%**, below *Ascochyta* / *P. nodorum*, yet they
  co-occur in 191 wheat runs. The chapter's largest co-infection signal is therefore not
  k-mer cross-talk.

## One representative assembly per species

Chosen on BUSCO completeness (`03_kraken/utilities/saturation/data/sketch_index.tsv`,
`complete_pct`), tie-broken on scaffold N50. Intraspecific pangenome redundancy splits reads
across near-identical transcripts and deflates the unique-read count that is the main signal.

The saturation module's `pangenome_ratio` says when one assembly is not enough. Most flagged
species sit at 1.1 to 1.8, so one assembly holds 60 to 90% of species k-mer space. Two do not:
*Rhizoctonia solani* at 14.3 (one assembly covers **7%**) and *Fusarium oxysporum* at 5.25
(**19%**), both known species complexes. These get `MULTI_N` assemblies instead of one. That is
a compromise, not a fix: their pangenomes are open and `n_for_95pct` is 48 and 56.

## What it feeds back

This is the ground truth the `filter` thresholds (`score`, `coverage`) were never able to
calibrate against (see `03_kraken/05_filter/README.md` and memory `crypt-kmer-fraction-criterion`).
Once alignment labels a stratified set of detections real/artefact, the score weights and the
coverage floor can be set objectively (ROC / cost-weighted), rather than at the current
provisional operating point.

## The stratified first pass

Do not quant the whole cohort. `01_select/kallisto_select.py` writes
`01_select/data/strata.tsv`: **200 detections over 174 runs and all 19 flagged hosts**, six
strata of 40.

| stratum | runs | species | hosts | expectation |
|---|---|---|---|---|
| `strong` (cov >= 0.10, score >= 0.90, uncontested) | 40 | 20 | 10 | confirm |
| `borderline` (cov < 0.02, just over the floor) | 38 | 23 | 14 | calibrates the floor |
| `congeneric` (>=2 of a genus in one run) | 27 | 25 | 6 | collapse if cross-mapping |
| `off_host` (species' runs sit on other hosts) | 34 | 13 | 9 | reject |
| `discordant` (cov >= 0.10 but score < 0.90) | 40 | 18 | 14 | criteria conflict |
| `mid` (0.02 <= cov < 0.10) | 40 | 20 | 17 | mid-range |

The last two exist because the first four label only **45%** of the 3,980 flagged detections.
`discordant` is the largest gap at 1,630 detections, dominated by *P. striiformis* (698): the
coverage says the genome is largely present while the score says the accumulation shape is
marginal. That is exactly where the two filter criteria disagree, so it is the most informative
class for calibrating them, and omitting it would have meant calibrating on the extremes only.

Selection is deterministic (no RNG) and reproducible from the committed inputs. Each stratum
takes one detection per qualifying host first, then fills round-robin over host x species,
spread along coverage. Without that first pass a flat quota silently drops whole hosts: at 40
per stratum it lost *Zea mays*, the second largest with 421 flagged detections.

The sharpest single test in the set: *Puccinia striiformis* on sweet potato at **0.476
coverage, score 0.737**, a confident-looking wheat rust on a host it cannot infect.

## Layout

Numbered = a pipeline step, in order. Unnumbered = support.

```
04_kallisto/
  _common.py                       shared paths and helpers
  01_select/  kallisto_select.py       the stratified set      -> data/strata.tsv      (tracked)
  02_build/   kallisto_neighbours.py   measured containment    -> data/neighbours.tsv  (tracked)
              kallisto_build.py        tagged CDS + index      -> data/idx/            (gitignored)
  03_quant/   kallisto_quant.py        kallisto quant per run  -> data/{host}/{run}/   (gitignored)
  04_confirm/ kallisto_confirm.py      three facts, calls      -> data/calls.tsv
  slurm/      kallisto_neighbours.slurm  kallisto_build_index.slurm
```

`select` precedes `build` because the selection defines the scope: only hosts in the stratified
set need an index. The selection guarantees all 19 flagged hosts are represented, so nothing is
dropped by that dependency.

Index sizes, from `kallisto_build.py --all --dry-run` on Setonix:

| host | species | assemblies | CDS |
|---|---|---|---|
| *Triticum aestivum* | 125 | 138 | 3.38 GB |
| *Mangifera indica* | 64 | 77 | 2.08 GB |
| *Ipomoea batatas* | 63 | 76 | 2.05 GB |
| *Hordeum vulgare* | 60 | 73 | 1.76 GB |
| *Eleusine coracana* | 9 | 8 | 0.26 GB |

## Reads and the scratch purge

The reads live only on Setonix scratch, under a 21-day purge keyed on **atime**. They are not
backed up: 11.9 TB for the flagged runs against ~930 GB free on Acacia, and they are public SRA
accessions, so `03_kraken/03_select/data/run_list.tsv` plus the committed `manifest.tsv` is the
record.

`03_kraken/03_select/freshen.py` keeps what is in active use alive. `/scratch` is Lustre without
`noatime`, so reading **4 KB** of a file resets its clock as well as reading all of it: a full
pass over 5,990 files is 24 MB of I/O. All 3,225 runs were reset on 2026-10-02, so the deadline
is **2026-10-23**. `crontab` is blocked for this account, so run `freshen.py --report` at the
start of a session rather than scheduling it.

## What the analysis is looking for

Per detection, three read-level facts. Two work now; `survival` needs a solo index per species,
which `02_build` does not yet produce, so it is blank in the output rather than approximated.

| fact | real infection | cross-mapping |
|---|---|---|
| `unique_reads` (EC maps to one species) | many | ~0 |
| `gene_breadth` (distinct genes with reads) | broad | a few conserved genes |
| `survival` (with competitors vs alone) | holds | collapses |

The six strata are predictions, and the analysis is whether they separate:

| stratum | expected | role |
|---|---|---|
| `strong` | high unique, broad genes | positive control |
| `off_host` | **real sequence, confirms** | provenance, not cross-mapping |
| `congeneric` | collapses under competition | the cross-mapping test |
| `discordant` | unknown | **the open question** |
| `mid`, `borderline` | unknown | the calibration band |

Flat coverage cuts inherit a known, unfixed bias: the k-mer denominator varies ~10x with
database representation (*Puccinia* ~8.9M minimizers per species against *Alternaria* 913k,
because 24 close *Alternaria* species collapse shared k-mers to the genus node). *Alternaria*
clears 1% at a tenth of the reads *Puccinia* needs and is 2,985 of wheat's 5,506
above-criterion detections, largely as a saprophyte. So `discordant` and `mid` are enriched for
densely-represented genera. See memory `crypt-kmer-fraction-criterion`; the fix is Leon's call
and not yet made.

**`off_host` is not a negative control.** Off-host signal in this cohort is usually real
organism DNA, not a k-mer artefact: the obligate-rust-on-impossible-host null reaches 13%
coverage at its 75th percentile and up to 74%, and the Phakopsora-on-wheat case has 94% of
248,504 reads mapping to *Phakopsora* against 18 to *P. striiformis*. Expect these detections to
confirm as genuine sequence. What they test is provenance, active infection versus deposited
spores or a mis-called host, and that is answered by **which** genes are expressed rather than
how many reads are unique.

**The real negative control is already in every index, for free.** Each host index holds 9 to
153 species but only a handful are flagged in any given run. The unflagged species are a null
distribution: they should show near-zero unique reads in that run, and they bound what the
method calls when nothing is there. Read that distribution before interpreting any positive.

**Then read `strong` against `congeneric`.** `strong` should hold up and `congeneric` should
show collapse where cross-mapping is real; that pair, not strong-versus-off_host, is the
method's internal check.

Because every library is RNA-seq, gene identity carries more than gene count: an actively
growing pathogen expresses in-planta-induced genes, while deposited spores give a
spore/dormancy profile. That distinction is the part of "DNA presence is not proof of active
infection" this module can actually address.

**The deliverable is the thresholds.** With a labelled set, `03_kraken/05_filter`'s `score`
weights and the 1% coverage floor can be fit by ROC or a cost-weighted criterion, replacing the
provisional 0.7 and 0.01 that nothing has justified. The cost asymmetry has to be stated: a
missed co-infection and a false one are not equally bad for the chapter's claim. `discordant`
carries the most weight here, being 1,630 detections where coverage says the genome is largely
present and score says the shape is marginal, dominated by *P. striiformis*.

## The limit: cross-mapping yes, database absence no

For 11 flagged species, including all three *Puccinia*, *Zymoseptoria tritici* and
*P. nodorum*, nothing in the index shares >=1% k-mer space. For those, `unique_frac` near 1
means "matches this reference and nothing else present in the index", which is also exactly
what an absent relative produces. The first validation result is this case, not a finding:
*Z. tritici* on *Hordeum vulgare*, 0.9958 unique fraction over 10,153 genes, with zero
*Zymoseptoria* congeners in db_v3 and *Z. passerinii* being the barley pathogen.

So `calls.tsv` resolves cross-mapping and not database absence. The `n_neighbours` column flags
every detection where the distinction applies. Closing it needs one of:

- **predicted CDS for the missing congeners.** 23 congener species across the nine affected
  genera have assemblies but no annotation, including *Puccinia hordei* (barley brown rust,
  chromosome-level, N50 10.9 Mb) and all five *Zymoseptoria* congeners. Homology-guided
  prediction only, used as competitors and never for positive calls about the predicted species,
  because predicted gene sets are biased toward the guide's conserved genes. See
  `03_kraken/utilities/annotation_gap/`.
- **allelic structure** (`05_hisat2`, not built). Align the same CDS references with hisat2,
  pileup, and ask whether the reads assigned to one species are one population or two. Needs no
  new reference and no predicted sequence. Discriminators, strongest first: divergence magnitude
  calibrated against our own conspecific versus congeneric assemblies; haplotype coherence
  within reads; uniformity across genes. Variant allele frequency is weak here because RNA-seq
  carries allele-specific expression. *Puccinia* are dikaryotic and their two nuclear haplotypes
  differ by ~1-2%, which mimics a two-species mixture on coherence alone, so magnitude is the
  discriminator that matters. Run over all 174 runs, not a subset: the `strong` stratum is the
  null distribution, and depth filtering (~47% of detections reach 5x) belongs at analysis.
  Cost is ~15 node-hours with hisat2; stream to pileup and discard BAMs. Requires `samtools`
  and `bcftools` spack installs, neither currently present.

## Reference type: CDS only

The cohort is 100% RNA-Seq (ENA, 60-run sample: 95% transcriptomic, 5% metatranscriptomic), so
CDS is the correct reference throughout and genomic references are wrong rather than merely
larger. kallisto pseudoaligns to a transcriptome by construction; reads spanning exon-exon
junctions do not exist contiguously in genomic sequence. Scale makes the same point: wheat CDS
is 204 MB against a ~15 Gb genome. Do not mix reference types across the analysis.

## Status and next steps

Implemented: `01_select`, `02_build`, `03_quant`, `04_confirm`. The index build is validated end
to end and `kallisto inspect` confirms each index. `survival()` and the threshold fit in
`04_confirm` raise `NotImplementedError`; survival needs a solo index per species, and the
calibration waits on the full results.

**Quant is complete: 174/174 runs, abundance + bus, no failures.** Nothing is in flight.

`04_confirm` has never produced a full `calls.tsv`. Two problems were found and fixed: it
printed nothing for 12 minutes (computing silently, printing once at the end), so a stuck run
looked identical to a working one, and it re-parsed each host's million-line `transcripts.txt`
once per run. There is now per-run progress output and a per-host transcripts cache. Note that
`matrix.ec` is **not** cacheable across runs: kallisto bus writes the equivalence classes
actually observed in that run, so an ec id means different things in different runs, while
`transcripts.txt` is identical per host. Caching the ec-to-species map would corrupt the
unique/shared split.

The one data point so far, from a single-run smoke test: *Z. tritici* on *Hordeum*, unique_frac
0.9958, 10,153 genes, 0 neighbours. That is a textbook database-absence case (*Z. passerinii* is
the real barley pathogen and is absent from db_v3), not a confirmed co-infection.

1. Run `04_confirm --all` over the 174 runs to produce `calls.tsv`, then read the per-stratum
   table: does `strong` hold, and does `congeneric` collapse where cross-mapping is real.
   `off_host` is **not** a negative control; off-host signal is real DNA.
2. Speed up `04_confirm` (cache the per-index ec-to-species map) before scaling past 174.
3. Migrate `03_kraken/05_filter/detect.py` to the shared `detections()` predicate, verifying
   `_classified.tsv` byte-for-byte.
4. Back up `03_quant/data/` to Acacia. The indices are regenerable, the quant output is not, and
   the reads expire 2026-10-23.
5. Decide on `05_hisat2`.
6. Deferred: predicted congener CDS; the solo-index build `survival` needs; the threshold fit.

**Never run `kallisto_confirm.py` over a bare ssh call.** Doing so computes on a Setonix login
node, and when the ssh client times out the remote python keeps running as an orphan, racing the
real SLURM job for the same `calls.tsv`. Always `sbatch`.
