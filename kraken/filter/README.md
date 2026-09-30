# Detection filter: real infections vs classifier artefacts

Deciding which Kraken2 pathogen detections in the SRA screen are real and which are
artefacts. Short version: a real infection leaves a genome-wide signal that grows as we
sequence deeper; a false one does not.

## The problem

Kraken2 assigns reads to species by matching short genome fragments (k-mers). It is
sensitive but it over-calls. Reads from a conserved region shared across many fungi, or
from a close relative of something genuinely present, pile onto a species that was never
there. A raw read count cannot separate these from true detections: soybean rust was
"detected" in almost every wheat sample, which is not credible.

## The insight

Count how much of an organism's genome we actually see (distinct k-mers), not how many
reads landed on it. A real infection accumulates new k-mers as reads increase, so its
cloud rises. A false one stacks reads onto the same small fragment, so its cloud stays
flat no matter how much data arrives.

![examples](figures/example_detections.png)

*Zymoseptoria* and *Puccinia triticina* on wheat rise cleanly (real). *Melampsora*
(poplar rust) on wheat is flat (noise). *Phakopsora* (soybean rust) is the same organism
read on two hosts: a genuine rising infection on soybean, a flat artefact on wheat. The
classifier makes the same call on both; the shape of the cloud tells them apart.

## The signal is shape, not abundance

Two numbers frame it:

- **Level**: how much genome we saw (distinct k-mers). A pile-up sits on a floor of
  ~1,000; a real detection reaches tens of thousands to millions.
- **Slope**: whether coverage grows with depth. Real ~0.7 to 0.9; a pile-up ~0.

A detection needs to fail both (low level *and* flat) to be an artefact. Level alone
rescues saturated real infections (soybean rust is flat only because its genome is
already fully covered); slope alone rescues genuine but shallow infections.

![model](figures/model2d.png)

## It depends on the host

The same organism can be real on its own host and a cross-mapping artefact on another,
so the test is applied per host, not globally.

![slopes](figures/slopes.png)

Read across a row: *Sclerotinia* is real on canola and soybean (blue) but noise on wheat
(red); stripe rust (*Puccinia striiformis*) is real on wheat but noise on maize. The
colours track known host ranges, the check that the method measures real biology.

## formae speciales split the crowded clouds

*Puccinia graminis* and the other rusts on wheat show two overlapping clouds, one rising
and one flat, with a single fitted line running straight through both. The split is
biological: the reads separate by *forma specialis*. Broad, genome-wide coverage lands on
the species node (real); reads biased onto one f.sp. sliver crash to the floor (artefact
of LCA smearing). Treating each f.sp. as its own detection resolves the two clouds.

![graminis](figures/graminis_ratio_split.png)

The current classifier folds this in as an `unbiased` feature (1 minus the f.sp. fraction
of species reads), so a detection dominated by one sub-taxon is penalised.

## How the classifier scores a detection

`detect.py` scores every (host x species) cloud, then labels each point real or crash-out.
Five features, each 0 to 1, combined as a weighted mean; real if score >= 0.5:

| feature   | meaning                                   | weight |
|-----------|-------------------------------------------|:------:|
| tightness | R2 of the upper accumulation line         | 0.20   |
| rising    | slope of that line, clipped               | 0.25   |
| coverage  | distinct / genome ceiling, log-scaled     | 0.20   |
| on_line   | this point's residual vs the fitted line  | 0.15   |
| unbiased  | 1 minus f.sp. fraction                    | 0.20   |

![signal vs crash-out](figures/signal_vs_crashout.png)

## Status

The weights are a provisional operating point, not a calibrated model. Fitting them to
PHI-base labels made coverage dominate and crushed low-abundance real detections; dropping
coverage over-rejected obvious reals (striiformis to 0%). The declared-pathogen ground
truth is biased toward abundant infections, so no automated fit converges cleanly. These
hand weights sit between the failure modes and err toward sensitivity: crash-outs are kept
as a separate labelled group rather than deleted, so a wrong call can still be recovered
downstream. Real calibration waits on the planned read-level alignment, which is the
ground truth this filter should be tuned against. See the memory note
`crypt-kmer-fraction-criterion` for the full arc.

Alternative separation methods tried (a 1D line residual, a GMM cluster, a rarefaction
test) are compared in `figures/method_comparison.png`; none beat the scored per-host model.

## What this changes

The co-infection rate is built from the kept detections only. The filter removes roughly
a fifth of pathogen detections cohort-wide, and the survivors are dominated by the
pathogens each crop actually gets.

## Layout

```
detect.py            classifier: scores every host x species cloud -> data/_classified.tsv
charts.R             one chart per species (63), faceted by host, points coloured by score
charts/              the 63 per-species PNGs
data/                classifier output + figure inputs (large tables gitignored, rebuilt from reports)
figures/             the explanatory figures and their scripts
_scratch/            earlier exploratory versions of the idea, kept for reference
```

Run from the repository root (`phd/01-review`):

```
python3 kraken/filter/detect.py     # -> data/_classified.tsv
Rscript  kraken/filter/charts.R     # -> charts/*.png
```

Figure scripts under `figures/` regenerate the same way (prep `.py` then `.R`). Inputs are
the Kraken2 reports in `kraken/assign/data/reports/` (gitignored data, on Setonix) and the
read-based host calls in `kraken/assign/host/`.
