# 05_detect — is the organism there

One script. **State: BUILT and RUN**, 1,184 runs called. Writes
`data/silver/_classified.tsv` (one row per candidate), `data/gold/detections.tsv` and
`data/gold/coinfection.tsv`.

This step replaced `03_kraken/05_filter`. The method changed, not just the paths: the fitted
`score` is gone and **saturation** took its place.

## The call

All four must hold, tested in this order and always **positively**:

| | | |
|---|---|---|
| `distinct >= 10,000` | absolute floor | not optional |
| `coverage >= 0.01` | abundance floor | is the DNA there |
| `saturation >= 0.05` | shape test | are the reads spread or piled |
| `fsp_frac <= 0.80` | forma specialis | currently inert |

Order matters. The absolute floor is first because a relative floor alone is meaningless on a
small reference: **1% of db_v4's smallest species is 29 k-mers**, reachable by cross-mapping
alone, and 20 species have ceilings under 100,000. The shape test is last because it only
applies where there are enough reads to judge spread.

**Test positively.** The db_v3 table carried NaN in 13.8% of rows, and NaN compares False
against everything, so `if score < FLOOR: continue` *admits* those rows. That inflated a cereal
species count from 32 to 44 on 2026-10-13 before it was caught. `call_one()` is written
`dist is not None and dist >= MIN_DISTINCT` for that reason.

## Where the 2% coverage floor came from

This is the number most likely to be challenged, and CLAUDE.md's old account of it was the
threshold sweep, which reads like a filter tuned to move a result. The honest account:

**It is anchored to the ceiling of a negative control.** The control is obligate rusts on host
genera they carry no USDA record for: real organism DNA that cannot be infecting, so whatever
coverage it reaches is cross-mapping plus deposition. Over **10,314 such rows** it runs

    p75  0.018%    p95  0.046%    p99  0.095%    p99.9  0.221%    max  1.970%

**Raised from 1% to 2% on 2026-10-28.** 1% was defended as "roughly 30x the null's 75th
percentile", which was true but left the floor *inside* the control's range -- five null rows
cleared it. At 2% **zero of 10,314 clear it**. That is a categorically different sentence: not
a percentile ratio but the point at which the negative control is exhausted.

Measured cost of the change:

| floor | calls | runs | undeclared | biotroph calls |
|---|---|---|---|---|
| 1% | 3,150 | 1,370 | 1,873 | 1,617 |
| **2%** | **2,635** | **1,301** | **1,388** | **1,519** |
| 3% | 2,428 | 1,288 | 1,192 | 1,468 |
| 5% | 2,189 | 1,276 | 983 | 1,403 |

Calls fall 16% but runs only 5%, and the obligate biotrophs, which carry the strongest claim,
retain 94%. Past 2% the control is already exhausted, so more stringency buys no extra
justification and spends real evidence.

Note when rebuilding this null: **exclude Triticale.** It is x*Triticosecale, a hybrid whose
genus token matches nothing in USDA, so genuine wheat rusts on it score as off-host and
contaminate the tail (max 53.9% instead of 1.970%).

### What the floor does and does not exclude

The previous wording -- "coverage is not an artefact filter and we do not claim it is" -- is
now too modest in one direction and still correct in the other. Precisely:

- **Cross-mapping: excluded, and validated.** Zero of 10,314 negative-control rows clear 2%.
  *Phakopsora pachyrhizi* on wheat, the worst db_v3 offender, tops out at 0.07% over 1,510
  db_v4 rows and is never called. (CLAUDE.md's "residual Phakopsora is a real finding" entry
  describes db_v3; db_v4 closed it.)
- **Drift and deposition: NOT excluded, and not excludable.** Deposited material is real DNA
  off a real genome and spreads k-mers normally, so no abundance or shape statistic separates
  it from infection. That is answered by the host axis and by obligate biotrophy, not by a
  threshold. Disclose it; do not threshold it away.
- **Saturation is not the artefact filter either.** On the same control its values run p95
  1.57, p99 5.22, max 30.3, against p05 0.069 / median 0.245 for genuinely called host
  pathogens -- the artefacts score *higher*, because at low N the denominator collapses toward
  N*L. It is a shape test that works only where there is depth to judge spread.

Two limits that remain:

- **Coverage is not comparable across database versions.** The denominator is each build's own
  `db_inspect.tsv`. A db_v3 rate and a db_v4 rate are different measurements, and the db_v3
  wheat 47% cannot be restated against db_v4 by arithmetic.
- **The control is built from rusts**, which have large, well-represented genomes in db_v4, so
  it may not generalise to small or thinly-represented taxa. The axenic controls are what
  would generalise it.

## Saturation, and its limits

`saturation = distinct / M(1 - exp(-N·L/M))`: the distinct k-mers observed over what random
sampling at that depth would give. Normalising both axes by `M` makes the null **one universal
curve** rather than one per species, which is what the old fitted `score` could not do.

Three constraints, all measured:

- **It is not scale-free, and must be read against redundancy.** For small x,
  `1 - exp(-x) ~= x`, so coverage approximates redundancy mechanically and everything scores
  near 1. Median saturation falls 0.67 -> 0.23 -> 0.053 -> 0.019 -> 0.002 across redundancy
  decades. That is the mechanism, not biology. Below `MIN_REDUNDANCY` the statistic is silent
  and the call rests on the floors alone.
- **`SAT_FLOOR = 0.05` is provisional, not calibrated.** It sits below the 25th percentile
  (0.091) of detections db_v3's score called real, and above the 460 calls that clear coverage
  with piled reads. Real calls span 0.07 to 0.56 rather than clustering near 1, because the null
  assumes uniform sampling and RNA-seq expression is not uniform. If expression skew varies by
  organism or infection stage then this floor is not constant across taxa. **Calibrating it
  needs the axenic controls.** Until then run `--sensitivity` and report the range.
- **`reads / distinct` is a dead end.** It looked like a 20x discriminator on one run and decayed
  20.8x -> 7.5x -> 3.1x -> 1.0x across four. Recorded so it is not re-walked.

Prior art exists and should be cited rather than claimed: Zebra, KrakenUniq, Metapresence,
YACHT. None of them is this statistic, but the idea of judging a detection by k-mer-space
occupancy rather than read count is not novel.

## Two features that are deliberately not doing what they look like

**The forma-specialis test is INERT** and the script prints `f. sp. feature: INERT` at startup.
`MAX_FSP` would reject a detection dominated by one forma specialis, but db_v4 tags nothing
below species (db_v3 had 222 S1 nodes, db_v4 has zero). It is left in place with a notice rather
than deleted, so a future build with f. sp. nodes cannot silently change behaviour.

**Congener clusters are a flag, never a filter**, and the reason is that collapsing them is
sometimes right and sometimes destroys a real result:

- **Right:** ERR1607767 calls fourteen *Alternaria* species, every one at coverage 5.6 to 7.3%
  and saturation 0.070 to 0.080. Fourteen genuinely co-infecting congeners would differ wildly
  in abundance. One organism's reads spread across fourteen near-identical references looks
  exactly like this. **The uniformity is the diagnostic.**
- **Wrong:** all 96 runs that a collapse rule would flip are *Puccinia*, 82 of them
  *striiformis* plus *triticina*, which differ 46.9% against 6.9% coverage and are separable by
  exact k-mers at a 16-point ANI gap. Collapsing those discards a real co-infection.

Saturation cannot catch congener inflation, and that is not a failure of the statistic: the
reads genuinely do spread across each reference, because the references are nearly the same
sequence. Nor can any per-detection threshold, because each call is individually
well-supported. It is a property of the database's congener density, and db_v4 deliberately has
a lot of that, because breadth was added to give the classifier competitors.
`04_align` is what resolves it.

## Outputs and how to read them

`coinfection.tsv` carries `n_pathogens`, `n_pathogen_genera`, `coinfected`, `coinfected_genus`
and `congener_clusters`. **`coinfected` is species-level and `coinfected_genus` is the
genus-level read; quote both or neither.**

Pathogen status is a **label applied after detection**, from FUNGuild guilds via holobase.
Never a filter applied before: filtering on who is a known pathogen of this host is circular,
because the study exists to find unreported associations.

**No rate out of this step is a prevalence.** The 59.2% over 1,184 runs has four nested
selections in it. See CLAUDE.md, "What is NOT quotable yet".
