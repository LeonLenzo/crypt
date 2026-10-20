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

## Where the 1% coverage floor came from

This is the number most likely to be challenged, and CLAUDE.md's old account of it was the
threshold sweep, which reads like a filter tuned to move a result. The honest account:

**1% is an abundance floor, not an artefact filter.** It asks whether enough of the organism's
k-mer space is present to say its DNA is in the library. It was never meant to establish that
the organism is actively infecting.

Under db_v3 the floor could not be anchored. The cleanest negative control available is
obligate rusts on hosts they cannot infect, and that null reached **13% coverage at its 75th
percentile and up to 74%**, far above the floor, because the signal is real organism DNA from
spore drift or a mis-called host. So the sweep was all the justification there was, and the
sweep is not a justification: wheat went 0% -> 80%, 1% -> 47%, 2% -> 40%, 5% -> 30%,
10% -> 23%, and any of those could be argued for.

**db_v4 changed that.** The same off-host null now sits at **0.032% at the 75th percentile,
1.97% maximum, with 0 of 6,285 rows called.** 1% is therefore roughly 30x the null's 75th
percentile, which is a basis rather than a preference.

Two things this still does not give:

- **Coverage is not comparable across database versions.** The denominator is each build's own
  `db_inspect.tsv`. A db_v3 rate and a db_v4 rate are different measurements, and the db_v3
  wheat 47% cannot be restated against db_v4 by arithmetic.
- **Coverage is not proof of active co-infection.** Disclose that; do not threshold it away.

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
