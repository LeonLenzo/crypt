# validation — the supplementary prediction figure

Four panels, four claims, all answering "why should a predicted gene model be trusted at
all". They do not answer "did prediction help", which is the MAIN figure and is blocked on
the db2 build.

    a  fragmentation does not degrade recovery   n=60   benchmark by assembly level
    b  predictions are as complete as downloads  n=135  BUSCO vs assembly quality
    c  a congener does not cover a species       n=59   novelty against db_v4
    d  large genomes are over-predicted          n=53   the limitation, stated

Run `python3 prep.py`, then `Rscript validation.R`, then `python3 caption.py`. The caption
is generated from `validation_numbers.tsv` so a quoted number cannot drift from its panel.

## Where the data comes from, and why that is temporary

    a, c, d   research/02_classify/data/gold/predict_benchmark.tsv
    b         research/predict/data/{all135_busco_results.tsv, worklist.tsv}
    b ref     research/superseded/03_kraken/.../busco_scores.tsv   db_v4's own BUSCO record
    b, d      release/02_reference/data/curation/screen_assemblies.tsv   genome sizes

This reads research tables directly, which a `release/` figure may not do: the contract
there is that a figure reads an outputs table and nothing else. When
`supp_data_s4_predictions.tsv` exists, panel b should read it in place of the worklist and
BUSCO pair, and the benchmark needs an outputs table of its own before a, c and d can move.
Until then this lives in `research/`.

## Panel b's quality axis is size deviation, not contig count

Both were available and contig count alone does not separate the failures. Three of the four
sub-90% targets are within normal fragmentation and anomalous in SIZE (0.36x to 1.8x the
genus median: a partial assembly, a MAG and a haplotig-inflated one), and only *P.
brachypodii* fails on fragmentation, at 23,880 contigs against a next-worst 12,392. Plotting
size deviation on x and contig count as marker area lets each failure declare which kind it
is. The genus median needs at least 3 assemblies to mean anything, which drops 8 of 135.

## Panel c's split is db_v4 membership, not reference scope

Novelty was measured against db_v4's sketches, so "is the genus represented" has to mean
represented IN db_v4, taken from the selected rows of its BUSCO record. Splitting on
`supp_data_s3_reference.tsv` instead (reference scope, which is wider) gives 0.851 against
0.980 on 43/16 rather than 0.644 against 0.972 on 22/37. Both say the same thing and the
db_v4 split is the one that matches the measurement.

Neither reproduces the 0.829 / 0.984 pair quoted in the November notes. Those numbers do not
come out of `predict_benchmark.tsv` under any split of it, so they are not used here.

## Panel d is drawn on a linear size axis deliberately

A log axis looked better and would have been dishonest: `geom_smooth` fits on the
transformed scale, so the line would have been a log-x fit while the caption quoted the
linear r (+0.76 against +0.47). Linear axis, linear fit, one r.
