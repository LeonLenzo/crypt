# Complexity by host — per-detection real vs crash-out

One image per pathogen species, faceted by host. Each detection is labelled:

- **real** (blue): coverage accumulates with depth — on the fitted accumulation line,
  or has seen ≥10% of the genome (saturated).
- **crash-out** (orange): reads stack on a fixed sliver — fallen below the line, no
  broad coverage.

Classification is per (host × species): fit the real accumulation line to the upper
cloud, require it to rise (else the whole cloud is floor), then split each detection by
residual, with a saturation rescue for the high-depth tip. Hosts from the read-based
call. Cells with too few runs to fit are omitted.

Per-detection labels: `_classified.tsv` (species, host, reads, distinct, cls). Crash-outs
are a separate group, not discarded — kept for alignment (which genes they hit).

Regenerate: `kraken/assign/figures/prep_by_host_individual.py` then
`complexity_by_host_individual.R`.
