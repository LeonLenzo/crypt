# 02_build — db_v4

Two scripts: `fetch_cds.py` pulls the CDS for everything `01_refs` selected, `build_db.py`
builds the Kraken2 index. **State: DONE.**

## db_v4, as built

```
data/silver/db/db_v4
  hash.k2d  43.0 GB    root minimizers 7,539,317,601   (db_v3: 4,800,752,701)
  k=35  l=31           built in 5h 40m, peak RSS 106 GB
  db_inspect.tsv       the coverage denominator for every 05_detect figure
```

43 GB, not the ~97 GB a linear extrapolation gives: minimizer dedup produced 1.57x db_v3's
minimizer count from 3.5x the sequence. A classify job needs ~43 GB resident.

## Decisions

- **CDS only, never genomic.** The cohort is 100% RNA-seq, so the transcriptome is what reads
  come from. This is not a size argument: a genomic reference is the wrong reference here.
  Wheat CDS is ~0.2 GB against ~15 GB genomic.
- **`db_inspect.tsv` is an output, not a by-product.** Every coverage figure downstream divides
  by a number from it. That is also why **coverage is not comparable across database versions**:
  the denominator is per-build.
- **db_v4 tags nothing below species.** db_v3 had 222 S1 nodes; db_v4 has zero. `05_detect`'s
  forma-specialis feature is therefore inert, and says so at startup rather than being deleted.

## Validated against the off-host null

The cleanest available negative control is obligate rusts on hosts they cannot infect. db_v3
reached **13% coverage at its 75th percentile and up to 74%**, which meant the 1% abundance
floor could not be anchored to anything. db_v4 gives **0.032% at the 75th percentile, 1.97%
maximum, and 0 of 6,285 rows called.**

That collapse is the single strongest result this module produced, and it is quotable, because
it measures the method rather than the world.

## Cost lessons, both learned the expensive way

- **Do not quote a per-genome cost before an array drains.** Early means over n=3 and n=5 gave
  13.4 then 16.7 SU; the true figure over 120 was **51.1**. Short jobs finish first, so an
  early mean over a size-heterogeneous sample is biased low.
- **The MEMORY request sets the CPU bill**, not `--cpus-per-task`, at roughly 1.9 GB per
  physical core. `AllocCPUS` double-counts because of SMT.
