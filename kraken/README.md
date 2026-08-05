# Kraken2 orthogonal classification

## Rationale

The primary co-infection detection pipeline uses NCBI STAT pre-computed k-mer taxonomy
(see `02_filter_runs.py`). A kallisto pilot on 45 rust-dominated MAL runs revealed that
STAT fails to distinguish *Puccinia striiformis* f. sp. *tritici* (PST) from other wheat
rusts — the diagnostic k-mer sets are too similar at the SRA database scale. This raised
the question of whether other co-infection detections share the same reliability problem.

Kraken2 with a purpose-built, curated database provides an orthogonal check:

- STAT is pre-computed across the entire SRA; Kraken2 is run on-demand against a
  controlled reference set
- The Kraken2 DB contains only PHI-base eukaryotic pathogen CDS + PHI-base plant host
  CDS — no environmental microbes, no human, no ambiguous sequences
- CDS sequences (not whole genome) are used to match the spliced nature of RNA-seq reads
  and keep the database small enough to build and query locally

The primary paper figure compares STAT detections against Kraken2 detections for the
same runs, quantifying where the two methods agree and where they diverge.

## Database design

| Component | Source | Species | Sequences |
|---|---|---|---|
| Pathogens | PHI-base fungi + oomycetes | 101 annotated seeds | CDS from genomic |
| Hosts | PHI-base plant hosts | 92 annotated seeds | CDS from genomic |

Unannotated assemblies (no gene model) are excluded — without CDS coordinates there
is no `cds_from_genomic.fna` file to use.

Seeds with annotation but absent from the SRA data (`runs.tsv`) are still included for
completeness, consistent with the approach taken for pathogens.

Known gaps: *Nicotiana benthamiana* has no annotated assembly in NCBI Datasets as of
2026-08 despite having published chromosome-level assemblies (Bally et al. 2022); the
annotated NbBMZ assembly is hosted externally on Sol Genomics Network.

All FASTA headers are retagged to `>kraken:taxid|TAXID|original_header` before adding
to the library, which bypasses the 15 GB accession2taxid lookup and corrects stale
PHI-base taxids using the `organism.taxId` field from `assembly_data_report.jsonl`.

## Scripts

```
screen_refs.py   Query NCBI Datasets for the best assembly per seed taxid.
                 Ranks by: annotation present > gene count tier > assembly level >
                 RefSeq (tie-breaker) > release date.
                 --include-hosts also screens the 180 PHI-base plant host seeds.
                 Output: ref_screen.tsv

build.py         Download CDS FASTAs and build the Kraken2 database.
                 Reads ref_screen.tsv for accessions. Skips seeds without annotation
                 (fasta_type=genome) and deduplicates seeds sharing the same assembly.
                 --download-only   fetch FASTAs without building
                 --build-only      build from already-downloaded FASTAs
                 Output: kraken/db/  (gitignored)
                         kraken/cds_from_genomic/  (gitignored)

classify.py      Stream 500k reads per run from ENA FTP and classify with Kraken2.
                 Reads mal_runs.json / hal_runs.json from output/01_fetch_runs/data/.
                 Resumable via append-only kraken_cache.jsonl + index file.
                 --reports-dir   save raw kraken2 report.txt per run (reproducibility)
                 --limit N       process at most N runs (testing)
                 Output: output/kraken_classify/data/kraken_cache.jsonl
```

## Parameters

| Parameter | Value | Rationale |
|---|---|---|
| Reads per run | 500,000 | Sufficient for detection; limits ENA bandwidth |
| Confidence | 0.1 | Standard; reduces spurious species-level calls |
| Workers | 8 | Respects ENA concurrent connection limit |
| Kraken2 threads/worker | 4 | workers × threads ≤ available cores |

## Running

```bash
# 1. Screen references (already done; re-run if phibase_db.json is rebuilt)
python kraken/screen_refs.py --include-hosts

# 2. Download CDS FASTAs
python kraken/build.py --download-only

# 3. Build database (delete kraken/db/ first if rebuilding)
rm -rf kraken/db/
conda run -n kraken2 python kraken/build.py --build-only

# 4. Classify (test run first)
conda run -n kraken2 python kraken/classify.py \
    --mode mal --db kraken/db --limit 100 \
    --reports-dir output/kraken_classify/data/reports

# Full MAL run
conda run -n kraken2 python kraken/classify.py \
    --mode mal --db kraken/db \
    --reports-dir output/kraken_classify/data/reports
```

## Output format

`kraken_cache.jsonl` — one JSON object per run:

```json
{
  "run": "SRR2079621",
  "layout": "PAIRED",
  "error": null,
  "ts": "2026-08-05T12:11:13",
  "pct_classified": 37.0,
  "pct_unclassified": 63.0,
  "n_reads": 500000,
  "species": [
    {"taxid": 318829, "name": "Pyricularia oryzae", "pct": 35.78, "reads": 178924},
    {"taxid": 5507,   "name": "Fusarium oxysporum", "pct": 0.14,  "reads": 699}
  ]
}
```

Species are all taxa at rank `S` in the Kraken2 report with `pct > 0`, sorted by
percentage descending. Cross-reference `taxid` against `phibase_db.json` to distinguish
host species from pathogen species.

Raw report files (saved with `--reports-dir`) preserve the full Kraken2 taxonomy tree
including genus, family, and root counts — sufficient to recompute any summary statistic
without re-running classification.
