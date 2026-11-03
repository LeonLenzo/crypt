# 01_refs — who belongs in the database, and what exists for them

Four scripts, in order. Each writes a table the next one filters, and none of them re-query
what an earlier one established.

| | |
|---|---|
| `routes.py` | membership: who is in scope, and under which route |
| `survey.py` | what NCBI actually holds for them, measured per species |
| `select_refs.py` | the reference set: every assembly carrying a gene model |
| `compare_refs.py` | the diff against db_v3, so a database change is reviewable |
| `predict/` | the gene-prediction benchmark, for the species NCBI has not annotated |

**State: DONE.** 4,271 annotated assemblies over 1,547 species.

## The rule, and what it replaced

**Take every assembly that carries a gene model. No floor, no cap, no quality gate.**

That one line replaced three mechanisms, each of which had to be defended separately and one
of which failed open:

- **The depth floor and its no-cap exception.** db_v3 took `min(floor, available)` for most
  species and unlimited depth for the 38 genera STAT had detected, which made STAT a gate on
  database *composition*. When `01_stat` moved out of the pipeline the loader read a missing
  file as "all genera detected" and lifted the cap everywhere. Taking everything is simpler and
  removes the gate. It is affordable because the floor was never the real constraint on
  breadth: the median species has ONE assembly, and pangenomes are open with Heaps gamma ~0.16,
  so extra assemblies on a deep species add sublinear k-mer content.
- **The assembly-level gate.** db_v3 required `scaffold_plus AND has_annotation`, and for six
  resident genera those never co-occur, which excluded *Papiliotrema* (rice's most abundant
  genus at 14.5%) and *Dioszegia* (a wheat core yeast). A reference exists to give a read
  somewhere nearer to park, and that needs real k-mers rather than a good assembly.
- **BUSCO as a threshold.** It is a QC record, kept in `utilities/`, and a per-taxid fallback
  means the screen can never remove a taxon. 50% is a junk filter, not a quality bar.

## Why breadth, stated once

Kraken2's LCA can only decline to resolve a read when two or more database taxa share its
k-mers. db_v2 held 326 taxa, so most reads had exactly one plausible match and got a confident
species call whether or not the species was present. That is why soybean rust was "detected" in
every wheat sample. **Breadth buys competitors, not coverage of more pathogens.**

Assembly count does not inflate a taxon's pull, because Kraken2 stores each distinct minimizer
once. The *Phakopsora* artefact was caused by the **absence** of 37 rust taxa that have
assemblies, not by imbalance among those present.

## Scope is residents, not pathogens

Measured 2026-10-13 against Colvin & Chopra's field maize and sorghum abundances: **32% and 39%
of species-level fungal abundance sat on taxa absent from db_v3.** A read from an organism with
no representation is not left unclassified, it goes to the nearest taxon present, and no
detection threshold can catch that, because the reads and the organism are real and only the
assignment is wrong.

This is why `routes.py` has five routes rather than a pathogen list, and why PHI-base is not the
scope source: `pathogenic` resolves to only 340 fungi.

## Annotation is the constraint, not availability

59% of target species have no gene model anywhere, and NCBI gives no CDS without annotation.
`predict/` measured whether ab-initio prediction can close that: **recall 0.963 / precision
0.945** against NCBI annotation (n=60), **novelty 0.963** against db_v4 (n=59), and
fragmentation does not degrade it (contig-level 0.958 against complete 0.976).

**A production run is not decided, and should be prioritised rather than run over all 2,138.**
The phyllosphere residents score ~0.99 novelty; *Penicillium* and *Alternaria* score under
0.08. Cost at the measured 51.1 SU/genome is ~27% of the allocation for all 2,138.

**Predicted taxa have coverage that is NOT comparable to downloaded taxa.** Over-prediction on
repeat-rich genomes inflates the minimizer denominator (*B. graminis* 3.91x at 0.513 precision),
which deflates every coverage figure computed against it. This is a rule for `05_detect`, not a
caveat about the benchmark.

## Trap

`_expand_taxids()` must pass `intermediate_nodes=True`. The default returns only leaf taxa and
silently drops species nodes that have named strains as children, which hits rust formae
speciales and many plant pathogens directly.
