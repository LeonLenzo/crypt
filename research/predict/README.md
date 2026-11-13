# predict

Gene prediction for the species `release/02_reference` reaches but cannot build from.

Unnumbered on purpose. CLAUDE.md's rule is that a numbered directory is a pipeline step in a
fixed order and an unnumbered one is support; this is a method evaluation feeding
`release/02_reference`, not a step between `02_classify` and `03_align`. `02a_predict` was
considered and rejected, on the precedent of `01a_Literature` becoming `01_literature`.

Moved here from `02_classify/01_refs/predict/` on 2026-11-09, with `sample.py`,
`fetch_genomes.py` and `compare.py`. It sat under `01_refs` because gene prediction is
reference construction; it moves out because reference construction now lives in
`release/02_reference` and this is the research half that feeds it.

## The problem

`release/02_reference` reaches 5,855 species with an assembly. Only 2,319 carry a gene model,
so 3,536 have nothing to extract CDS from. They are not equal:

    genus already has annotated assemblies   2,604   predicting here mostly adds competitors
    genus has NO annotation at all             932   the gap that costs calls

The 932 is the worklist. It is a selection on
`release/02_reference/data/curation/screen_assemblies.tsv`, not a copied file, so it stays
correct when the assembly summary is refreshed:

    species_annotated == no  AND  genus_annotated == 0,  best assembly per species

751 of them (80.6%) have exactly one assembly. Best level is Scaffold 620, Contig 230,
Chromosome 64, Complete Genome 18, and 238 carry more than 5,000 contigs.

## Why the 2,604 are excluded and the 932 are not

An absent species whose genus is already annotated has a k-mer neighbour, so its reads land on
a near relative: wrong, but bounded and visible. An absent species whose whole genus is
missing has no neighbour, and its reads go somewhere arbitrary. The worked case is *Puccinia
hordei*, absent from db_v3 and read as *P. triticina*.

The same point in reverse is why the 2,604 are not worth predicting: the unannotated mass
concentrates in *Fusarium* (189 species), *Aspergillus* (144), *Penicillium* (93),
*Colletotrichum* (51) and *Alternaria* (45), genera the database is already deep in.

## Which tool, and why the split is forced

**miniprot** projects an annotated congener's proteins onto the target. It was chosen over ab
initio deliberately: in a competitive index an invented gene corrupts every species'
apportionment, so precision beats recall. It needs an annotated congener as a guide.

**vipsania** is an unsupervised ab-initio finder trained on raw sequence, needing only a FASTA.
Fungi locus-level F1 is 0.660, which reads badly until you ask the right question. A Kraken2
CDS reference does not need correct gene models, it needs CDS-like sequence from the right
organism: a misplaced exon boundary still yields real sequence whose k-mers point at the right
taxon. On *P. nodorum* GCA_016801405.1 vipsania produced 13,634 CDS against the true 17,764
(77%) but 19.29 Mb against the true 20.55 Mb (94%). Mediocre by gene count, excellent by bulk.

**All 932 are in genera with no annotation, so none has a congener to guide miniprot.** This
whole worklist is vipsania's case. miniprot stays relevant only for targets inside the 2,604.

## Predicted CDS is a competitor, never evidence

Homology projection biases toward genes conserved with the guide, and an ab-initio model is
not a species claim either. Predicted CDS may absorb reads; it must never support a positive
call about the species it came from. Every predicted species therefore dilutes species-level
resolution for the species that CAN be reported and buys nothing reportable itself, which is
the argument for the smallest set that closes real gaps rather than for coverage.

## Two databases, measured against each other

The design is a baseline of submitter annotations only and a second adding predicted CDS, not
one database with both mixed in. `kraken2-build --add-to-library` is incremental, so the
second reuses the first's library and only the predicted CDS is added. Without the baseline,
a change in the calls cannot be attributed to prediction rather than to the scope change.

## State

`sample.py`, `fetch_genomes.py` and `compare.py` carry the completed benchmark: recall 0.963 /
precision 0.945 (n=60) and novelty 0.963 (n=59). A production run was never decided.

**They still read `02_classify/_layout.py` and its silver tables**, which is why `MODULE`
points back there. Their inputs must be repointed at `release/02_reference` before a
production run: the worklist is now `screen_assemblies.tsv`, not `ref_candidates.tsv`.
