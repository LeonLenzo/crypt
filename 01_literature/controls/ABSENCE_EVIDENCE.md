# Evidence that no usable control library exists, per species

An absence claim needs more than a database count. "No RNA-seq runs in ENA" is a statement
about two archives on one date; "no study has produced such data" is a statement about the
field, and only the second is worth putting in a paper. This file records what was checked, how,
and what the honest claim is for each species.

Method, 2026-10-26, three independent checks:

1. **Archive counts.** ENA portal `tax_tree(taxid) AND library_strategy="RNA-Seq"`, species node
   not strain node, cross-checked against NCBI SRA via eutils. Both agreed everywhere tested.
2. **Synonymy.** Every superseded teleomorph and anamorph name queried against ENA taxonomy to
   confirm it resolves to the same node rather than a separate one that the count would miss.
3. **Literature.** Undermind deep search over the published record, asking specifically whether
   transcriptome work exists that was never deposited, went to a regional archive such as GSA
   which neither ENA nor SRA indexes, or remains embargoed.

## Synonymy is not the explanation, checked explicitly

`Cochliobolus carbonum`, `Helminthosporium carbonum`, `Drechslera zeicola`, `Bipolaris carbonum`,
`Cochliobolus victoriae`, `Helminthosporium victoriae`, `Drechslera tritici-repentis`,
`Helminthosporium tritici-repentis`, `Mycosphaerella zeae-maydis` and
`Fusarium tricinctum var. sporotrichioides` all return **no separate ENA taxonomy node**: NCBI
folds them into the current species, so a `tax_tree` count on the current name already covers
the old literature's names. The zeros are not a naming artefact.

**One real separation was found.** *Cercospora zeina* is a distinct species (taxid 348901) with
**7 RNA-seq runs**, not a synonym of *C. zeae-maydis*. It is the sister grey leaf spot pathogen
and is the declared target in two of the cohort's own maize projects, so the *Cercospora*
cluster is better served than the *C. zeae-maydis* count alone suggests.

## *Cercospora zeae-maydis*: data were generated, and are not retrievable

This is the case that justifies doing the literature check at all, because the bare archive
count would have been misleading in both directions.

- **Axenic transcript resources exist, but predate RNA-seq.** 768 cDNA-subtraction clones from
  cercosporin-inducing against suppressing cultures (Shi02), and 27,551 Sanger ESTs from
  vegetative and sporulating cultures on V8 agar and 0.2x PDA at 24 C (Blu08), deposited in
  GenBank dbEST. Confirmed from the full texts: no plant tissue in the EST library material.
  These are real axenic fungal transcript data and they are **unusable here**: Sanger ESTs carry
  no depth, no library chemistry comparable to the cohort, and nothing that a k-mer saturation
  statistic can be computed from.
- **Modern RNA-seq exists and its accession does not resolve.** Li et al. 2026 (J Fungi 12:650,
  `10.3390/jof12090650`) performed time-course RNA-seq of *C. zeae-maydis* and *C. zeina* during
  maize B73 infection, growth chamber, 28 C, 14/10 h, single-isolate inoculation by needleless
  syringe, sampled 0/1/3/6 dpi, poly(A) selected, paired 150 bp, **and in vitro cultures of each
  isolate alongside**. Its data availability statement reads: "The raw in vitro and in planta
  RNA-seq reads generated for this study are available in the NCBI Sequence Read Archive under
  BioProject accession number PRJNA1463572."
  **PRJNA1463572 returns zero records in NCBI BioProject, zero runs in SRA, and zero runs in
  ENA** (checked 2026-10-26). The accession is private, withdrawn or mistyped.

So the defensible claim for this species is not "no transcriptome work exists". It is: *the only
modern RNA-seq of this pathogen is published with an accession that does not resolve, and the
only retrievable fungal transcript data are Sanger-era ESTs.* That is a stronger sentence than a
database count and it is checkable by a reviewer.

**Worth one email.** The study is otherwise an ideal control: single isolate, growth chamber,
defined timepoints, poly(A), and it carries both the in planta and the axenic arm for the same
isolate. If the authors release PRJNA1463572 it would cover the species we are least able to
control for. The absence of a mock-inoculated arm is the one weakness; the 0 dpi sample was
taken 30 minutes after inoculation and is a baseline, not a pathogen-free control.

## *Bipolaris zeicola*

Archive: 0 RNA-seq runs, 5 runs total, all genomic. The literature hit that looks like a
counter-example is not one: the cited transcriptome study profiles **maize leaf responses** to
infection rather than fungal transcripts (Liu15b), and the other fungal resource is a genome
(Xia22h). No fungal transcriptome study is established.

This matters more than the others: *B. zeicola* is called as an undeclared pathogen in **68
cohort runs**, at 22% median coverage in one maize project where every run carries it. It is the
most frequently detected species for which no control can currently be built.

## *Bipolaris victoriae*, *Fusarium sporotrichioides*

Archive: 0 RNA-seq runs each (3 and 5 runs total, genomic). The literature search found no
qualifying fungal transcriptome study, **and explicitly could not establish field-wide absence**.
State these as "no public RNA-seq data and none found in the literature", never as "none exists".

## The limit of this evidence, stated plainly

For every species here the retrieved records lack repository accessions, so absence from ENA and
SRA cannot by itself distinguish non-deposition from regional deposition, embargo or restricted
access. GSA in China indexes a substantial amount of fungal transcriptome work and is indexed by
neither archive, and it cannot be swept programmatically the way ENA can. Any claim in the
manuscript should be scoped to retrievability, which is what actually matters for a control:
**no control library can be built for these species from public data**, rather than no such data
were ever produced.
