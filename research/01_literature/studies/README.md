# studies/ — one directory per paper, named by DOI

    studies/doi_<DOI with / replaced by _>/

So `10.1038/s41467-019-09287-7` lives at `studies/doi_10.1038_s41467-019-09287-7/`.

**Always the DOI, never the paper_ref.** Undermind-style refs (`Cai21b`, `Tom26`, `Man24b`)
are convenient in a worklist but they are not stable identifiers: the same paper can carry
different refs across search reports, and a ref tells a reader nothing. Directories named by
ref caused a real failure on 2026-10-06, when a DOI-keyed audit reported `Cai21b`'s source as
missing while `fgene-12-716821.pdf` was sitting inside it. Migrated to DOI naming the same
day, on Leon's call.

The `dir` column of `data/joins.tsv` points at these names, so **renaming a directory means
updating that column**. Check before moving anything:

    awk -F'\t' 'NR>1 && $3!=""{print $3}' data/joins.tsv | sort -u

Each directory should hold:

    SOURCE.txt     the citation, the exact retrieval URL, the stored filename, and which
                   evidence_id the file backs. Tracked in git.
    <the source>   full text and supplements. NOT tracked: studies/** is gitignored and the
                   data policy keeps published material local.

`apply_curation.py` warns on every run if a DOI carrying evidence rows has no directory with
content, so a quote can always be checked against the paper it is attributed to.

For preprints, record the VERSION in SOURCE.txt and store the version you actually read.
Versions differ materially: bioRxiv 10.1101/2025.05.29.656841 v1 and v2 contain no field arm
at all, while v3 adds one plus eight accessions.
