#!/usr/bin/env python3
"""scaffold.py — a working directory per paper, and a worklist saying what each one needs.

The manual step in this module is unavoidable: supplements are shaped differently in every
journal and no amount of code makes them uniform. What code CAN do is make the manual step
small, ordered, and recorded. This builds the place where that work happens.

For every paper it creates `studies/<paper_ref>/` holding:

    NEEDED.md       what is missing, which BioProjects it affects, how many runs are at
                    stake, and which SRA identifier columns are available to join on
    adapter.yaml    the declarative config to fill in. Writing one is the manual act; it
                    names a file, a sheet, a key column, and which columns carry tissue,
                    site and date. No code is written per study.

Then `worklist.tsv`, ordered by runs affected, so the expensive attention goes where the data
is. A paper's need is derived from the triage state of the BioProjects linked to it:

    none            every linked BioProject is sra-complete. Nothing to do.
    supplement      at least one is sra-partial, few-biosamples or one-biosample.
    accession       at least one resolves to no runs, so the accession is wrong and must be
                    fixed before a supplement can help. Undermind cited PRJNA42371, which
                    NCBI resolves to *Amphibacillus xylanus*.
    no-accession    the paper names no BioProject at all.

**The join key is reported, not guessed.** For each BioProject, NEEDED.md lists which SRA
identifier columns are unique per run and therefore usable as a join target. PRJNA383416's
supplement joined on `LibraryName` at 1,960/1,960, and nothing in BioSample would have
revealed that. Knowing the candidates up front is what turns "work out how this table links"
into "check which column matches".

Directories are a drop zone, not the record. The authoritative tables stay flat in `data/`,
because a question like "how many field runs across all studies" must not become a loop over
300 directories. Nothing here is committed: these hold publisher PDFs and spreadsheets.

Usage:

    python 01a_Literature/curate/scaffold.py              # create or refresh every paper directory
    python 01a_Literature/curate/scaffold.py --needs-only # only papers that need something
"""
import argparse, collections, csv, json, re, sys
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _paths import ROOT
from _layout import DATA, MODULE, STUDIES


# Triage states that a supplement could fix, and the ones it cannot.
_NEEDS_SUPP = {"sra-partial", "few-biosamples", "one-biosample", "per-run-biosamples"}
_COMPLETE   = {"sra-complete"}

ADAPTER_TEMPLATE = """\
# Adapter for {ref}
#
# Fill this in once, by hand, after reading the paper's methods and opening its supplement.
# Nothing here is code: it names a file, a sheet, and which columns carry what.
#
# Delete `blocked:` when you fill the rest in. Leave a reason there if the paper genuinely
# cannot be resolved, so it stops reappearing in the worklist.
blocked: not yet attempted

paper_ref: {ref}
doi: {doi}
bioprojects: [{bps}]

# The supplement file, relative to this directory, and the sheet holding one row per sequenced
# sample. Drop the file in beside this adapter.
supplement:
  file:
  sheet:

# The column in that sheet whose values match an SRA identifier, and which identifier it is.
# Candidates measured for these BioProjects (unique per run, so usable as a join target):
{candidates}
#
# CHECK THIS BEFORE ANYTHING ELSE. If no column matches any SRA identifier, the table is not
# describing the sequenced samples and must not be merged. PRJNA1217477's PNAS supplement
# tabulates the origins of the inoculum ISOLATES, not the glasshouse plants that were
# sequenced; merging it would have stamped Californian vineyard coordinates onto Arabidopsis.
join:
  column:
  matches:          # one of: Run, Experiment, LibraryName, SampleName, BioSample

# Columns carrying the fields we need. Leave blank where the supplement does not have it.
columns:
  tissue:
  location:
  collection_date:
  setting:

# Where a column holds codes rather than values, map them here, and say where the mapping came
# from. An entry without evidence is a guess, and guesses here are indistinguishable from data.
# Worked example, PRJNA383416 (Kremling et al. 2018, doi:10.1038/nature25966):
#   value_map:
#     tissue:
#       GRoot:  {{value: "germinating seedling root",  setting: "growth chamber",
#                evidence: "Methods, Tissue collection: walk-in growth chamber, vermiculite"}}
#       LMAD:   {{value: "mature leaf, day",           setting: "field",
#                evidence: "Methods: Musgrave Research Farm, collected 11:00-13:00"}}
value_map:
"""


def dir_name(paper: dict) -> str:
    """Filesystem-safe directory name for a paper. ALWAYS the DOI where one exists.

    The DOI comes first on purpose, changed 2026-10-06 (leon's call). This function used to
    prefer the Undermind reference because it is short and readable, but a ref is not a stable
    identifier: the same paper carries different refs across search reports, and a ref tells a
    reader nothing. It also broke an audit - a DOI-keyed check for "does this evidenced paper
    have a local source" reported `Cai21b` as missing while its PDF sat inside it. See
    studies/README.md.

    Never return an empty string: `STUDIES / ""` resolves to `studies/` itself, and 182 papers
    without a reference each overwrote a NEEDED.md there instead of getting a directory.
    """
    doi = (paper.get("doi") or "").strip()
    if doi:
        return "doi_" + re.sub(r"[^A-Za-z0-9._-]", "_", doi)
    ref = (paper.get("paper_ref") or "").strip()
    if ref:
        return re.sub(r"[^A-Za-z0-9._-]", "_", ref)
    key = (paper.get("paper_key") or "").strip()
    return re.sub(r"[^A-Za-z0-9._-]", "_", key) or "_unkeyed"


def read(name: str) -> list[dict]:
    path = DATA / name
    if not path.exists():
        sys.exit(f"{path} missing; run triage.py first")
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def join_candidates(runs: list[dict], attrs: dict | None = None) -> list[str]:
    """Columns a supplement could join on, with how well each discriminates.

    Two kinds, and missing the second is a real failure mode. The obvious kind is an SRA
    identifier. The other is a BIOLOGICAL ATTRIBUTE carried in the BioSample: genotype,
    cultivar, ecotype, isolate, line. Supplements are written by biologists about plants, so
    they key on the plant far more often than on an archive accession.

    PRJNA746402 is the case that forced this. Its supplement matched 0 of 450 runs on every
    SRA identifier, which by the usual rule means "this table does not describe the sequenced
    samples" and stop. It in fact joins perfectly on `ecotype`, 6 of 6 genotypes, and the run
    identifiers are an anonymous ZEMAYS1_38_1 scheme that could never have matched anything.

    An attribute is reported with its distinct count, because a genotype column is a
    many-to-one join: it carries study design down to groups of runs, not to single runs.
    """
    out = []
    for col in ("Run", "Experiment", "LibraryName", "SampleName", "BioSample"):
        vals = [r[col] for r in runs if r.get(col)]
        if not vals:
            continue
        n = len(set(vals))
        tag = "unique" if n == len(runs) else f"{n} distinct of {len(runs)}"
        out.append(f"{col} ({tag})")

    # Attribute values, from whatever BioSamples were fetched. Single-valued attributes are
    # dropped: a column that is "Spain" on every row joins everything to everything.
    for key in ("ecotype", "cultivar", "genotype", "isolate", "strain", "host",
                "tissue", "dev_stage", "treatment", "collected_by"):
        vals = [a[key] for a in (attrs or {}).values() if a.get(key)]
        n = len(set(vals))
        if n > 1:
            out.append(f"BioSample.{key} ({n} distinct values, many-to-one)")
    return out


def inbox_files() -> list[Path]:
    """Loose files sitting in studies/ itself, not yet filed against a paper.

    An inbox, deliberately. While a project is being worked, PDFs and supplements arrive
    faster than their owning paper is known, and guessing wrong is worse than waiting: five
    papers ended up in doi_10.1038_srep37302/ because that was the directory open at the
    time, and four of them had nothing to do with it. Drop files here, file them once the
    paper is identified.
    """
    if not STUDIES.is_dir():
        return []
    return sorted(f for f in STUDIES.iterdir()
                  if f.is_file() and ":Zone.Identifier" not in f.name
                  and not f.name.startswith("~$"))


def build_project_views(links: list[dict], dirs: dict) -> int:
    """studies/by-project/<PRJ>/ -> symlinks to every paper directory touching that project.

    The relationship is a graph and a directory tree cannot hold one: Ada21 names ten
    BioProjects, PRJNA306542 is cited by eleven papers. Files belong to PAPERS (a supplement
    is a property of a paper) and curation rules belong to PROJECTS, so the files live under
    the paper and this gives the other lookup without duplicating a byte.

    Regenerated each run, so a stale symlink cannot outlive the link it came from.
    """
    root = STUDIES / "by-project"
    if root.exists():
        for d in sorted(root.iterdir()):
            if d.is_dir():
                for l in d.iterdir():
                    l.unlink()
                d.rmdir()
    made = 0
    by_bp = collections.defaultdict(set)
    for l in links:
        if l["BioProject"] and l["paper_ref"] in dirs:
            by_bp[l["BioProject"]].add(dirs[l["paper_ref"]])
    for bp, paper_dirs in sorted(by_bp.items()):
        live = [d for d in sorted(paper_dirs) if (STUDIES / d).is_dir()]
        if not live:
            continue
        (root / bp).mkdir(parents=True, exist_ok=True)
        for d in live:
            link = root / bp / d
            if not link.exists():
                link.symlink_to(Path("../..") / d)
            made += 1
    return made


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--needs-only", action="store_true",
                    help="skip papers whose BioProjects are all sra-complete")
    ap.add_argument("--materialise", action="store_true",
                    help="also create the per-paper directories. Off by default: 423 stub "
                         "directories holding nothing but a blank adapter are noise, and "
                         "review.py creates each one on sight when you reach that paper.")
    ap.add_argument("--key", nargs="*", metavar="K", default=None,
                    help="materialise only these papers (by key, ref or directory name)")
    args = ap.parse_args()

    papers = read("papers.tsv")
    bps    = {r["BioProject"]: r for r in read("bioprojects.tsv")}
    links  = read("paper_bioproject.tsv")
    runs   = read("runs.tsv")

    runs_by_bp = collections.defaultdict(list)
    for r in runs:
        runs_by_bp[r["BioProject"]].append(r)
    attrs_by_bp = {}
    for f in (DATA / "biosample_attrs").glob("*.json") if (DATA / "biosample_attrs").is_dir() else []:
        try:
            obj = json.loads(f.read_text())
        except Exception:
            continue
        attrs_by_bp[f.stem] = obj.get("attrs", obj) if isinstance(obj, dict) else {}
    bps_by_ref = collections.defaultdict(list)
    for l in links:
        bps_by_ref[l["paper_ref"]].append(l)

    work = []
    for p in papers:
        # papers.tsv is keyed on paper_key (a DOI where one exists). The link table uses the
        # same key, so joining on paper_ref silently orphans every cohort paper: 302 of 423
        # read as no-accession when their BioProjects were right there.
        ref = p["paper_key"]
        mine = [bps[l["BioProject"]] for l in bps_by_ref.get(ref, [])
                if l["BioProject"] in bps]
        states = {b["triage"] for b in mine}
        if not mine:
            need = "no-accession"
        elif "no-runs" in states:
            need = "accession"
        elif states <= _COMPLETE:
            need = "none"
        elif states & _NEEDS_SUPP:
            need = "supplement"
        else:
            need = "review"

        affected = [b for b in mine if b["triage"] in _NEEDS_SUPP]
        runs_affected = sum(int(b["sra_runs"] or 0) for b in affected)
        unexamined = sum(int(b["runs_unexamined"] or 0) for b in mine)

        work.append(dict(
            paper_key=ref, paper_ref=p.get("paper_ref", ""), dir=dir_name(p),
            need=need, doi=p.get("doi", ""),
            n_bioprojects=len(mine), runs_affected=runs_affected,
            runs_unexamined=unexamined,
            bioprojects=";".join(b["BioProject"] for b in mine),
            states=";".join(sorted(states)),
            undermind_decision=p.get("undermind_decision", ""),
            relations=";".join(sorted({l["relation"] for l in bps_by_ref.get(ref, [])})),
            adapter="", supplement_present="",
        ))

    work.sort(key=lambda r: (-r["runs_affected"], -r["runs_unexamined"]))

    made = 0
    for w in work:
        if args.needs_only and w["need"] in ("none", "no-accession"):
            continue
        if args.key is not None and not ({w["paper_key"], w.get("paper_ref"), w["dir"]}
                                         & set(args.key)):
            continue
        # The directory is a drop zone for files a human downloads. Creating 423 of them
        # before anyone has downloaded anything just buries the one that matters.
        if not (args.materialise or args.key) and not (STUDIES / w["dir"]).exists():
            continue
        d = STUDIES / w["dir"]
        d.mkdir(parents=True, exist_ok=True)
        affected = [b for b in (bps[a] for a in w["bioprojects"].split(";") if a)
                    if b["triage"] in _NEEDS_SUPP]

        lines = [f"# {w['paper_key']} — what this study needs", "",
                 f"- **directory:** {w['dir']}",
                 f"- **need:** {w['need']}",
                 f"- **doi:** {w['doi'] or '(none recorded)'}",
                 f"- **Undermind decision:** {w['undermind_decision'] or '(none)'}",
                 f"- **runs affected:** {w['runs_affected']:,}",
                 f"- **runs never examined:** {w['runs_unexamined']:,}", ""]
        if not affected:
            lines += ["Nothing to resolve by hand. Every linked BioProject has per-run "
                      "BioSamples carrying both geography and a collection date.", ""]
        for b in affected:
            rr = runs_by_bp.get(b["BioProject"], [])
            lines += [f"## {b['BioProject']} — {b['triage']}", "",
                      f"- organism: {b['organism']}",
                      f"- runs: {b['sra_runs']} in SRA, {b['sra_biosamples']} BioSamples, "
                      f"{b['gate_passed']} passed the old STAT gate",
                      f"- BioSample coverage: geography {b['geo_pct'] or '?'}%, "
                      f"date {b['date_pct'] or '?'}%", "",
                      "Join targets available in SRA for this project:"]
            lines += [f"  - {c}" for c in join_candidates(rr, attrs_by_bp.get(b['BioProject']))] or ["  - (no runs cached)"]
            if b["triage"] == "one-biosample":
                lines += ["", "> All runs share a single BioSample, so SRA holds no per-run "
                          "metadata at all. A supplement is mandatory, and the join will have "
                          "to go through LibraryName or Experiment."]
            lines.append("")
        (d / "NEEDED.md").write_text("\n".join(lines))

        ad = d / "adapter.yaml"
        if not ad.exists():
            cand = []
            for b in affected:
                for c in join_candidates(runs_by_bp.get(b["BioProject"], [])):
                    cand.append(f"#   {b['BioProject']}: {c}")
            ad.write_text(ADAPTER_TEMPLATE.format(
                ref=w["paper_key"], doi=w["doi"],
                bps=", ".join(w["bioprojects"].split(";")) if w["bioprojects"] else "",
                candidates="\n".join(cand) or "#   (none; no BioProject needs a supplement)"))
        w["adapter"] = str(ad.relative_to(MODULE))
        w["supplement_present"] = "yes" if any(
            f.suffix.lower() in (".xls", ".xlsx", ".csv", ".tsv", ".txt")
            for f in d.iterdir() if f.is_file()) else ""
        made += 1

    dirs = {w["paper_key"]: w["dir"] for w in work}
    n_links = build_project_views(links, dirs)
    if n_links:
        print(f"  by-project/: {n_links} symlinks across "
              f"{len(list((STUDIES / 'by-project').iterdir()))} projects")
    loose = inbox_files()
    if loose:
        print(f"\n  INBOX: {len(loose)} unfiled file(s) in studies/ — "
              f"file them once their paper is known")
        for f in loose[:10]:
            print(f"    {f.name[:70]}")

    out = DATA / "worklist.tsv"
    with open(out, "w", newline="") as fh:
        w8 = csv.DictWriter(fh, fieldnames=list(work[0]), delimiter="\t")
        w8.writeheader()
        w8.writerows(work)

    note = "" if (args.materialise or args.key) else \
        "  (existing only; --materialise or --key to create more)"
    print(f"{len(work)} papers; {made} directories written under "
          f"{STUDIES.relative_to(ROOT)}/{note}")
    for need, n in collections.Counter(w["need"] for w in work).most_common():
        ra = sum(w["runs_affected"] for w in work if w["need"] == need)
        print(f"  {need:<14}{n:>4} papers{ra:>9,} runs affected")
    print(f"\ntop of the queue:")
    print(f"  {'paper':<24}{'need':<13}{'affected':>9}{'unexam':>8}  bioprojects")
    for w in work[:12]:
        if w["need"] == "none":
            continue
        print(f"  {w['dir'][:22]:<24}{w['need']:<13}{w['runs_affected']:>9,}"
              f"{w['runs_unexamined']:>8,}  {w['bioprojects'][:48]}")
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
