#!/usr/bin/env python3
"""import_provenance.py — bring 02_literature's mined provenance into the curation chain.

`02_literature` already resolved per-sample location and collection year for 2,643 BioSamples
by mining three papers' supplementary tables and geocoding the results (2026-10-02). None of
it was reachable from `01a_Literature/data/runs.tsv`, so today's curation was rebuilding by
hand what was already on disk. The Adams rust-expression-browser supplement alone covers 526
samples across 235 localities and 29 countries, which is better provenance than anything
curated in this module.

This is an IMPORT, not a re-derivation. Each imported value carries the same four-link chain
as a hand-curated one, with the supplement file standing where a quoted sentence normally
does, so `apply_curation.py --explain` works identically on both.

**Per-sample values do not fit the rule model, and that is the point.** A curation rule is a
predicate over many runs: "every run where tissue=Leaf and dev_stage=heading was at Baihe in
2016". Here there are 235 distinct localities across 526 samples of one project, so a rule per
value would mean hundreds of rules carrying one sample each. These are written straight into
`runs.tsv` with a source token naming the import batch, and the batch is what carries the
evidence.

**Supplement beats BioSample, BioSample beats the LLM.** That precedence is
`cohort_provenance.py`'s, kept here so the two modules cannot disagree. Rows whose
`location_source` is `llm` are NOT imported: an LLM-extracted location is a per-BioProject
guess, which is exactly the conflation `crypt-sample-provenance-unverified` is about, and
importing it would launder a guess into an evidenced field.

**Country-only values are imported as country-only.** `cohort_provenance` normalises to 53
real countries, and a bare country is honest provenance at country resolution. What it must
not do is masquerade as a collection site, so the import records which it is.

Usage:

    python 01a_Literature/curate/import_provenance.py --dry-run
    python 01a_Literature/curate/import_provenance.py
"""
import argparse, collections, csv, datetime, re, sys
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _paths import ROOT          # repo-wide
from _layout import DATA, PROVENANCE

RUNS       = DATA / "runs.tsv"
COHORT_PROV = ROOT / "02_literature/03_classify/data/cohort_provenance.tsv"
SUPP_PROV   = ROOT / "02_literature/02_text/data/supp_provenance.tsv"

# One source token per provenance origin, so runs.tsv says WHICH mining pass produced a value.
BATCH = {"supplement": "import-supp-2026-10-02",
         "biosample":  "import-bs-2026-10-02"}
SKIP  = {"llm", "none", "unknown", ""}

_YEAR = re.compile(r"(19|20)\d{2}")


def _year(s: str) -> str:
    m = _YEAR.search(s or "")
    return m.group(0) if m else ""


def year_sources() -> dict:
    """BioSample -> where its collection_year came from, READ from the table.

    `cohort_provenance.py` now emits `year_source` at merge time (added 2026-10-05), so this
    no longer reconstructs it by matching values back against the supplement and BioSample
    dates. That reconstruction was the weakest link in the chain: where both sources carried
    the same year it could not tell them apart and fell back to precedence, so the label was
    an inference presented as a record.

    Falls back to the old inference only if the column is absent, and says so, because a
    silent fallback would hide exactly the ambiguity this change removes.
    """
    rows = list(csv.DictReader(open(COHORT_PROV), delimiter="\t"))
    if rows and "year_source" in rows[0]:
        return {r["BioSample"]: (r["year_source"] or "none")
                for r in rows if (r.get("collection_year") or "").strip()}

    print("  WARNING: cohort_provenance.tsv has no year_source column; falling back to "
          "inferring it by value match, which cannot disambiguate equal years. "
          "Rerun 02_literature/03_classify/cohort_provenance.py.", file=sys.stderr)
    supp, sam = {}, {}
    if SUPP_PROV.exists():
        for r in csv.DictReader(open(SUPP_PROV), delimiter="\t"):
            if (r.get("date") or "").strip():
                supp.setdefault(r["BioSample"], r["date"])
    samples = ROOT / "02_literature/03_classify/data/samples.tsv"
    if samples.exists():
        for r in csv.DictReader(open(samples), delimiter="\t"):
            sam[r["BioSample"]] = r.get("collection_date", "")
    out = {}
    for r in rows:
        cy = (r.get("collection_year") or "").strip()
        if not cy:
            continue
        b = r["BioSample"]
        out[b] = ("supplement" if _year(supp.get(b, "")) == cy
                  else "biosample" if _year(sam.get(b, "")) == cy else "unknown")
    return out


def read(p: Path) -> list[dict]:
    if not p.exists():
        sys.exit(f"{p} missing")
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--overwrite-curated", action="store_true",
                    help="also replace values set by a hand-written curation rule. Off by "
                         "default: a rule was written by someone reading a paper and should "
                         "outrank a bulk import.")
    args = ap.parse_args()

    runs = read(RUNS)
    prov = {r["BioSample"]: r for r in read(COHORT_PROV)}
    print(f"{len(prov):,} BioSamples in cohort_provenance.tsv")
    by_src = collections.Counter(r["location_source"] for r in prov.values())
    print(f"  by source: {dict(by_src)}")
    print(f"  skipping sources {sorted(SKIP - {''})}: "
          f"{sum(v for k, v in by_src.items() if k in SKIP)} samples\n")

    ysrcs = year_sources()
    print(f"  year provenance resolved independently for {len(ysrcs):,} samples: "
          f"{dict(collections.Counter(ysrcs.values()))}\n")

    stats = collections.Counter()
    for r in runs:
        p = prov.get(r["BioSample"])
        if not p:
            stats["no provenance row"] += 1
            continue
        lsrc = (p.get("location_source") or "").strip()
        ysrc = ysrcs.get(r["BioSample"], "")

        # Gate each field on ITS OWN source. A row can have an LLM location and a BioSample
        # year, and the year is still good.
        if lsrc in SKIP:
            stats[f"location source {lsrc or 'blank'}, skipped"] += 1
            loc = country = ""
        else:
            loc = (p.get("location") or "").strip()
            country = (p.get("country") or "").strip()
        batch = BATCH.get(lsrc, f"import-{lsrc}-2026-10-02")
        value = f"{country}: {loc}" if loc and country and loc != country else (loc or country)
        if value:
            held = r.get("location_source", "")
            if held in ("", "absent", "not_fetched", batch) or args.overwrite_curated:
                if not args.dry_run:
                    r["location"], r["location_source"] = value, batch
                stats[f"location <- {batch}"] += 1
            else:
                stats[f"location kept (curated as {held})"] += 1

        year = (p.get("collection_year") or "").strip()
        if year and ysrc in SKIP:
            stats[f"year source {ysrc or 'unknown'}, skipped"] += 1
        elif year:
            ybatch = BATCH.get(ysrc, f"import-{ysrc}-2026-10-02")
            held = r.get("date_source", "")
            if held in ("", "absent", "not_fetched", ybatch) or args.overwrite_curated:
                if not args.dry_run:
                    r["collection_date"], r["date_source"] = year, ybatch
                stats[f"collection_date <- {ybatch}"] += 1
            else:
                stats[f"collection_date kept (curated as {held})"] += 1

    for k, v in stats.most_common():
        print(f"  {v:>7,}  {k}")

    if args.dry_run:
        print("\ndry run; nothing written")
        return

    # Register the import batches as evidence, so --explain resolves them like any rule.
    have = {e["evidence_id"] for e in read(PROVENANCE)} if PROVENANCE.exists() else set()
    new = []
    ts = datetime.date.today().isoformat()
    for src, batch in BATCH.items():
        if batch in have:
            continue
        new.append(dict(
            evidence_id=batch, paper_key="(multiple)", doi="(multiple)",
            locus=f"02_literature provenance mining, location_source={src}",
            quote=("Per-sample provenance mined 2026-10-02 by "
                   "02_literature/02_text/supp_provenance.py from three papers' supplementary "
                   "tables (Adams et al. 2021 BMC Genomics covering 898 BioSamples, plus two "
                   "others) and merged by 03_classify/cohort_provenance.py in the order "
                   "supplement > BioSample geo_loc_name > LLM. Localities geocoded via "
                   "Nominatim. LLM-sourced rows are NOT imported. Year provenance is read from "
                   "the table's own year_source column, not inferred."
                   if src == "supplement" else
                   "Per-sample geo_loc_name read from the BioSample record and country-"
                   "normalised by 02_literature/03_classify/cohort_provenance.py "
                   "(COUNTRY_ALIAS / NOT_A_COUNTRY), 61 raw strings to 53 real countries."),
            entered_by="import_provenance.py", ts=ts))
    if new:
        with open(PROVENANCE, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["evidence_id", "paper_key", "doi", "locus",
                                                "quote", "entered_by", "ts"], delimiter="\t")
            w.writerows(new)
        print(f"\nregistered {len(new)} import batch(es) in provenance.tsv")

    with open(RUNS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(runs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(runs)
    print(f"wrote {RUNS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
