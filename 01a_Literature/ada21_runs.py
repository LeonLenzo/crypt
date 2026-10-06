#!/usr/bin/env python3
"""Flatten Ada21's Supplementary Table S1 onto run accessions.

Table S1 (`studies/doi_10.1186_s12864-021-07488-3/12864_2021_7488_MOESM1_ESM.xlsx`) is the
richest per-sample table in the corpus: 1,025 rows covering eleven studies, carrying a
town-level `Location`, a `Country`, a date, the host species and variety, and crucially a
`Sample Type` of Field or Lab. It is the only source we have that states, per sample,
whether a plant was collected from a field or inoculated in a lab.

It cannot be joined directly. Its key is `Accession Code`, which is an ENA SECONDARY sample
accession (ERS...) for the five PRJEB studies and a run accession (SRR...) for the four
PRJNA ones. `runs.tsv` holds neither ERS; it holds `Run` and `BioSample` (SAMEA...). So this
script resolves ERS -> run through the ENA portal and emits a table keyed on `Run`, which
`joins.tsv` can then use with no regex at all.

Why a derived table rather than a cleverer join rule: the ERS -> SAMEA mapping lives in the
archive, not in the supplement, so any single-step rule would have to smuggle a network
lookup into the join engine. Writing the mapping down once, reproducibly, keeps the join
engine dumb and keeps the derivation auditable. `SOURCE.txt` in the study directory records
where both inputs came from.

Two things this script deliberately does NOT do:

  * It does not treat `Date Isolated` as a collection date for Lab rows. The sheet footnotes
    that column: pre-2013 isolates were inoculated onto wheat in the lab and the date shown
    is the initial isolation. Lab rows therefore get `date=""`. The split is checked, not
    assumed: every Lab row must be pre-2013 and every Field row 2013 or later, and the script
    refuses to write if that fails.
  * It does not decide what a row means by itself. `Setting`, `SamplingSelection` and
    `Tissue` are emitted as columns because the discriminator is `SampleType`, which lives
    only in this table and so cannot be written as a `curation.tsv` predicate over runs.tsv.
    The mapping from `SampleType` to those values is SAMPLE_TYPE_MAP below, and the quote
    justifying it is the `evidence_id` on each join rule in `joins.tsv`.

Run:  python 01a_Literature/ada21_runs.py            (uses the cached ENA map if present)
      python 01a_Literature/ada21_runs.py --refetch   (re-queries the ENA portal)
"""

from __future__ import annotations

import argparse, csv, re, sys, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
STUDY = HERE / "studies" / "doi_10.1186_s12864-021-07488-3"
S1 = STUDY / "12864_2021_7488_MOESM1_ESM.xlsx"
ENA_CACHE = STUDY / "ena_run_sample_map.tsv"
OUT = STUDY / "ada21_runs.csv"

# The five ENA studies key Table S1 by ERS; the four SRA studies key it by run accession.
ENA_PROJECTS = ["PRJEB39201", "PRJEB33109", "PRJEB31334", "PRJEB15280", "PRJEB12497"]
PORTAL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
          "?accession={}&result=read_run&format=tsv"
          "&fields=run_accession,sample_accession,secondary_sample_accession,sample_title")


# What a Sample Type means for the cohort's own vocabulary: (setting, sampling_selection).
#
# `disease-selected` is the point of this table. These leaves were collected BECAUSE they
# were showing yellow rust, which is exactly the condition the cryptic co-infection
# hypothesis wants to test (the authors reported Pst and nothing else, so any second
# pathogen is unreported). These are therefore the cohort's PRIMARY co-infection samples:
# the question asked of them is not "is anything here" but "is anything here besides Pst".
# They are not a prevalence sample, which is a different question, not a lesser one.
# Neither `setting` nor `tissue` can express that, which is why the field exists.
#
# `inoculated` is a third case, not a shade of the other two: a Lab row is an archived
# isolate deliberately put onto wheat, so disease is present by the experimenter's design
# rather than by the sampler's choice. Same evidence standard as `disease-selected` - a
# finding is a pathogen beyond the inoculum - but with a susceptibility confound on top.
SAMPLE_TYPE_MAP = {
    "Field": ("field", "disease-selected"),
    "Lab": ("growth chamber (lab inoculation)", "inoculated"),
}


def fetch_ena_map() -> list[dict]:
    """One request per project, not per sample. See crypt-ena-portal-fetch."""
    rows = []
    for bp in ENA_PROJECTS:
        with urllib.request.urlopen(PORTAL.format(bp), timeout=120) as fh:
            text = fh.read().decode()
        lines = text.splitlines()
        hdr = lines[0].split("\t")
        for line in lines[1:]:
            if not line.strip():
                continue
            d = dict(zip(hdr, line.split("\t")))
            d["BioProject"] = bp
            rows.append(d)
        print(f"  {bp}: {len(lines) - 1} runs", file=sys.stderr)
    return rows


def load_ena_map(refetch: bool) -> list[dict]:
    if ENA_CACHE.exists() and not refetch:
        return list(csv.DictReader(ENA_CACHE.open(), delimiter="\t"))
    rows = fetch_ena_map()
    cols = ["BioProject", "run_accession", "sample_accession",
            "secondary_sample_accession", "sample_title"]
    with ENA_CACHE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"cached {len(rows)} rows -> {ENA_CACHE.relative_to(ROOT)}", file=sys.stderr)
    return rows


def load_s1() -> list[dict]:
    import openpyxl
    ws = openpyxl.load_workbook(S1, read_only=True)["Sheet1"]
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    hdr = [str(c).strip() if c else "" for c in rows[0]]
    # The sheet's final row is the footnote, which lands in the Study column; it is the only
    # row whose Study text runs past 80 characters.
    return [dict(zip(hdr, r)) for r in rows[1:] if r[0] and len(str(r[0])) < 80]


def norm_date(v) -> str:
    """Table S1 writes 2018/05/24, or 2012/11, or 2006. Keep whatever precision is given."""
    s = str(v).strip()
    if not s or s.lower() in ("none", "unknown", "na"):
        return ""
    s = s.split(" ")[0].replace("/", "-")
    m = re.match(r"^(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?$", s)
    if not m:
        return ""
    y, mo, d = m.groups()
    if d:
        return f"{y}-{int(mo):02d}-{int(d):02d}"
    return f"{y}-{int(mo):02d}" if mo else y


def year(s: str) -> int | None:
    m = re.match(r"^(\d{4})", s or "")
    return int(m.group(1)) if m else None


def place(rec: dict) -> str:
    """`Country: Location`, matching the cohort's existing location spelling. Where Location
    repeats the country or is unknown, the country alone is the honest answer."""
    loc = str(rec.get("Location") or "").strip()
    cty = str(rec.get("Country") or "").strip()
    if not cty or cty.lower() in ("none", "unknown"):
        return ""
    if not loc or loc.lower() in ("none", "unknown", "na") or loc.lower() == cty.lower():
        return cty
    return f"{cty}: {loc}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refetch", action="store_true", help="re-query the ENA portal")
    args = ap.parse_args()

    s1 = load_s1()
    print(f"Table S1: {len(s1)} sample rows", file=sys.stderr)

    # The date column means different things per Sample Type, so verify the footnote's claim
    # before relying on it. A violation means the sheet changed and the gate below is unsafe.
    bad = [(str(r["Sample IDs"]), str(r["Sample Type"]), norm_date(r["Date Isolated1"]))
           for r in s1
           if (y := year(norm_date(r["Date Isolated1"]))) is not None
           and ((str(r["Sample Type"]) == "Lab" and y >= 2013)
                or (str(r["Sample Type"]) == "Field" and y < 2013))]
    if bad:
        sys.exit(f"REFUSED: {len(bad)} rows contradict the footnote's pre-2013 = Lab split, "
                 f"e.g. {bad[:4]}. The date column can no longer be gated on Sample Type.")
    print("  footnote check: every Lab row pre-2013, every Field row 2013+", file=sys.stderr)

    by_acc = {str(r["Accession Code"]).strip(): r for r in s1}

    ena = load_ena_map(args.refetch)
    out, unmatched = [], []
    for e in ena:
        rec = by_acc.get(e["secondary_sample_accession"].strip())
        if rec is None:
            unmatched.append(e["run_accession"])
            continue
        out.append((e["run_accession"], rec, "ERS"))
    # SRA studies: Table S1 keys those rows by run accession directly.
    for acc, rec in by_acc.items():
        if acc.startswith(("SRR", "ERR")):
            out.append((acc, rec, "run"))

    rows = []
    for run, rec, via in out:
        st = str(rec["Sample Type"]).strip()
        dt = norm_date(rec["Date Isolated1"])
        setting, selection = SAMPLE_TYPE_MAP.get(st, ("", ""))
        host = str(rec.get("Host Species") or "").strip()
        rows.append(dict(
            Run=run,
            SampleID=str(rec["Sample IDs"]).strip(),
            AccessionCode=str(rec["Accession Code"]).strip(),
            JoinedVia=via,
            Study=str(rec["Study"]).strip(),
            SampleType=st,
            Setting=setting,
            SamplingSelection=selection,
            # Every sample in the compendium is infected leaf: "Pst-infected wheat leaf
            # samples were collected and initially stored in RNAlater". Left blank where the
            # sheet names no host, which happens only on Lab rows.
            Tissue="leaf" if host and host.lower() not in ("none", "unknown") else "",
            # Gated on Sample Type, per the footnote. Lab rows keep their isolation date in
            # IsolationDate so nothing is lost, but CollectionDate stays empty.
            CollectionDate=dt if st == "Field" else "",
            IsolationDate=dt if st != "Field" else "",
            Location=place(rec),
            Country=str(rec.get("Country") or "").strip(),
            HostSpecies=str(rec.get("Host Species") or "").strip(),
            HostVariety=str(rec.get("Host Variety") or "").strip(),
        ))
    rows.sort(key=lambda r: r["Run"])

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    n_f = sum(1 for r in rows if r["SampleType"] == "Field")
    locs = {r["Location"] for r in rows if r["SampleType"] == "Field" and r["Location"]}
    print(f"\nwrote {OUT.relative_to(ROOT)}: {len(rows)} runs "
          f"({n_f} Field, {len(rows) - n_f} other)", file=sys.stderr)
    print(f"  field localities: {len(locs)} distinct, "
          f"{len({r['Country'] for r in rows if r['SampleType'] == 'Field'})} countries",
          file=sys.stderr)
    if unmatched:
        print(f"  {len(unmatched)} ENA runs absent from Table S1, e.g. {unmatched[:5]}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
