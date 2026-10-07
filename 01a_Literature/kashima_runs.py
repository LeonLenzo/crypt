#!/usr/bin/env python3
"""Rebuild Kashima's collection dates at full precision, in local time.

PRJDB7234 is the largest study in the field cohort (1,557 runs) and the only one that
repeatedly samples a single field across a season: sixteen sets of bihourly sampling over
24 h, May to September 2015 at Takatsuki plus a 2016 set at Kizugawa, 35 calendar days in
all. The date is the independent variable for anything asking how a community accumulates,
so collapsing it loses the point of the dataset.

It was collapsed. The 2026-10-05 join rules `kashima-year` read the supplement's `year`
column alone, writing "2015" or "2016" over an archive field that already held a full
timestamp for all 1,557 runs. The 26 runs the join missed are the only ones that kept it,
which is how the loss became visible at all.

Two sources, and the supplement is the better one:

    archive (ENA portal)   2015-08-03T13:00:00Z        UTC, complete, hour resolution
    supplement TableS2     year/month/day/hour, JST    local, complete, hour resolution

They agree (the supplement's JST hour converts exactly to the SRA UTC stamp), but UTC is the
wrong frame for a diel series: 23:00 UTC is 08:00 the NEXT day in JST, so a UTC stamp both
shifts the hour away from its solar meaning and can move a sample to the wrong calendar day.
This script emits the local time with an explicit +09:00 offset, which keeps the solar hour
readable and stays convertible.

Output is keyed on `Run` so `joins.tsv` needs no regex, the same shape as `ada21_runs.py`.
The supplement is keyed on `sampleID`, a 5-digit integer that appears at the end of the
BioSample description, so the sampleID -> Run mapping comes from the ENA portal's
`sample_title`.

Run:  python 01a_Literature/kashima_runs.py
      python 01a_Literature/kashima_runs.py --refetch
"""

from __future__ import annotations

import argparse, csv, re, sys, urllib.request
from pathlib import Path

from _layout import STUDIES

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
STUDY = STUDIES / "doi_10.1093_pcp_pcab088"
SUPP = STUDY / "Copy of pcp-2021-e-00065-File009.xlsx"
SHEET, HEADER_ROW = "TableS2", 7          # zero-based: the header sits on sheet row 8
ENA_CACHE = STUDY / "ena_run_sample_map.tsv"
OUT = STUDY / "kashima_runs.csv"

PORTAL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
          "?accession=PRJDB7234&result=read_run&format=tsv"
          "&fields=run_accession,sample_accession,sample_title,collection_date")
JST = "+09:00"


def load_portal(refetch: bool) -> list[dict]:
    if ENA_CACHE.exists() and not refetch:
        return list(csv.DictReader(ENA_CACHE.open(), delimiter="\t"))
    with urllib.request.urlopen(PORTAL, timeout=180) as fh:
        lines = fh.read().decode().splitlines()
    hdr = lines[0].split("\t")
    rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]
    with ENA_CACHE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"cached {len(rows)} rows -> {ENA_CACHE.relative_to(ROOT)}", file=sys.stderr)
    return rows


def sample_id(title: str) -> str | None:
    """The supplement's sampleID is the trailing integer of the sample title. Two spellings
    exist for the same value and both end the same way:
        EBI BioSamples : "... Sample ID: 20001"
        ENA portal     : "The youngest fully expanded leaf Oryza sative (SL1207):20001"
    """
    m = re.search(r"(\d{4,6})\s*$", title or "")
    return m.group(1) if m else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refetch", action="store_true")
    args = ap.parse_args()

    import openpyxl
    ws = openpyxl.load_workbook(SUPP, read_only=True)[SHEET]
    rows = list(ws.iter_rows(min_row=HEADER_ROW + 1, values_only=True))
    hdr = [str(c).strip() if c else "" for c in rows[0]]
    supp = {}
    for r in rows[1:]:
        rec = dict(zip(hdr, r))
        if rec.get("sampleID") in (None, ""):
            continue
        sid = str(rec["sampleID"]).split(".")[0].strip()
        supp[sid] = rec
    print(f"supplement TableS2: {len(supp)} samples", file=sys.stderr)

    portal = load_portal(args.refetch)
    print(f"portal: {len(portal)} runs", file=sys.stderr)

    def ival(v):
        try:
            return int(float(str(v).strip()))
        except (TypeError, ValueError):
            return None

    out, miss, disagree = [], 0, []
    for p in portal:
        sid = sample_id(p.get("sample_title", ""))
        rec = supp.get(sid) if sid else None
        if rec is None:
            miss += 1
            continue
        y, mo, d, h = (ival(rec.get(k)) for k in ("year", "month", "day", "hour"))
        if None in (y, mo, d):
            miss += 1
            continue
        local = f"{y:04d}-{mo:02d}-{d:02d}"
        if h is not None:
            local += f"T{h:02d}:00:00{JST}"
        # Cross-check the supplement's local day against the archive's UTC stamp. A JST day
        # and a UTC day differ by at most one, so anything further apart means the join is
        # wrong, not merely offset.
        arc = (p.get("collection_date") or "")[:10]
        if arc and abs((int(arc[:4]) - y) * 365 + (int(arc[5:7]) - mo) * 30
                       + (int(arc[8:10]) - d)) > 1:
            disagree.append((p["run_accession"], arc, local))
        out.append(dict(
            Run=p["run_accession"], SampleID=sid,
            CollectionDateLocal=local,
            CollectionDateArchiveUTC=p.get("collection_date", ""),
            Site=str(rec.get("field") or "").strip(),
            LineName=str(rec.get("LineName") or "").strip(),
            TransPlantSet=str(ival(rec.get("TransPlantSet")) or "").strip(),
            QC=str(rec.get("QC") or "").strip(),
        ))

    if disagree:
        sys.exit(f"REFUSED: {len(disagree)} runs where the supplement's local day and the "
                 f"archive's UTC day differ by more than one day, e.g. {disagree[:4]}. "
                 f"The sampleID join is not aligning the right samples.")
    print("  day cross-check: every supplement day within 1 of the archive's UTC day",
          file=sys.stderr)

    out.sort(key=lambda r: r["Run"])
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    days = {r["CollectionDateLocal"][:10] for r in out}
    hours = {r["CollectionDateLocal"][11:16] for r in out if len(r["CollectionDateLocal"]) > 10}
    print(f"\nwrote {OUT.relative_to(ROOT)}: {len(out)} runs", file=sys.stderr)
    print(f"  {len(days)} calendar days, {len(hours)} distinct local hours", file=sys.stderr)
    print(f"  sites: {sorted({r['Site'] for r in out})}", file=sys.stderr)
    if miss:
        print(f"  {miss} portal runs had no matching sampleID", file=sys.stderr)


if __name__ == "__main__":
    main()
