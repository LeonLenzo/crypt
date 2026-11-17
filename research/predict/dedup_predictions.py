#!/usr/bin/env python3
"""dedup_predictions.py — repair prediction outputs that two concurrent jobs wrote twice.

    python research/predict/dedup_predictions.py --check          # report only, touch nothing
    python research/predict/dedup_predictions.py --apply

Reads  research/predict/data/predicted/<acc>/{coding.fa,proteins.fa,genes.gff3}
Writes the same files, reduced to one copy, plus data/dedup_report.tsv
                                                   data/dedup_accessions.txt

## What happened, because the cause decides whether the repair is safe

Two array jobs ran the same list CONCURRENTLY on 2026-11-16: 53228691 submitted at 14:19
and 53231426 at 14:25, while the first was still running. `predict125.slurm` skips a target
whose `coding.fa` already exists, but vipsania writes that file only at the END of a task,
so for the 36 genomes still in flight at 14:25 the guard saw nothing and both jobs
predicted the same genome into the same directory. Vipsania opens its outputs for append,
so each of the three files came out holding two complete copies. BUSCO then refused all 36
with `ERROR: Duplicate of sequence >...`, which is a BUSCO input constraint and says
nothing about the quality of the predictions.

THE FILES ARE NOT ALWAYS TWO CONCATENATED COPIES, AND ASSUMING SO BREAKS 9 OF THE 36. Where
the first job finished writing before the second began, the file is exactly A+A and a
byte-halves test finds it. Where the two write phases OVERLAPPED, the copies are interleaved
in buffer-sized chunks: GCA_028456385.1's `proteins.fa` holds 27,614 records over 13,807
distinct ids, every id exactly twice, but no single split point separates them. So the
repair works at the level of RECORDS, not bytes, which covers both cases with one rule.

LOSSLESSNESS IS VERIFIED, NOT ASSUMED. A file is only rewritten when every record appears
exactly twice AND the two occurrences are byte-identical. Anything else — an odd
multiplicity, or two occurrences of one id whose content differs, which is what a write
spliced mid-record would look like — is reported and LEFT ALONE. There is one prediction
here, stored twice; this is not a dedup that has to choose between competing records.

Worth recording separately: two independent GPU runs over the same genome produced
byte-identical records, so vipsania's inference is deterministic on this hardware. The
October CPU-versus-GPU comparison saw that at the level of sequence COUNT; this confirms it
per record.

## Why not simply re-predict the 36

Nothing would be gained. The duplicate billing (~480 SU on fl3-gpu) is already spent, and
re-running would spend it a third time to arrive at the file this script produces. Both
copies completed exit 0; the prediction was never what failed.

## GFF3 record order

Records are emitted in order of first occurrence. For an A+A file that is exactly copy A.
For an interleaved file it is a mix of the two copies' positions, so the output is not
necessarily in contig order. Nothing downstream depends on GFF3 ordering — `coding.fa` is
what enters the Kraken2 library and the GFF3 is kept as a record — but it is stated rather
than left to be discovered.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

MODULE = Path(__file__).resolve().parent
PREDICTED = MODULE / "data" / "predicted"
REPORT = MODULE / "data" / "dedup_report.tsv"
ACCESSIONS = MODULE / "data" / "dedup_accessions.txt"

TARGETS = ("coding.fa", "proteins.fa", "genes.gff3")


def parse_fasta(raw: bytes) -> tuple[list[bytes], list[tuple[bytes, bytes]]]:
    """Return (leading junk, [(header, sequence block)]). Sequence keeps its newlines."""
    records: list[tuple[bytes, bytes]] = []
    header: bytes | None = None
    chunks: list[bytes] = []
    for line in raw.splitlines(keepends=True):
        if line.startswith(b">"):
            if header is not None:
                records.append((header, b"".join(chunks)))
            header, chunks = line, []
        elif header is not None:
            chunks.append(line)
    if header is not None:
        records.append((header, b"".join(chunks)))
    return [], records


def parse_gff(raw: bytes) -> tuple[list[bytes], list[tuple[bytes, bytes]]]:
    """Return (leading comment lines, [(line, b"")]). The whole line is the record key."""
    head: list[bytes] = []
    records: list[tuple[bytes, bytes]] = []
    seen_body = False
    for line in raw.splitlines(keepends=True):
        if line.startswith(b"#") and not seen_body:
            head.append(line)
            continue
        if line.strip():
            seen_body = True
            records.append((line, b""))
    return head, records


def classify(raw: bytes, name: str):
    """Decide what `raw` is. Returns (status, repaired bytes or None, n_before, n_after).

    single     every record occurs once — nothing to do
    doubled    every record occurs exactly twice with identical content — repairable
    partial    mixed multiplicities — a partial second write, left alone
    conflict   a record occurs twice with DIFFERENT content — left alone, inspect by hand
    empty      no records
    """
    parse = parse_gff if name.endswith(".gff3") else parse_fasta
    head, records = parse(raw)
    if not records:
        return "empty", None, 0, 0

    multiplicity = Counter(key for key, _ in records)
    distinct = len(multiplicity)
    counts = set(multiplicity.values())

    if counts == {1}:
        return "single", None, len(records), len(records)
    if counts != {2}:
        return "partial", None, len(records), distinct

    # Every record twice. Require the two occurrences to be byte-identical before trusting
    # that one of them can be dropped: a write spliced mid-record would show up here.
    bodies: dict[bytes, bytes] = {}
    for key, body in records:
        if key in bodies:
            if bodies[key] != body:
                return "conflict", None, len(records), distinct
        else:
            bodies[key] = body

    out = [*head]
    emitted: set[bytes] = set()
    for key, body in records:
        if key in emitted:
            continue
        emitted.add(key)
        out.append(key)
        if body:
            out.append(body)
    return "doubled", b"".join(out), len(records), distinct


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true", help="report only, write nothing")
    g.add_argument("--apply", action="store_true", help="rewrite verified doublings")
    args = ap.parse_args()

    if not PREDICTED.is_dir():
        print(f"no prediction tree at {PREDICTED}", file=sys.stderr)
        return 1

    rows: list[tuple] = []
    counts: Counter[str] = Counter()
    touched: set[str] = set()
    flagged: list[tuple] = []

    accessions = sorted(p for p in PREDICTED.iterdir() if p.is_dir())
    for d in accessions:
        for name in TARGETS:
            f = d / name
            if not f.is_file() or f.stat().st_size == 0:
                counts["missing"] += 1
                rows.append((d.name, name, "missing", 0, 0, 0, 0))
                continue
            raw = f.read_bytes()
            status, repaired, n_before, n_after = classify(raw, name)
            counts[status] += 1
            before, after = len(raw), len(raw)
            if status == "doubled":
                touched.add(d.name)
                after = len(repaired)
                if args.apply:
                    tmp = f.with_name(f.name + ".tmp")
                    tmp.write_bytes(repaired)
                    tmp.replace(f)
            elif status in ("partial", "conflict"):
                flagged.append((d.name, name, status, n_before, n_after))
            rows.append((d.name, name, status, before, after, n_before, n_after))

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT.open("w") as fh:
        fh.write("accession\tfile\tstatus\tbytes_before\tbytes_after\t"
                 "records_before\trecords_after\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    verb = "repaired" if args.apply else "would repair"
    print(f"{len(rows)} files over {len(accessions)} accessions")
    for k in ("single", "doubled", "partial", "conflict", "empty", "missing"):
        if counts[k]:
            print(f"  {k:9s} {counts[k]}")
    print(f"{verb} {counts['doubled']} files across {len(touched)} accessions")
    if flagged:
        print("  LEFT ALONE — inspect by hand:")
        for a, n, s, nb, na in flagged:
            print(f"    {a}  {n}  {s}  records {nb} -> {na} distinct")
    print(f"report -> {REPORT}")
    if touched:
        ACCESSIONS.write_text("\n".join(sorted(touched)) + "\n")
        print(f"affected accessions -> {ACCESSIONS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
