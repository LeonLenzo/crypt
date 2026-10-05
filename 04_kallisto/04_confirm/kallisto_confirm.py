#!/usr/bin/env python3
"""kallisto_confirm.py — turn pseudoalignment into a real/artefact call per detection.

Three read-level facts per (run, species), from 03_quant's two passes:

  unique_reads   reads whose equivalence class contains transcripts of ONE species only. A real
                 infection keeps many; a pure cross-map keeps ~0 because every k-mer it matched
                 is shared with the species it is bleeding from.
  gene_breadth   distinct genes carrying reads. A real infection spreads across the gene set;
                 a cross-map concentrates on a few conserved ones.
  survival       est_counts with the species' measured neighbours in the index, over est_counts
                 with the species quantified alone. Not computable from one index, so it needs
                 a second solo quant; see survival() below. It cannot run at all for the 11
                 flagged species with no neighbour above 1% containment (all three Puccinia,
                 Zymoseptoria tritici and P. nodorum among them), which are judged on the first
                 two facts only.

Reading equivalence classes without bustools
--------------------------------------------
bustools is not installed on Setonix and is not in Pawsey's spack repo, so the BUS file is
parsed directly. The format is small, documented and stable:

  header   magic "BUS\\0" (4B), version u32, bc_len u32, umi_len u32, tlen u32, tlen bytes text
  records  32 bytes each: barcode u64, UMI u64, ec i32, count u32, flags u32, pad u32

`matrix.ec` maps each ec id to transcript indices; `transcripts.txt` gives the target names in
index order, which are the `{species}|{accession}|{transcript}` tags 02_build wrote. So species
membership of an equivalence class is a lookup, and an EC whose species set has one member
contributes uniquely-assigned reads.

    python 04_kallisto/04_confirm/kallisto_confirm.py --run SRR27326517 --host "Arabidopsis thaliana"
    python 04_kallisto/04_confirm/kallisto_confirm.py --all
"""
import argparse
import collections
import csv
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CALLS, NEIGHBOURS, QUANT_DIR, STRATA, num, slug, step_log

BUS_MAGIC = b"BUS\x00"
BUS_RECORD = struct.Struct("<QQiIII")   # barcode, umi, ec, count, flags, pad


def read_transcripts(path):
    """Target names in index order. Index position is what matrix.ec refers to."""
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def species_of(target):
    """The species tag 02_build wrote as the first "|" field."""
    return target.split("|", 1)[0]


def read_ec_species(ec_path, transcripts):
    """ec id -> frozenset of species tags it is compatible with.

    Storing species sets rather than transcript lists collapses the big ECs early: the wheat
    index has ECs of up to 326 transcripts, but most resolve to a handful of species.
    """
    out = {}
    with open(ec_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ec = int(parts[0])
            idxs = (int(i) for i in parts[1].split(",") if i)
            out[ec] = frozenset(species_of(transcripts[i]) for i in idxs
                                if 0 <= i < len(transcripts))
    return out


def read_bus_counts(bus_path):
    """ec id -> read count, by streaming the BUS records."""
    counts = collections.Counter()
    with open(bus_path, "rb") as fh:
        magic = fh.read(4)
        if magic != BUS_MAGIC:
            raise ValueError(f"{bus_path}: not a BUS file (magic {magic!r})")
        _version, _bc_len, _umi_len, tlen = struct.unpack("<IIII", fh.read(16))
        fh.read(tlen)                      # free-text header, not needed
        size = BUS_RECORD.size
        while True:
            chunk = fh.read(size * 65536)
            if not chunk:
                break
            # a truncated tail means the bus pass was cut off; treat as corrupt rather than
            # silently undercounting, same lesson as the 2026-09-14 gzip incident
            if len(chunk) % size:
                raise ValueError(f"{bus_path}: truncated, {len(chunk) % size} trailing bytes")
            for off in range(0, len(chunk), size):
                _bc, _umi, ec, count, _flags, _pad = BUS_RECORD.unpack_from(chunk, off)
                counts[ec] += count
    return counts


def unique_and_shared(quant_dir):
    """species -> (unique_reads, shared_reads) from the equivalence classes."""
    bus = quant_dir / "bus"
    transcripts = read_transcripts(bus / "transcripts.txt")
    ec_species = read_ec_species(bus / "matrix.ec", transcripts)
    counts = read_bus_counts(bus / "output.bus")

    unique = collections.Counter()
    shared = collections.Counter()
    for ec, n in counts.items():
        sps = ec_species.get(ec)
        if not sps:
            continue
        if len(sps) == 1:
            unique[next(iter(sps))] += n
        else:
            for s in sps:
                shared[s] += n
    return unique, shared


def abundance_by_species(quant_dir):
    """species -> (est_counts, tpm, gene_breadth).

    Gene breadth counts distinct targets with non-zero est_counts. These are CDS records, one
    per annotated gene product, so a target is a gene for this purpose.
    """
    est = collections.Counter()
    tpm = collections.Counter()
    genes = collections.Counter()
    with open(quant_dir / "abundance.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = float(r["est_counts"])
            if c <= 0:
                continue
            s = species_of(r["target_id"])
            est[s] += c
            tpm[s] += float(r["tpm"])
            genes[s] += 1
    return est, tpm, genes


def neighbour_count():
    """species -> how many measured >=1% containment neighbours it has.

    Zero means the competitive-survival arm cannot run for that species at all.
    """
    n = collections.Counter()
    if NEIGHBOURS.exists():
        for r in csv.DictReader(open(NEIGHBOURS), delimiter="\t"):
            n[r["species"]] += 1
    return n


def survival(host, run, species):
    """est_counts with neighbours present vs the species alone.

    Needs a second quant against a neighbour-free index, which 02_build does not build yet.
    Deliberately not approximated: assigning every shared-EC read to the species would give an
    upper bound that looks like a measurement, and the whole point of this module is to stop
    treating an upper bound as evidence.
    """
    raise NotImplementedError("needs a solo index per species; see 02_build")


def facts_for(host, run):
    """The per-species facts for one quanted run, or None if its output is absent."""
    qd = QUANT_DIR / slug(host) / run
    if not (qd / "abundance.tsv").exists():
        return None
    est, tpm, genes = abundance_by_species(qd)
    unique, shared = ({}, {})
    if (qd / "bus" / "output.bus").exists():
        unique, shared = unique_and_shared(qd)
    out = {}
    for s in set(est) | set(unique):
        u, sh = unique.get(s, 0), shared.get(s, 0)
        out[s] = {"est_counts": round(est.get(s, 0.0), 2),
                  "tpm": round(tpm.get(s, 0.0), 2),
                  "gene_breadth": genes.get(s, 0),
                  "unique_reads": u,
                  "shared_reads": sh,
                  "unique_frac": round(u / (u + sh), 4) if (u + sh) else 0.0}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host")
    ap.add_argument("--run")
    ap.add_argument("--all", action="store_true", help="every quanted run in the strata")
    ap.add_argument("--out", type=Path, default=CALLS)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(STRATA), delimiter="\t"))
    want = [r for r in rows
            if (not args.host or r["host"] == args.host)
            and (not args.run or r["run"] == args.run)]
    if not want:
        sys.exit("no detections match")

    nb = neighbour_count()
    log = step_log(__file__, "kallisto_confirm")
    try:
        cache, out_rows, skipped = {}, [], 0
        for r in want:
            key = (r["host"], r["run"])
            if key not in cache:
                cache[key] = facts_for(*key)
            f = cache[key]
            if f is None:
                skipped += 1
                continue
            d = f.get(slug(r["species"]), {})
            out_rows.append({
                "host": r["host"], "run": r["run"], "species": r["species"],
                "stratum": r["stratum"], "expectation": r["expectation"],
                "kraken_coverage": r["coverage"], "kraken_score": r["score"],
                "est_counts": d.get("est_counts", 0.0), "tpm": d.get("tpm", 0.0),
                "gene_breadth": d.get("gene_breadth", 0),
                "unique_reads": d.get("unique_reads", 0),
                "shared_reads": d.get("shared_reads", 0),
                "unique_frac": d.get("unique_frac", 0.0),
                "n_neighbours": nb.get(r["species"], 0),
                "survival": "",   # requires a solo index; see survival()
            })

        if not out_rows:
            sys.exit(f"nothing to report: {skipped} detections have no quant output yet")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        cols = list(out_rows[0])
        with open(args.out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in out_rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"{args.out}: {len(out_rows)} detections"
              + (f" ({skipped} not yet quanted)" if skipped else ""))

        print(f"\n  {'stratum':12s} {'n':>4s} {'median unique':>14s} {'median breadth':>15s}")
        by = collections.defaultdict(list)
        for r in out_rows:
            by[r["stratum"].split(",")[0]].append(r)
        for st, rs in sorted(by.items()):
            u = sorted(r["unique_reads"] for r in rs)
            g = sorted(r["gene_breadth"] for r in rs)
            print(f"  {st:12s} {len(rs):4d} {u[len(u)//2]:14,d} {g[len(g)//2]:15,d}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
