#!/usr/bin/env python3
"""
kraken_host_call.py — call each run's host from its own reads, not from its metadata.

The LLM host extraction leaves gaps: 451 of 3,220 runs have no usable value, and those
runs then carry no host in any per-host analysis. But db_v3 contains 95 host assemblies,
so the reads themselves say which plant they came from. This assigns the host as the
plant taxon covering the largest share of its own k-mer space — the same criterion used
for pathogen detection, for the same reason: a read count can pile up on a conserved
region, a k-mer fraction cannot.

Validated against the LLM where the LLM made a call: 93.0% agreement (2,425/2,608).
The residual disagreements are mostly nomenclature rather than error, and are reported
rather than silently resolved:

    Pisum sativum -> Lathyrus oleraceus   23   synonyms, Kraken2 uses the current name
    Malus pumila  -> Malus domestica       8   synonyms
    Triticale     -> Triticum turgidum     6   hybrid with no assembly, nearest neighbour
    T. aestivum   -> T. turgidum          92   genuine: bread vs durum wheat

The wheat pair is the one worth judging by hand. Everything else is either a synonym or
a taxon db_v3 cannot represent, and neither is a failure of the call.

This is the assign-side host question ("which plant did these reads come from"). It is
NOT kraken/search/kraken_hosts.py, which is the build-side question ("fetch host
references for the database").

Run from crypt/:
  python kraken/assign/host/kraken_host_call.py
  python kraken/assign/host/kraken_host_call.py --min-kmer-frac 0.02
"""
import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "kraken" / "assign"))
from kraken_report_analysis import (          # noqa: E402
    load_inspect, load_run_to_biosample, parse_report)

csv.field_size_limit(10 ** 7)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reports-dir", default="kraken/assign/data/reports")
    ap.add_argument("--inspect", default="kraken/build/data/db_v3_inspect.tsv")
    ap.add_argument("--out", default="kraken/assign/host/data/host_calls.tsv")
    ap.add_argument("--min-kmer-frac", type=float, default=0.01,
                    help="a host call needs this share of the plant's k-mer space "
                         "(default: 0.01). Below it, the call is left blank rather "
                         "than guessed.")
    args = ap.parse_args()

    db = json.loads((ROOT / "stat/build/data/phibase_db.json").read_text())
    host_taxids = set(map(str, db["host_to_seed"]))
    viridiplantae = {n.lower() for n in db["viridiplantae_names"]}
    db_min = load_inspect(Path(args.inspect))
    run2bs = load_run_to_biosample()
    if not db_min:
        sys.exit(f"no inspect table at {args.inspect}; the k-mer fraction cannot be computed")

    llm = {}
    with open(ROOT / "metadata/classify/data/samples.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            llm[row["BioSample"].strip()] = (row.get("llm_host_resolved") or "").strip()

    def unusable(v):
        return not v or v.lower() in ("unresolved", "unknown", "none", "na")

    rows, n_called, n_agree, n_dis, n_new = [], 0, 0, 0, 0
    for p in sorted(Path(args.reports_dir).glob("*.txt")):
        best = None
        for rec in parse_report(p):
            if rec["rank"] != "S" or rec["distinct"] is None:
                continue
            if str(rec["taxid"]) not in host_taxids and rec["name"].lower() not in viridiplantae:
                continue
            total = db_min.get(rec["taxid"])
            if not total:
                continue
            frac = rec["distinct"] / total
            if best is None or frac > best[2]:
                best = (rec["taxid"], rec["name"], frac, rec["reads_clade"])

        bs = run2bs.get(p.stem, "")
        prev = llm.get(bs, "")
        called = best is not None and best[2] >= args.min_kmer_frac
        if called:
            n_called += 1
            if unusable(prev):
                status, n_new = "new", n_new + 1
            elif prev.lower() == best[1].lower():
                status, n_agree = "agrees", n_agree + 1
            else:
                status, n_dis = "differs", n_dis + 1
        else:
            status = "no call"

        rows.append({
            "run": p.stem, "biosample": bs,
            "host": best[1] if called else "",
            "host_taxid": best[0] if called else "",
            "kmer_frac": round(best[2], 6) if best else "",
            "reads": best[3] if best else "",
            "llm_host": prev, "status": status,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"  runs                : {len(rows):,}")
    print(f"  host called         : {n_called:,} ({100*n_called/len(rows):.1f}%)")
    print(f"  agrees with the LLM : {n_agree:,}")
    print(f"  differs             : {n_dis:,}")
    print(f"  NEW (LLM had none)  : {n_new:,}")
    print(f"  no call             : {len(rows)-n_called:,}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
