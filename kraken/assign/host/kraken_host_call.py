#!/usr/bin/env python3
"""
kraken_host_call.py — call each run's host from its own reads, not from its metadata.

The LLM host extraction leaves gaps: 451 of 3,220 runs have no usable value, and those
runs then carry no host in any per-host analysis. But db_v3 contains 95 host assemblies,
so the reads themselves say which plant they came from. This assigns the host as the plant taxon with the most reads, and it deliberately does
NOT apply the k-mer-fraction criterion used for pathogen detection. Detection asks
whether an organism is present at all and needs an absolute bar; host assignment is a
ranking question, because every plant RNA-seq library has a host, so withholding a call
below a threshold just manufactures gaps. Measured both ways over 3,220 runs: ranking
by reads agrees with the LLM 91.2% of the time against 90.2% for k-mer fraction, and
either way every run gets a call. Reads is also the sounder key here, since the k-mer
fraction normalises by how much of that plant sits in db_v3 and so favours hosts with
smaller CDS sets, which has nothing to do with which plant the library came from.

The k-mer fraction is still reported per call, as a confidence signal rather than a
gate, and it earns its place: it is the only thing separating "host identified" from
"host absent from db_v3, here is its nearest relative". Fragaria is not in db_v3 at
all, so the 30 runs the LLM called strawberry are called Populus x canadensis (17),
Arachis hypogaea (7) and Glycine max (6). Ranking alone cannot detect that; low
coverage can. A call resting on under 1% of the host's k-mer space is flagged
low_coverage, and those flags are diagnostic of database gaps rather than weak data.

A second caveat, and the reason a correct host can still score low: polyploids lose
k-mer space to the LCA when a relative sharing a subgenome is also in the database.
db_v3 holds hexaploid T. aestivum (ABD) and tetraploid T. turgidum (AB), so 14.4M of
the 51.8M minimizers in the Triticum clade (27.7%) sit at the genus node rather than
on either species. That also drives the 92 T. aestivum -> T. turgidum calls. For
polyploid hosts, read the genus-level signal alongside the species call.

Validated against the LLM where the LLM made a call: 91.2% agreement.
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
    ap.add_argument("--rank-by", choices=("reads", "kmer"), default="reads",
                    help="which signal ranks the candidate hosts (default: reads)")
    ap.add_argument("--low-coverage-below", type=float, default=0.01,
                    help="flag (do not drop) calls resting on less than this share of "
                         "the host's k-mer space (default: 0.01)")
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
            key = rec["reads_clade"] if args.rank_by == "reads" else frac
            if best is None or key > best[4]:
                best = (rec["taxid"], rec["name"], frac, rec["reads_clade"], key)

        bs = run2bs.get(p.stem, "")
        prev = llm.get(bs, "")
        called = best is not None
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
            "confidence": ("low_coverage" if called and best[2] < args.low_coverage_below
                           else "ok" if called else ""),
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
    low = sum(1 for r in rows if r["confidence"] == "low_coverage")
    print(f"  no call             : {len(rows)-n_called:,}")
    print(f"  flagged low_coverage: {low:,} (called, but on <1% of the host's k-mer space)")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
