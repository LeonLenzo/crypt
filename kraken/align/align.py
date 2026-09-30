#!/usr/bin/env python3
"""align: competitive read assignment to confirm co-infections (SKETCH).

Pipeline (see README for the rationale):

  1. candidates(host)   flagged pathogens per host from filter/data/_classified.tsv
  2. build_index(host)  kallisto index of {candidates + congeners + host} CDS, tagged by species
  3. quant(run, host)   kallisto quant of a sample's reads against its host index
  4. summarise(run)     per species: uniquely-assigned reads, gene breadth, competitive survival

Runs on Setonix (kallisto via spack; reads on pawsey1168 scratch; CDS in
kraken/search/data/cds/pathogen/). This file is a skeleton: control flow and I/O contracts
are fixed, the cluster specifics (module load, SLURM array like kraken/assign) are TODO.
"""
import argparse, csv, collections, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CDS_DIR   = ROOT / "kraken/search/data/cds/pathogen"          # per-accession CDS FASTA (Setonix)
CLASSIFIED = ROOT / "kraken/filter/data/_classified.tsv"
DATA = Path(__file__).resolve().parent / "data"
SCORE_MIN, COV_MIN = 0.7, 0.01                                # match the filter's presence call


def candidates(host):
    """Species flagged present in `host` (score/coverage floors), from the filter output."""
    out = set()
    for r in csv.DictReader(open(CLASSIFIED), delimiter="\t"):
        if r["host"] != host:
            continue
        try:
            if float(r["score"]) >= SCORE_MIN and float(r["coverage"]) >= COV_MIN:
                out.add(r["species"])
        except ValueError:
            pass
    return sorted(out)


def reference_set(species):
    """Species to put in the index: the candidates + their congeners (same genus).

    Congeners make within-genus reads compete, which is where cross-mapping hides. The host
    CDS is added separately by build_index(). Returns species -> [CDS fasta paths].
    TODO: map species -> accession dir(s) under CDS_DIR (via the db build manifest /
    assembly_data_report taxids), expand to all genera present in CDS_DIR.
    """
    raise NotImplementedError


def build_index(host):
    """Build one kallisto index per host from tagged CDS of {refs + host}.

    Tag every transcript header with its species so abundances aggregate up (mirror the
    .tagged.fna convention in kraken/build). Then: kallisto index -i data/idx/{host}.idx.
    TODO: concat tagged CDS; run kallisto index.
    """
    raise NotImplementedError


def quant(run, host, reads):
    """kallisto quant one sample against its host index -> data/quant/{host}/{run}/.

    Paired where mates exist, single otherwise (same read layout as kraken/assign).
    TODO: kallisto quant -i {host}.idx -o out [reads]; capture run_info.json.
    """
    raise NotImplementedError


def summarise(run, host):
    """Per species in the sample: uniquely-assigned reads, gene breadth, competitive survival.

    From abundance.tsv (est_counts, tpm) aggregated by the species tag, plus the equivalence
    classes for the unique-vs-shared split. Competitive survival = est_counts with congeners
    in the index vs the species quantified alone. Emits data/align_calls.tsv.
    TODO: parse abundance.tsv + run_info; aggregate by tag; join to the flagged co-infections.
    """
    raise NotImplementedError


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--host", help="host species to build index / quant for")
    ap.add_argument("--stratified", action="store_true",
                    help="first pass: quant only a stratified sample of flagged detections")
    args = ap.parse_args()
    print("align.py is a sketch; see README for the plan.", file=sys.stderr)
    if args.host:
        print(f"{args.host}: {len(candidates(args.host))} flagged candidate pathogens")


if __name__ == "__main__":
    main()
