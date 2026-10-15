#!/usr/bin/env python3
"""kraken_annotate.py — predict CDS for an unannotated assembly by protein homology.

Why this exists: 73% of assemblies in our pathogen genera carry no NCBI annotation, so they
have no `cds_from_genomic.fna` and are invisible to the CDS-based pipeline. The 04_kallisto
blind spot is the concrete cost: Puccinia hordei (barley leaf rust, a distinct species) has a
chromosome-level assembly but no CDS, so barley leaf rust is missed or misattributed to
P. triticina. See memory crypt-annotation-gap and crypt-align-crosstalk-measurement.

Method: homology-guided, not de novo. The annotated congeners' proteins are projected onto the
unannotated target genome with miniprot, and the aligned gene models' CDS are extracted. This
is deliberate: in a competitive index an invented gene corrupts every species' apportionment,
so precision matters more than recall, and homology projection is high-precision where ab
initio is not. The cost is a known bias: predicted CDS is enriched for genes conserved with the
guide, so a PREDICTED reference is valid ONLY as a competitor (does a present species' signal
survive when this one is in the index), NEVER for positive calls about the predicted species.
A marker file records this so downstream code can tell predicted CDS from downloaded.

Guides are the annotated congeners whose CDS are already on disk; their proteins come from
translating that CDS, so no protein download is needed.

Output mirrors the download layout exactly, so 02_build treats it like any other accession:
    cds/pathogen/{target}/ncbi_dataset/data/{target}/cds_from_genomic.fna
    cds/pathogen/{target}/PREDICTED                      (marker: guide accessions, date, tool)

Run on Setonix (builds miniprot if absent; needs the datasets CLI, both wired in the SLURM job):
    sbatch 03_kraken/slurm/kraken_annotate.slurm
Direct (from crypt/, inside a job — downloads a genome and runs miniprot, never on a login node):
    python 03_kraken/utilities/annotate/kraken_annotate.py \
        --target GCA_062009444.1 \
        --guides GCA_977970215.1,GCA_025169535.1,GCA_008522505.1 \
        --threads 16
"""
import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CDS_DIR = ROOT / "03_kraken/01_search/data/cds/pathogen"
WORK = Path(__file__).resolve().parent / "data"

# Standard genetic code, for translating guide CDS to the protein evidence miniprot needs.
_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(_CODONS, _AA))


def sh(cmd, **kw):
    print(f"  $ {' '.join(str(c) for c in cmd)}", flush=True)
    r = subprocess.run(cmd, **kw)
    if r.returncode != 0:
        raise SystemExit(f"command failed (rc={r.returncode}): {' '.join(map(str, cmd))}")
    return r


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def translate(nt):
    """Translate an in-frame CDS to protein, stopping at the first stop codon."""
    aa = []
    for i in range(0, len(nt) - 2, 3):
        a = CODON_TABLE.get(nt[i:i + 3].upper(), "X")
        if a == "*":
            break
        aa.append(a)
    return "".join(aa)


def cds_fasta_for(accession):
    return CDS_DIR / accession / "ncbi_dataset/data" / accession / "cds_from_genomic.fna"


def make_guide_proteome(guides, out_faa, min_aa=40):
    """Translate the guides' CDS into one deduplicated protein FASTA.

    Deduplication matters: several congeners share near-identical housekeeping proteins, and
    feeding miniprot the same protein many times wastes alignment time without adding evidence.
    """
    seen, n = set(), 0
    with open(out_faa, "w") as fh:
        for acc in guides:
            cds = cds_fasta_for(acc)
            if not cds.exists():
                print(f"  !! guide {acc} has no CDS at {cds}", file=sys.stderr)
                continue
            for hdr, nt in read_fasta(cds):
                prot = translate(nt)
                if len(prot) < min_aa or prot in seen:
                    continue
                seen.add(prot)
                fh.write(f">{acc}|{hdr.split()[0]}\n{prot}\n")
                n += 1
    print(f"  guide proteome: {n} distinct proteins from {len(guides)} assemblies")
    return n


def ensure_miniprot(bindir):
    """Build miniprot into bindir if it is not already there. ~2 min, no deps beyond zlib."""
    exe = bindir / "miniprot"
    if exe.exists():
        return exe
    bindir.mkdir(parents=True, exist_ok=True)
    src = WORK / "miniprot_src"
    if not src.exists():
        sh(["git", "clone", "--depth", "1", "https://github.com/lh3/miniprot", str(src)])
    sh(["make", "-C", str(src), "-j4"])
    # copy, not rename: /scratch and /software are different filesystems, so os.replace fails
    # with EXDEV (Invalid cross-device link)
    shutil.copy2(src / "miniprot", exe)
    exe.chmod(0o755)
    print(f"  built miniprot -> {exe}")
    return exe


def download_genome(accession, dest):
    """Fetch an assembly's genomic FASTA via the datasets CLI, matching kraken_search's pattern."""
    gfna = dest / f"{accession}.genome.fna"
    if gfna.exists() and gfna.stat().st_size:
        print(f"  genome present: {gfna.name}")
        return gfna
    dest.mkdir(parents=True, exist_ok=True)
    zp = dest / "genome.zip"
    sh(["datasets", "download", "genome", "accession", accession,
        "--include", "genome", "--filename", str(zp)])
    with zipfile.ZipFile(zp) as z:
        fna = [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta"))
               and "cds" not in n.lower()]
        if not fna:
            raise SystemExit(f"no genomic fasta in {zp}")
        with z.open(fna[0]) as s, open(gfna, "wb") as o:
            o.write(s.read())
    zp.unlink()
    print(f"  downloaded genome: {gfna.name} ({gfna.stat().st_size/1e6:.0f} MB)")
    return gfna


def run_miniprot(exe, genome, proteome, out_gff, threads):
    if out_gff.exists() and out_gff.stat().st_size:
        print(f"  miniprot output present, reusing {out_gff.name}")
        return out_gff
    tmp = out_gff.with_suffix(".gff.tmp")        # write to tmp so a killed run does not look done
    with open(tmp, "w") as fh:
        sh([str(exe), "-t", str(threads), "--gff", "-u", str(genome), str(proteome)], stdout=fh)
    tmp.replace(out_gff)
    return out_gff


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def extract_cds(gff, genome, out_cds, target):
    """Pull each mRNA's CDS from the genome, splicing its CDS features in order.

    miniprot writes one mRNA per alignment with ordered CDS children; concatenating the CDS
    substrings (reverse-complemented on the minus strand) reconstructs the coding sequence.
    Headers mirror the NCBI cds_from_genomic convention closely enough for the 02_build tagger.
    """
    # GFF references the bare contig id (first token), so key on that, not the full header
    seqs = {hdr.split()[0]: seq for hdr, seq in read_fasta(genome)}
    # bucket CDS lines by their parent mRNA, preserving strand and contig
    mrna, cds = {}, {}
    for line in open(gff):
        if line.startswith("#") or "\t" not in line:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ctg, _, feat, start, end, _, strand, _, attr = f[:9]
        at = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)
        if feat == "mRNA":
            mrna[at.get("ID", "")] = (ctg, strand, at)
            cds.setdefault(at.get("ID", ""), [])
        elif feat == "CDS":
            cds.setdefault(at.get("Parent", ""), []).append((int(start), int(end)))

    # build one record per model, then collapse overlapping models to one gene per locus.
    # Each guide protein projects independently, so three congeneric guides give ~3 models at
    # the same P. hordei locus; left in, a real read there is split across the duplicates by the
    # EM and the present species looks like a cross-map. Keep the highest-identity model per
    # cluster of overlapping same-strand CDS.
    models = []
    for mid, (ctg, strand, at) in mrna.items():
        parts = sorted(cds.get(mid, []))
        if not parts or ctg not in seqs:
            continue
        seq = "".join(seqs[ctg][s - 1:e] for s, e in parts)
        if strand == "-":
            seq = revcomp(seq)
        if len(seq) < 150:                          # drop fragments below ~50 aa
            continue
        try:
            ident = float(at.get("Identity", "0"))
        except ValueError:
            ident = 0.0
        models.append({"ctg": ctg, "strand": strand, "beg": parts[0][0], "end": parts[-1][1],
                       "ident": ident, "mid": mid, "seq": seq})

    models.sort(key=lambda m: (m["ctg"], m["strand"], m["beg"], m["end"]))
    kept, cluster_end, best = [], None, None
    for m in models:
        if (best is not None and m["ctg"] == best["ctg"] and m["strand"] == best["strand"]
                and m["beg"] <= cluster_end):       # overlaps the open cluster
            cluster_end = max(cluster_end, m["end"])
            if m["ident"] > best["ident"]:
                best = m
        else:
            if best is not None:
                kept.append(best)
            best, cluster_end = m, m["end"]
    if best is not None:
        kept.append(best)

    with open(out_cds, "w") as fh:
        for m in kept:
            fh.write(f">lcl|{target}_pred_{m['mid']} [gene={m['mid']}] "
                     f"[identity={m['ident']:.4f}] "
                     f"[location={m['ctg']}:{m['beg']}-{m['end']}]\n")
            for i in range(0, len(m["seq"]), 80):
                fh.write(m["seq"][i:i + 80] + "\n")
    print(f"  {len(mrna)} models -> {len(models)} valid -> {len(kept)} after locus dedup "
          f"-> {out_cds.name}")
    return len(kept)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True, help="unannotated assembly accession to annotate")
    ap.add_argument("--guides", required=True,
                    help="comma-separated annotated congener accessions (CDS already on disk)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bindir", type=Path,
                    default=Path("/software/projects/pawsey1168/llenzo/bin"),
                    help="where to build/find miniprot")
    args = ap.parse_args()

    guides = [g.strip() for g in args.guides.split(",") if g.strip()]
    stage = WORK / args.target
    stage.mkdir(parents=True, exist_ok=True)
    print(f"annotating {args.target} from {len(guides)} guides: {', '.join(guides)}")

    exe = ensure_miniprot(args.bindir)
    proteome = stage / "guides.faa"
    if make_guide_proteome(guides, proteome) == 0:
        raise SystemExit("no guide proteins; check the guide accessions have CDS on disk")
    genome = download_genome(args.target, stage)
    gff = run_miniprot(exe, genome, proteome, stage / "miniprot.gff", args.threads)

    # write CDS into the download tree so 02_build cannot tell it apart, plus a provenance marker
    dest = CDS_DIR / args.target / "ncbi_dataset/data" / args.target
    dest.mkdir(parents=True, exist_ok=True)
    out_cds = dest / "cds_from_genomic.fna"
    written = extract_cds(gff, genome, out_cds, args.target)
    if written == 0:
        raise SystemExit("miniprot produced no usable gene models")

    (CDS_DIR / args.target / "PREDICTED").write_text(
        f"predicted by kraken_annotate.py (miniprot homology projection)\n"
        f"guides: {','.join(guides)}\n"
        f"n_cds: {written}\n"
        f"USE AS A COMPETITOR ONLY, not for positive calls about this species "
        f"(homology bias; see docstring).\n")
    print(f"done: {written} predicted CDS for {args.target}, marked PREDICTED")


if __name__ == "__main__":
    main()
