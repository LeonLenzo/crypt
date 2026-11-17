#!/bin/bash
#SBATCH --job-name=predict_congener
#SBATCH --account=pawsey1168
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/congener_%j.out
#SBATCH --error=/scratch/pawsey1168/llenzo/crypt/research/predict/logs/congener_%j.err

# PRECISION CHECK, the gap BUSCO cannot touch. BUSCO scores completeness of the CONSERVED
# CORE; it says nothing about whether predicted sequence is from the right ORGANISM. For a
# competitive k-mer index that is the load-bearing property: a false gene corrupts every
# species apportionment. These targets have no annotation of their own, so the nearest
# available truth is an ANNOTATED CONGENER.
#
# Containment of the prediction in the congener CDS, FracMinHash k=31 scaled=200, the SAME
# implementation and parameters as compare.py / the saturation tables, so the numbers are
# directly comparable to the n=60 benchmark (precision 0.945 against a species own CDS).
#
# READ IT AS A FLOOR, NOT A SCORE. Congener containment is bounded ABOVE by real between-
# species divergence, so a genuine prediction scores BELOW its benchmark self-containment,
# by an amount that grows with how far the congener sits. A HIGH value means predicted
# k-mers land in the right clade; a LOW value is ambiguous (could be a bad prediction OR a
# distant congener) and must be read against the pair divergence, not alone. The honest
# signal is the CONTRAST across targets at similar congener distance, and the floor against
# the null below.
#
# NULL: each prediction is ALSO compared to an OUT-OF-GENUS CDS (the next target s congener),
# so cross-genus containment gives the floor a real prediction must clear.

set -uo pipefail
module load cray-python/3.12.12   # NOT python/3.11.6, which has no numpy; the saturation
                                 # and compare jobs use cray-python for exactly this reason
R=/scratch/pawsey1168/llenzo/crypt
P=$R/research/predict
CDS=$R/research/02_classify/data/bronze/cds_db0/pathogen   # cds_db0, NOT cds: the plain cds/
                                                         # tree holds empty dirs; the real CDS are here
PRED=$P/data/predicted
cd "$P/../02_classify" || exit 1     # so compare.py-style imports resolve _layout

python3 - "$P/data/congener_map.tsv" "$PRED" "$CDS" "$P/data/congener_check.tsv" <<PY
import sys, csv, glob, importlib.util
from pathlib import Path
import numpy as np

mapf, pred_dir, cds_dir, outf = sys.argv[1:5]
pred_dir, cds_dir = Path(pred_dir), Path(cds_dir)

# reuse the one validated sketcher
sp = importlib.util.spec_from_file_location(
    "_sat", "utilities/saturation/saturation.py")
sat = importlib.util.module_from_spec(sp); sp.loader.exec_module(sat)
sketch = sat.sketch_fasta
K, SCALED = 31, 200

def truth_cds(acc):
    d = cds_dir / acc
    if not d.is_dir(): return None
    c = [f for f in d.glob("**/*.fna")
         if f.is_file() and not f.name.endswith(("_clean.fna","_combined.fna",".tagged.fna"))
         and f.stat().st_size >= 100_000]
    return max(c, key=lambda f: f.stat().st_size) if c else None

def contain(a, b):
    if a.size == 0: return float("nan")
    return float(np.intersect1d(a, b, assume_unique=True).size) / float(a.size)

rows = list(csv.DictReader(open(mapf), delimiter="\t"))
# pre-sketch predictions and congener CDS
pS, cS = {}, {}
for r in rows:
    coding = pred_dir / r["pred_acc"] / "coding.fa"
    if coding.exists(): pS[r["pred_acc"]] = sketch(str(coding), k=K, scaled=SCALED)
    t = truth_cds(r["congener_acc"])
    if t is not None: cS[r["congener_acc"]] = sketch(str(t), k=K, scaled=SCALED)

out = []
for i, r in enumerate(rows):
    pa, ca = r["pred_acc"], r["congener_acc"]
    if pa not in pS or ca not in cS:
        out.append({**r, "containment": "NA", "null_containment": "NA",
                    "note": "missing sketch"}); continue
    cong = contain(pS[pa], cS[ca])
    # null: this prediction vs a DIFFERENT target s congener (out of genus)
    others = [rows[j]["congener_acc"] for j in range(len(rows))
              if rows[j]["congener_acc"] in cS and rows[j]["congener_acc"] != ca]
    null = max((contain(pS[pa], cS[o]) for o in others), default=float("nan"))
    out.append({**r, "containment": f"{cong:.4f}", "null_containment": f"{null:.4f}",
                "note": ""})

cols = ["pred_acc","taxon","congener_acc","congener","containment","null_containment","note"]
with open(outf, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(out)
hdr = ("taxon","congener","contain","null")
print("%-28s%-30s%9s%8s" % hdr)
for r in out:
    print("%-28s%-30s%9s%8s" % (r["taxon"][:27], r["congener"][:29], r["containment"], r["null_containment"]))
PY
echo "-> $P/data/congener_check.tsv"
