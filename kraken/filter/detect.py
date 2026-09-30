#!/usr/bin/env python3
"""Per-detection real/crash-out via a weighted SCORE of five features (per host x species).

Each feature is 0-1; the score is their weighted mean; real if score >= 0.5. A weighted
score, not hard gates, so a shallow slope no longer flips a whole cloud when the fit is
tight (fixes nicotianae/pseudograminearum), and each feature's contribution is inspectable.

  tightness  R2 of the upper accumulation line   (is there a clean real line at all)   w .25
  rising     slope of that line, clipped          (coverage grows with depth)           w .15
  coverage   distinct / genome ceiling, log       (breadth of genome seen; saturation)  w .25
  on_line    this point's residual vs the line    (sits on the real line, not the floor) w .15
  unbiased   1 - f.sp. fraction                    (not biased to a sub-taxon sliver)    w .20
"""
import csv, sys, math, collections
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "kraken" / "assign"))
from kraken_report_analysis import load_reference, parse_report, load_inspect  # noqa
csv.field_size_limit(10 ** 7)
pathogens, hosts = load_reference()
CEIL = load_inspect(ROOT / "kraken/build/data/db_v3_inspect.tsv")
SUBS = collections.defaultdict(list); cur=None
for line in open(ROOT / "kraken/build/data/db_v3_inspect.tsv"):
    f=line.rstrip("\n").split("\t")
    if len(f)<6: continue
    if f[3].strip()=="S": cur=int(f[4]); SUBS.setdefault(cur,[])
    elif f[3].strip() in ("S1","S2") and cur is not None: SUBS[cur].append(int(f[4]))
host_of={}; bs_of={}
for r in csv.DictReader(open(ROOT/"kraken/assign/host/data/host_calls.tsv"),delimiter="\t"):
    host_of[r["run"].strip()]=(r.get("host") or "").strip()
    bs_of[r["run"].strip()]=(r.get("biosample") or "").strip()

# each point carries its run so the per-run/biosample call survives to the output
pts=collections.defaultdict(lambda: collections.defaultdict(list))
for p in sorted((ROOT/"kraken/assign/data/reports").glob("*.txt")):
    h=host_of.get(p.stem,"") or "unresolved"
    recs=list(parse_report(p)); sub_reads={rec["taxid"]:rec["reads_clade"] for rec in recs}
    for rec in recs:
        if rec["rank"]!="S" or rec["distinct"] is None: continue
        tid=int(rec["taxid"])
        if str(tid) in hosts or str(tid) not in pathogens or rec["reads_clade"]<100 or rec["distinct"]<1: continue
        subr=sum(sub_reads.get(s,0) for s in SUBS.get(tid,[]))
        pts[(rec["name"],tid)][h].append((rec["reads_clade"],rec["distinct"],subr/rec["reads_clade"],p.stem))

def fit(P):
    xs=[math.log10(a) for a,_ in P]; ys=[math.log10(b) for _,b in P]; n=len(P)
    mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); syy=sum((y-my)**2 for y in ys); sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    if sxx==0: return None
    b=sxy/sxx; return my-b*mx, b, (sxy**2/(sxx*syy) if syy>0 else 1.0)

# Balanced feature weights (weighted mean, real if >=0.5). NOT the fitted logistic:
# fitting on PHI-base labels made coverage dominate (declared pathogens are abundant),
# which crushed low-abundance real detections. Dropping coverage over-rejected obvious
# reals (striiformis 0%). These hand weights sit between the two: tightness + rising +
# moderate coverage. This is a provisional operating point pending alignment ground truth.
W=dict(tight=.20, rising=.25, cov=.20, on_line=.15, unbiased=.20)
def clip(v,lo=0.0,hi=1.0): return max(lo,min(hi,v))

def score_cloud(P, ceil, band=0.4):
    rd=[(e[0],e[1]) for e in P]; keep=list(rd)
    for _ in range(6):
        f=fit(keep)
        if not f: break
        a,b,_=f; nk=[(r,d) for r,d in rd if math.log10(d)-(a+b*math.log10(r))>-band]
        if len(nk)==len(keep) or len(nk)<8: break
        keep=nk
    f=fit(keep)
    if not f: return [(0.0,"crash-out")]*len(P)
    a,b,r2=f
    F_tight=clip(r2); F_rising=clip(b/0.6)
    has_sub=max(e[2] for e in P)>0
    out=[]
    for (r,d,fr,_run) in P:
        cov = clip((math.log10(d/ceil)+3)/3) if ceil else 0.3   # frac 0.001->0 .. 1->1
        resid = math.log10(d)-(a+b*math.log10(r))
        F_on = clip(1+resid/band)                               # 1 on/above line, 0 a band below
        F_unb = (1-clip(fr)) if has_sub else 1.0
        s = (W["tight"]*F_tight + W["rising"]*F_rising + W["cov"]*cov +
             W["on_line"]*F_on + W["unbiased"]*F_unb)
        out.append((round(s,3), "real" if s>=0.5 else "crash-out"))
    return out

keep_sp=set(); rows=[]
for (sp,tid),hs in pts.items():
    ceil=CEIL.get(tid)
    for h,P in hs.items():
        span=(max(math.log10(e[0]) for e in P)-min(math.log10(e[0]) for e in P)) if P else 0
        if len(P)>=15 and span>=0.7:
            res=score_cloud(P,ceil); keep_sp.add(sp)
        else:
            res=[(float("nan"),"unclassified")]*len(P)
        for (r,d,fr,run),(s,c) in zip(P,res):
            cov = round(d/ceil, 5) if ceil else ""      # fraction of genome k-mer space seen in this run
            rows.append([sp,h,run,bs_of.get(run,""),r,d,round(fr,3),cov,s,c])
out=ROOT/"kraken/filter/data/_classified.tsv"
with open(out,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["species","host","run","biosample","reads","distinct","fsp_frac","coverage","score","cls"]); w.writerows(rows)
print(f"  {len(keep_sp)} species scored, {len(rows)} detections")
