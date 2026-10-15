import csv, sys, math, statistics, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _common import detections          # the shared detection predicate; see 03_kraken/_common.py
from _paths import FILTER_DATA
pts=collections.defaultdict(list)
for run, h, rec in detections(require_host=True):
    pts[(rec['name'],h)].append((rec['reads_clade'],rec['distinct']))
def ols(ps):
    xs=[math.log10(a) for a,_ in ps]; ys=[math.log10(b) for _,b in ps]
    mx=statistics.mean(xs); my=statistics.mean(ys); sxx=sum((x-mx)**2 for x in xs)
    if sxx==0: return None,None
    b=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sxx; return my-b*mx,b
def classify(ps,band=0.4):
    keep=list(ps)
    for _ in range(6):
        a,b=ols(keep)
        nk=[(r,d) for r,d in ps if math.log10(d)-(a+b*math.log10(r))>-band]
        if len(nk)==len(keep) or len(nk)<8: break
        keep=nk
    a,b=ols(keep); rising=b>=0.5
    out=[]
    for r,d in ps:
        real = rising and (math.log10(d)-(a+b*math.log10(r))>-band)
        out.append((r,d,'real' if real else 'crash-out'))
    return out,a,b,rising
show=[('Bipolaris sorokiniana','Triticum aestivum'),('Puccinia graminis','Triticum aestivum'),
      ('Puccinia striiformis','Triticum aestivum'),('Phakopsora pachyrhizi','Triticum aestivum')]
with open(str(FILTER_DATA / 'splitfig.tsv'),'w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(['panel','reads','distinct','cls','a','b'])
    for sp,h in show:
        ps=pts[(sp,h)]; out,a,b,rising=classify(ps)
        lab=f'{sp} on {h.split()[0]}'
        for r,d,c in out: w.writerow([lab,r,d,c,round(a,4),round(b,4)])
print('wrote 03_kraken/05_filter/data/splitfig.tsv')
