import csv, sys, math, collections
import numpy as np
sys.path.insert(0,'kraken/assign')
from kraken_report_analysis import load_reference, load_inspect, parse_report
from pathlib import Path
csv.field_size_limit(10**7)
pathogens,hosts=load_reference()
db=load_inspect(Path('kraken/build/data/db_v3_inspect.tsv'))
name2tid={}
hostcall={}
for r in csv.DictReader(open('kraken/assign/host/data/host_calls.tsv'),delimiter='\t'):
    hostcall[r['run'].strip()]=r['host'].strip()
pts=collections.defaultdict(list)
for p in Path('kraken/assign/data/reports').glob('*.txt'):
    h=hostcall.get(p.stem)
    if not h: continue
    for rec in parse_report(p):
        if rec['rank']=='S' and rec['distinct'] and str(rec['taxid']) in pathogens and str(rec['taxid']) not in hosts and rec['reads_clade']>=100:
            pts[(rec['name'],h)].append((rec['reads_clade'],rec['distinct'],rec['taxid']))

def ols(x,y):
    b=np.polyfit(x,y,1); return b[1],b[0]  # intercept, slope

# ---- Method A: iterative upper line + rising-slope gate ----
def method_A(x,y,band=0.4):
    keep=np.ones(len(x),bool)
    for _ in range(6):
        a,b=ols(x[keep],y[keep])
        nk=(y-(a+b*x))>-band
        if nk.sum()==keep.sum() or nk.sum()<8: break
        keep=nk
    a,b=ols(x[keep],y[keep])
    if b<0.5: return np.zeros(len(x),bool)   # pure floor
    return (y-(a+b*x))>-band

# ---- Method B: 1D GMM on residual-from-global-line, pick k by BIC (2 or 3) ----
def gmm1d(v,k,iters=100):
    v=np.asarray(v); n=len(v)
    qs=np.quantile(v,np.linspace(0.1,0.9,k)); mu=qs.copy()
    var=np.full(k,v.var()/k+1e-6); w=np.full(k,1/k)
    for _ in range(iters):
        R=np.array([w[j]*np.exp(-(v-mu[j])**2/(2*var[j]))/np.sqrt(2*np.pi*var[j]) for j in range(k)])
        R/= R.sum(0)+1e-12
        Nk=R.sum(1)+1e-9
        mu=(R*v).sum(1)/Nk; var=(R*(v-mu[:,None])**2).sum(1)/Nk+1e-6; w=Nk/n
    ll=np.log(np.array([w[j]*np.exp(-(v-mu[j])**2/(2*var[j]))/np.sqrt(2*np.pi*var[j]) for j in range(k)]).sum(0)+1e-12).sum()
    bic=-2*ll+(3*k-1)*math.log(n)
    lab=np.argmax(R,0)
    return lab,mu,bic
def method_B(x,y):
    a,b=ols(x,y); resid=y-(a+b*x)
    best=None
    for k in (2,3):
        lab,mu,bic=gmm1d(resid,k)
        if best is None or bic<best[0]: best=(bic,lab,mu)
    _,lab,mu=best
    topcomp=np.argmax(mu)                      # highest-residual component = real
    return lab==topcomp

# ---- Method C: genome-ceiling rarefaction; fit distinct=G(1-exp(-cN/G)) to upper cloud ----
def method_C(x,y,tid):
    G=db.get(tid)
    if not G: return method_A(x,y)             # fall back
    N=10**x; D=10**y
    # estimate c from upper-cloud slope at low N (distinct~cN), robust: top decile ratio
    up=method_A(x,y)
    if up.sum()<5: return np.zeros(len(x),bool)
    c=np.median(D[up]/N[up]); c=min(c,1.0)
    exp=G*(1-np.exp(-c*N/G)); exp=np.maximum(exp,1)
    return (np.log10(D)-np.log10(exp))>-0.4

show=[('Puccinia graminis','Triticum aestivum'),('Bipolaris sorokiniana','Triticum aestivum'),
      ('Puccinia striiformis','Triticum aestivum'),('Phakopsora pachyrhizi','Triticum aestivum')]
rows=[]
for sp,h in show:
    ps=pts[(sp,h)]
    if len(ps)<15: continue
    x=np.log10(np.array([r for r,_,_ in ps])); y=np.log10(np.array([d for _,d,_ in ps]))
    tid=ps[0][2]
    A=method_A(x,y); B=method_B(x,y); C=method_C(x,y,tid)
    for i in range(len(ps)):
        for m,cls in (('A line',A[i]),('B cluster',B[i]),('C rarefaction',C[i])):
            rows.append([f'{sp} / {h.split()[0]}',m,ps[i][0],ps[i][1],'real' if cls else 'crash-out'])
with open('/tmp/three.tsv','w',newline='') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(['species','method','reads','distinct','cls']); w.writerows(rows)
# summary
import collections as C2
summ=C2.Counter()
for r in rows: summ[(r[0],r[1],r[4])]+=1
for sp,h in show:
    lab=f'{sp} / {h.split()[0]}'
    line=[f'{lab[:34]:<36}']
    for m in ('A line','B cluster','C rarefaction'):
        rl=summ[(lab,m,'real')]; tot=rl+summ[(lab,m,'crash-out')]
        line.append(f'{m}: {rl}/{tot} ({100*rl//max(tot,1)}%)')
    print('  '.join(line))
