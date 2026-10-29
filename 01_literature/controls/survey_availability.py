import json, time, urllib.parse, urllib.request
PATHOGENS = ["Puccinia striiformis","Puccinia triticina","Puccinia graminis","Puccinia sorghi",
 "Blumeria graminis","Zymoseptoria tritici","Parastagonospora nodorum","Exserohilum turcicum",
 "Bipolaris zeicola","Bipolaris maydis","Bipolaris sorokiniana","Cercospora zeae-maydis",
 "Pyrenophora tritici-repentis","Pyrenophora teres","Curvularia clavata","Fusarium graminearum",
 "Fusarium verticillioides","Alternaria alternata","Alternaria tenuissima"]
HOSTS = ["Triticum aestivum","Zea mays","Sorghum bicolor","Oryza sativa","Hordeum vulgare"]
def get(u):
    r = urllib.request.Request(u, headers={"User-Agent":"crypt-controls-survey"})
    with urllib.request.urlopen(r, timeout=90) as fh: return fh.read().decode()
print("kind\tspecies\ttaxid\trank\trna_seq_runs")
for kind, names in (("pathogen",PATHOGENS),("host",HOSTS)):
    for s in names:
        try:
            tx=json.loads(get("https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"+urllib.parse.quote(s)))
            tid,rank=tx[0]["taxId"],tx[0].get("rank","")
        except Exception:
            print(f"{kind}\t{s}\t-\t-\tTAXON_LOOKUP_FAILED",flush=True); time.sleep(0.5); continue
        time.sleep(0.5)
        q=urllib.parse.quote(f'tax_tree({tid}) AND library_strategy="RNA-Seq"')
        try:
            n=get(f"https://www.ebi.ac.uk/ena/portal/api/count?result=read_run&query={q}&dataPortal=ena").strip().split("\n")[-1]
        except Exception: n="COUNT_FAILED"
        print(f"{kind}\t{s}\t{tid}\t{rank}\t{n}",flush=True)
        time.sleep(0.5)
