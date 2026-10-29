#!/usr/bin/env python3
"""Stage 1 and 2 of the control curation cascade: ENA metadata screen, then title screen.

Produces a SHORTLIST, never a decision. Admission needs the methods quote (stage 3) and the
purity screen on classified output (stage 4). See DESIGN.md.

The cues below are recorded per run in the `cue` column so a curator can see why a run was
shortlisted and overrule it. An in planta cue VETOES, because a library with host tissue is
disqualified whatever else its title says.
"""
import csv, json, os, re, sys, time, urllib.parse, urllib.request

OUT_RUNS = "candidates_runs.tsv"
OUT_STUDIES = "candidates_studies.tsv"
RAW = "raw"          # one cached ENA pull per taxon, so re-screening never re-fetches

# taxid: (species, control class wanted)
TAXA = {
    27350:  ("Puccinia striiformis", "spore"),
    208348: ("Puccinia triticina", "spore"),
    5297:   ("Puccinia graminis", "spore"),
    34373:  ("Blumeria graminis", "spore"),
    1047171:("Zymoseptoria tritici", "culture"),
    13684:  ("Parastagonospora nodorum", "culture"),
    93612:  ("Exserohilum turcicum", "culture"),
    5016:   ("Bipolaris maydis", "culture"),
    45130:  ("Bipolaris sorokiniana", "culture"),
    53485:  ("Pyrenophora teres", "culture"),
    45151:  ("Pyrenophora tritici-repentis", "culture"),
    95742:  ("Curvularia clavata", "culture"),
    5518:   ("Fusarium graminearum", "culture"),
    117187: ("Fusarium verticillioides", "culture"),
    5599:   ("Alternaria alternata", "culture"),
    119927: ("Alternaria tenuissima", "culture"),
}

FIELDS = ("run_accession,study_accession,sample_accession,study_title,experiment_title,"
          "sample_title,sample_description,description,host,host_scientific_name,"
          "isolation_source,dev_stage,strain,library_source,library_selection,"
          "read_count,base_count,instrument_platform")

# An in planta cue disqualifies. Host plant names are included because the sample taxon is
# whatever the submitter chose and the host often appears only in free text.
IN_PLANTA = [r"in[ _-]?planta", r"infect", r"inoculat", r"\bleaf\b", r"\bleaves\b", r"seedling",
             r"\bdpi\b", r"\bhpi\b", r"\bdai\b", r"lesion", r"greenhouse", r"glasshouse",
             r"\bfield\b", r"spikelet", r"\bgrain\b", r"\broot\b", r"\bstem\b", r"\bculm\b",
             r"\bglume\b", r"\bhead\b", r"wheat", r"barley", r"maize", r"sorghum", r"\brice\b",
             r"triticum", r"hordeum", r"\bzea\b", r"oryza", r"haustori", r"biotroph.{0,12}phase",
             r"symptom", r"pathogenesis", r"compatible interaction"]
SPORE = [r"urediniospore", r"uredospore", r"teliospore", r"basidiospore", r"pycnidiospore",
         r"conidi", r"\bspore", r"germinat", r"appressori"]
CULTURE = [r"in[ _-]?vitro", r"axenic", r"\bcultur", r"mycelia", r"mycelium", r"\bbroth\b",
           r"potato dextrose", r"\bpda\b", r"\bpdb\b", r"\bv8\b", r"czapek", r"minimal medium",
           r"complete medium", r"\bypd\b", r"malt extract", r"\bagar\b", r"fries", r"shake flask",
           r"submerged", r"liquid medium", r"\bsaprotroph", r"vegetative growth"]
# A medium made from the host disqualifies even though it reads as culture.
PLANT_MEDIUM = [r"straw", r"leaf extract", r"plant extract", r"wheat.{0,10}medium", r"host extract"]

def rx(pats, text):
    for p in pats:
        m = re.search(p, text)
        if m:
            return m.group(0).strip()
    return None

def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "crypt-controls-screen"})
            with urllib.request.urlopen(req, timeout=180) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))

# Fields that describe THIS library's material. Everything else is context.
#   study_title and experiment_title are STUDY level: a study comparing in vitro with in planta
#   mentions both, and letting its title vote vetoed all 8 `invitro_SN15_*` runs of PRJNA632579
#   on the phrase "in planta". Measured 2026-10-23, which is how this was caught.
#   host and isolation_source describe the ISOLATE, not the flask: a culture of a wheat-derived
#   isolate records "wheat" and would be vetoed wrongly. Same trap as crypt-isolate-vs-sample.
SAMPLE_FIELDS = ("sample_title", "sample_description", "description", "dev_stage")
CONTEXT_FIELDS = ("study_title", "experiment_title", "host", "host_scientific_name",
                  "isolation_source", "strain", "culture_collection")

def classify(row, wanted):
    text = " ".join(row.get(k, "") or "" for k in SAMPLE_FIELDS).lower()
    if row.get("library_source") != "TRANSCRIPTOMIC":
        return "not_rnaseq", row.get("library_source", "")
    bad = rx(IN_PLANTA, text)
    if bad:
        return "in_planta", bad
    med = rx(PLANT_MEDIUM, text)
    if med:
        return "plant_medium", med
    sp, cu = rx(SPORE, text), rx(CULTURE, text)
    if wanted == "spore":
        if sp: return "spore_candidate", sp
        if cu: return "culture_candidate", cu
    else:
        if cu: return "culture_candidate", cu
        if sp: return "spore_candidate", sp
    return "unknown", ""


def study_context(row):
    """Flag a study whose own text mentions host material. Not a veto. It tells the curator the
    study is probably mixed, so the methods read has to establish which arm a run belongs to."""
    text = " ".join(row.get(k, "") or "" for k in ("study_title", "experiment_title")).lower()
    hit = rx(IN_PLANTA, text)
    return f"study mentions {hit}" if hit else ""

def pmc_accession(acc):
    """Europe PMC full-text search on one accession."""
    j = json.loads(get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
                       + urllib.parse.quote(f'"{acc}"') + "&format=json&resultType=core&pageSize=3"))
    return j.get("resultList", {}).get("result", [])


def paper_for(study):
    """Three routes, in order of how much they can be trusted. All three are needed: measured
    2026-10-23 over 56 studies, the XREF alone found 6, adding the BioProject full-text search
    found 32 more, and the SRP/GEO search recovered studies that cite the GEO series instead of
    the BioProject, which is most of the older Fusarium and Zymoseptoria work."""
    alts = []
    try:
        xml = get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{study}")
        m = re.search(r"<DB>PUBMED</DB>\s*<ID>(\d+)</ID>", xml)
        if m:
            return "ena_xref", m.group(1), doi_for_pmid(m.group(1))
        srp = re.search(r"<SECONDARY_ID>(SRP\d+)</SECONDARY_ID>", xml)
        geo = re.search(r'namespace="GEO">(GSE\d+)<', xml)
        alts = [x.group(1) for x in (geo, srp) if x]
    except Exception:
        pass
    for acc, route in [(study, "europepmc")] + [(a, "europepmc_" + ("geo" if a.startswith("GSE") else "srp")) for a in alts]:
        time.sleep(0.5)
        try:
            res = pmc_accession(acc)
        except Exception:
            continue
        if len(res) == 1:
            return route, res[0].get("pmid", ""), res[0].get("doi", "")
        if len(res) > 1:
            return f"{route}_ambiguous_{len(res)}", res[0].get("pmid", ""), res[0].get("doi", "")
    return "none", "", ""

def doi_for_pmid(pmid):
    try:
        j = json.loads(get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}&format=json&resultType=core&pageSize=1"))
        r = j["resultList"]["result"]
        return r[0].get("doi", "") if r else ""
    except Exception:
        return ""

def main():
    runs = []
    for taxid, (sp, wanted) in TAXA.items():
        os.makedirs(RAW, exist_ok=True)
        cache = os.path.join(RAW, f"{taxid}.tsv")
        if os.path.exists(cache) and os.path.getsize(cache) > 0:
            tsv = open(cache).read()
        else:
            q = urllib.parse.quote(f'tax_tree({taxid}) AND library_strategy="RNA-Seq"')
            tsv = get(f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query={q}"
                      f"&fields={FIELDS}&format=tsv&limit=0")
            open(cache, "w").write(tsv)
            time.sleep(0.5)
        rdr = csv.DictReader(tsv.splitlines(), delimiter="\t")
        n = 0
        for row in rdr:
            verdict, cue = classify(row, wanted)
            row["species"] = sp; row["wanted"] = wanted
            row["verdict"] = verdict; row["cue"] = cue
            row["study_flag"] = study_context(row)
            row["isolate_origin"] = " / ".join(
                v for v in (row.get("host"), row.get("isolation_source")) if v)
            runs.append(row); n += 1
        print(f"{sp:<30} {n:>5} runs", file=sys.stderr, flush=True)

    cols = ["species", "wanted", "verdict", "cue", "run_accession", "study_accession",
            "sample_title", "strain", "read_count", "base_count", "instrument_platform",
            "study_flag", "isolate_origin", "study_title"]
    with open(OUT_RUNS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for r in runs:
            w.writerow(r)

    # study level, shortlisting only studies that could supply a control
    studies = {}
    for r in runs:
        k = (r["species"], r["study_accession"])
        s = studies.setdefault(k, {"species": r["species"], "wanted": r["wanted"],
                                   "study_accession": r["study_accession"],
                                   "study_title": r.get("study_title", ""),
                                   "n_runs": 0, "culture_candidate": 0, "spore_candidate": 0,
                                   "in_planta": 0, "plant_medium": 0, "unknown": 0,
                                   "not_rnaseq": 0, "gbases": 0.0})
        s["n_runs"] += 1
        s[r["verdict"]] = s.get(r["verdict"], 0) + 1
        try:
            s["gbases"] += float(r.get("base_count") or 0) / 1e9
        except ValueError:
            pass
    rows = sorted(studies.values(), key=lambda s: (s["species"],
                  -(s["culture_candidate"] + s["spore_candidate"])))
    short = [s for s in rows if (s["culture_candidate"] + s["spore_candidate"]) >= 2
             and s["study_accession"]]
    seen = {}
    if os.path.exists(OUT_STUDIES):
        for old in csv.DictReader(open(OUT_STUDIES), delimiter="\t"):
            if old.get("paper_route") and old["paper_route"] != "none":
                seen[old["study_accession"]] = (old["paper_route"], old["pmid"], old["doi"])
    print(f"\nresolving papers for {len(short)} shortlisted studies "
          f"({len(seen)} already cached)", file=sys.stderr, flush=True)
    for s in rows:
        s["paper_route"] = s["pmid"] = s["doi"] = ""
    for i, s in enumerate(short, 1):
        if s["study_accession"] in seen:
            s["paper_route"], s["pmid"], s["doi"] = seen[s["study_accession"]]
            continue
        s["paper_route"], s["pmid"], s["doi"] = paper_for(s["study_accession"])
        print(f"  [{i}/{len(short)}] {s['study_accession']} {s['paper_route']} {s['doi']}",
              file=sys.stderr, flush=True)
        time.sleep(0.6)
    scols = ["species", "wanted", "study_accession", "n_runs", "culture_candidate",
             "spore_candidate", "in_planta", "plant_medium", "unknown", "not_rnaseq",
             "gbases", "paper_route", "pmid", "doi", "study_title"]
    with open(OUT_STUDIES, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=scols, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for s in rows:
            s["gbases"] = f"{s['gbases']:.1f}"
            w.writerow(s)
    print(f"\nwrote {OUT_RUNS} ({len(runs)} runs) and {OUT_STUDIES} ({len(rows)} studies)",
          file=sys.stderr)

main()
