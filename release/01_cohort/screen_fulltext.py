#!/usr/bin/env python3
"""Pull the study-level facts out of each paper's text: what was sequenced, from where, and
what it was infected with.

    python release/01_cohort/screen_fulltext.py

Reads  raw/fulltext.tsv          every text we hold, and its source
       inputs/pathogen_genera.txt a flag-list of cereal pathogen genera
Writes curation/screen_fulltext.tsv   one row per paper

This replaces bioproject_screen.py, which found only accessions. The corpus is now field
pathogenomic wheat studies by construction, so host, setting and disease-selection no longer
need curating away a mixed population; they are captured here ONLY to verify the search did
not let something through, and are expected to be constant.

WHAT EACH COLUMN IS FOR
  accessions, accession_types   the deposits to resolve. Kept apart by type: only PRJ* is a
                                BioProject; GSE/SRP/E-MTAB resolve to one; SRR/SAM imply one.
  host, host_ok                 wheat is assumed; host_ok is 0 when a non-target crop (barley,
                                the yellow-dwarf-virus surveys) dominates, i.e. a search miss.
  field                         open-field language present. Expected true throughout; a false
                                is a paper to look at.
  tissue                        leaf/spike/stem/root, for the foliar gate downstream.
  location, years               the survey's geography and seasons: the site x year structure
                                the whole analysis rests on.
  target_pathogens              the organism(s) the STUDY declares it is about.
  mentioned_pathogens           EVERY pathogen genus named anywhere in the text. The point of
                                reading full text: "undeclared" must mean "the authors never
                                named it", not "not the study's subject", or the central claim
                                is wrong on any paper that remarks a second organism in passing.
  rna_selection                 polyA / rRNA-depleted / total / small. Drives how much non-host
                                sequence a library holds, and catches small-RNA studies that do
                                not belong.
  healthy_sampled               whether asymptomatic plants were also sequenced: an internal
                                comparison group.
  axenic_control                whether a pure culture was sequenced: a free cross-mapping
                                control.
  n_samples_claimed, n_sites_claimed   the paper's own counts, to check the archive against.
  supp_for                      fields NOT found in the body, whose value the paper says is in
                                a named supplementary file. The flag to go mine the supplement.
  text_source                   pmc / biorxiv / manual / abstract. An abstract-only row is a
                                thin screen and is marked as such, not hidden.
"""

import collections
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULLTEXT = HERE / "raw" / "fulltext.tsv"
ABSTRACT_SCREEN = HERE / "curation" / "screen_abstract.tsv"
DECISIONS = HERE / "curation" / "decisions.tsv"
GENERA = HERE / "inputs" / "pathogen_genera.txt"
OUT = HERE / "curation" / "screen_fulltext.tsv"

TAG = re.compile(r"<[^>]+>")
ACCESSION = re.compile(
    r"(?P<bioproject>PRJ(?:NA|EB|DB)\d{4,})|(?P<sra_study>\b[SED]RP\d{5,})|"
    r"(?P<geo>\bGSE\d{3,})|(?P<arrayexpress>\bE-[A-Z]{4}-\d+)|"
    r"(?P<run>\b[SED]RR\d{5,})|(?P<biosample>\bSAM[NED][A-Z]?\d{5,})")
DEPOSIT_SEC = re.compile(r"data availab|availability of data|accession|data deposition|"
                         r"sequence read archive|deposited", re.I)

TARGET_HOST = re.compile(r"\bwheat\b|triticum|durum|spelt|triticale", re.I)
OTHER_HOST = re.compile(r"\bbarley\b|hordeum|\bmaize\b|\bcorn\b|zea mays|\brice\b|oryza|"
                        r"\boat\b|avena|\brye\b|secale|sorghum", re.I)
FIELD = re.compile(r"\bfield[- ](grown|collect|sampl|condition|site|trial|experiment|surve)|"
                   r"\bin the field\b|naturally[- ]infect|farmers?[’' ]? ?field|nurser|"
                   r"commercial (?:field|crop)|growing season", re.I)
GLASSHOUSE = re.compile(r"glasshouse|greenhouse|growth chamber|controlled[- ]environment|"
                        r"growth room|climate chamber", re.I)
TISSUE = {"leaf": r"\bleaf|leaves|flag[- ]leaf|foliar",
          "spike": r"\bspike|\bhead\b|\bear\b|glume|spikelet|panicle|grain|kernel|seed\b",
          "stem": r"\bstem\b|\bculm\b|\bstalk\b|internode|sheath",
          "root": r"\broot|rhizosphere|crown"}
RNA_SEL = [("small", r"small[- ]?rna|srna|micro[- ]?rna|mirna|srna[- ]seq"),
           ("rRNA-depleted", r"ribo[- ]?(?:zero|minus)|rrna[- ]?deplet|ribo[- ]?deplet|"
                             r"total rna.{0,30}deplet"),
           ("polyA", r"poly[- ]?a\b|polyadenylat|mrna[- ]?seq|oligo[- ]?dt|poly\(a\)"),
           ("total", r"\btotal rna\b")]
HEALTHY = re.compile(r"asymptomatic|healthy (?:leaf|leaves|plant|tissue|control)|"
                     r"non[- ]infected|uninfected|symptomless", re.I)
AXENIC = re.compile(r"axenic|pure cultur|in vitro cultur|single[- ]spore|isolate.{0,20}cultur|"
                    r"myceli.{0,30}cultur", re.I)
SUPP = re.compile(r"(?:sup(?:plementary|plemental|p)?\.?\s*(?:table|data|file|dataset|"
                  r"information|materials?)\s*S?\d+[A-Za-z]?)", re.I)
YEAR = re.compile(r"\b(19[89]\d|20[0-2]\d)\b")
SPECIES = re.compile(r"\b([A-Z][a-z]{3,})\s+([a-z]{3,})\b")
SEQ_METHOD = re.compile(
    r"rna[- ]?seq|transcriptom|metatranscriptom|total rna|poly[- ]?a|ribo[- ]?deplet|"
    r"small[- ]?rna|high[- ]?throughput sequenc|shotgun sequenc|illumina|nanopore|hiseq|"
    r"novaseq|miseq|cdna librar|rna[- ]?tag[- ]?seq|quant[- ]?seq", re.I)


def sections(xml):
    out = []
    for m in re.finditer(r"<sec\b[^>]*>(.*?)</sec>", xml, re.S | re.I):
        block = m.group(1)
        t = re.search(r"<title[^>]*>(.*?)</title>", block, re.S | re.I)
        title = TAG.sub(" ", t.group(1)) if t else ""
        out.append((" ".join(title.split()), TAG.sub(" ", block)))
    out.append(("", TAG.sub(" ", xml)))
    return out


def plain(text, is_html):
    if is_html:
        text = re.sub(r"<(script|style)\b.*?</\1>", " ", text, flags=re.S | re.I)
        text = TAG.sub(" ", text)
    return " ".join(text.split())


def first(rx, text, n=1):
    found = []
    pat = rx if isinstance(rx, re.Pattern) else re.compile(rx, re.I)
    for m in pat.finditer(text):
        v = m.group(0)
        if v not in found:
            found.append(v)
        if len(found) >= n:
            break
    return found


def scan_accessions(text, is_sectioned, doi, rows):
    seen = set()
    blocks = sections(text) if is_sectioned else [("", plain(text, "<" in text[:400]))]
    for title, body in blocks:
        for m in ACCESSION.finditer(" ".join(body.split())):
            acc, kind = m.group(0), m.lastgroup
            if acc in seen:
                continue
            seen.add(acc)
            a, b = max(0, m.start() - 80), m.end() + 80
            rows.append({"doi": doi, "accession": acc, "type": kind,
                         "section": title or "unsectioned",
                         "in_deposit_section": int(bool(DEPOSIT_SEC.search(title))),
                         "context": " ".join(body[a:b].split())[:160]})
    return seen


def main():
    genera = {l.strip() for l in GENERA.read_text().splitlines()
              if l.strip() and not l.startswith("#")}
    recommend = {}
    if ABSTRACT_SCREEN.exists():
        with open(ABSTRACT_SCREEN, newline="", encoding="utf-8") as fh:
            recommend = {r["doi"]: r["recommend"] for r in csv.DictReader(fh, delimiter="\t")}
    overrides = {}
    if DECISIONS.exists():
        with open(DECISIONS, newline="", encoding="utf-8") as fh:
            overrides = {r["doi"]: (r["decision"], r.get("reason", ""))
                         for r in csv.DictReader(fh, delimiter="\t") if r.get("decision")}

    rows, acc_rows = [], []
    with open(FULLTEXT, newline="", encoding="utf-8") as fh:
        papers = list(csv.DictReader(fh, delimiter="\t"))

    for r in papers:
        doi, source, path = r["doi"], r.get("source", ""), r.get("path", "")
        raw = (HERE / path).read_text(errors="ignore") if path and (HERE / path).is_file() else ""
        is_sectioned = "<sec" in raw[:200000]
        text = raw if is_sectioned else plain(raw, path.lower().endswith((".html", ".htm")))
        full = (r.get("title", "") + " " + text)

        accs = scan_accessions(raw or text, is_sectioned, doi, acc_rows) if text else set()
        bp = sorted(a for a in accs if a.startswith("PRJ"))
        own = sorted({x["accession"] for x in acc_rows if x["doi"] == doi
                      and x["accession"].startswith("PRJ") and x["in_deposit_section"] == 1})

        tissue = [k for k, rx in TISSUE.items() if re.search(rx, text, re.I)]
        rna = next((name for name, rx in RNA_SEL if re.search(rx, full, re.I)), "")
        # Every binomial that looks like a pathogen genus, deduplicated at genus.
        named = {}
        for m in SPECIES.finditer(text):
            if m.group(1) in genera:
                named.setdefault(m.group(1), f"{m.group(1)} {m.group(2)}")
        years = sorted(set(first(YEAR, text, 12)))

        # Which wanted fields are absent from the body but promised in a supplement.
        supp_for = []
        if text:
            if not accs and SUPP.search(text):
                supp_for.append("accessions")
            if not years and re.search(r"collect|season|sampl", text, re.I) and SUPP.search(text):
                supp_for.append("year")
            if not tissue and SUPP.search(text):
                supp_for.append("tissue")

        seq_in_body = bool(text and SEQ_METHOD.search(text))
        rows.append({
            "doi": doi, "text_source": source, "seq_in_body": int(seq_in_body),
            "bioprojects": ";".join(bp), "bioprojects_own": ";".join(own),
            "all_accessions": ";".join(sorted(accs)),
            "accession_types": ";".join(sorted({a.split("_")[0] for a in
                                                {x["type"] for x in acc_rows if x["doi"] == doi}})),
            "host": "wheat" if TARGET_HOST.search(full) else
                    (OTHER_HOST.search(full).group(0).lower() if OTHER_HOST.search(full) else ""),
            "host_ok": int(bool(TARGET_HOST.search(full))),
            "field": int(bool(FIELD.search(full))),
            "glasshouse_mentioned": int(bool(GLASSHOUSE.search(full))),
            "tissue": ";".join(tissue),
            "location_terms": "",   # placeholder; geographic parsing is its own step
            "years": ";".join(years),
            "target_pathogens": "",  # set below from title
            "mentioned_pathogens": ";".join(sorted(named.values())),
            "n_pathogen_genera_named": len(named),
            "rna_selection": rna,
            "healthy_sampled": int(bool(HEALTHY.search(text))),
            "axenic_control": int(bool(AXENIC.search(text))),
            "supp_for": ";".join(supp_for),
            "has_body": int(source in ("pmc", "biorxiv", "manual")),
        })
        # The target is whatever pathogen the TITLE names.
        tnamed = {m.group(1) for m in SPECIES.finditer(r.get("title", "")) if m.group(1) in genera}
        rows[-1]["target_pathogens"] = ";".join(sorted(
            v for g, v in named.items() if g in tnamed)) or ";".join(sorted(tnamed))

    for r in rows:
        doi = r["doi"]
        rec = recommend.get(doi, "keep")
        if doi in overrides:
            r["verdict"], r["verdict_reason"] = overrides[doi][0], overrides[doi][1]
        elif rec == "keep":
            r["verdict"], r["verdict_reason"] = "keep", "clear at the abstract"
        elif r["seq_in_body"] and (r["field"] == 1 or r["bioprojects"]):
            # A review-flagged paper whose BODY shows RNA-seq and either field material or a
            # deposit is a genuine study; the abstract was just terse.
            r["verdict"], r["verdict_reason"] = "keep", "body confirms RNA-seq of field material"
        elif r["has_body"] and not r["seq_in_body"] and not r["bioprojects"]:
            r["verdict"], r["verdict_reason"] = "drop", "body shows no RNA-seq and no deposit"
        else:
            r["verdict"], r["verdict_reason"] = "review", "body inconclusive, needs a human"

    COLS = ["doi", "text_source", "verdict", "verdict_reason", "seq_in_body",
            "has_body", "host", "host_ok", "field",
            "glasshouse_mentioned", "tissue", "years", "target_pathogens",
            "mentioned_pathogens", "n_pathogen_genera_named", "rna_selection",
            "healthy_sampled", "axenic_control", "bioprojects", "bioprojects_own",
            "accession_types", "all_accessions", "location_terms", "supp_for"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    v = collections.Counter(r["verdict"] for r in rows)
    print(f"-> {OUT.name}  {n} papers")
    print(f"   VERDICT: {dict(v.most_common())}  "
          f"({sum(1 for r in rows if r['doi'] in overrides)} by hand override)")
    still = [r["doi"] for r in rows if r["verdict"] == "review"]
    if still:
        print(f"   still need a human: {still}")
    print(f"   with a body (not abstract-only)   {sum(r['has_body'] for r in rows)}")
    print(f"   host verified wheat               {sum(r['host_ok'] for r in rows)}")
    print(f"   NON-wheat dominant (search miss?)  {[r['doi'] for r in rows if not r['host_ok']][:6]}")
    print(f"   field language present            {sum(r['field'] for r in rows)}")
    print(f"   names a BioProject                {sum(1 for r in rows if r['bioprojects'])}")
    print(f"   names >=1 pathogen genus          {sum(1 for r in rows if r['n_pathogen_genera_named'])}")
    print(f"   sampled healthy plants too        {sum(r['healthy_sampled'] for r in rows)}")
    print(f"   sequenced an axenic control       {sum(r['axenic_control'] for r in rows)}")
    rs = collections.Counter(r["rna_selection"] or "unstated" for r in rows)
    print(f"   RNA selection: {dict(rs.most_common())}")
    print(f"   need a supplement for some field  {sum(1 for r in rows if r['supp_for'])}")


if __name__ == "__main__":
    main()
