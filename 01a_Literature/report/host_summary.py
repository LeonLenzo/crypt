#!/usr/bin/env python3
"""Per-host summary of the field cohort, with the species taken per RUN.

Host is NOT read from `bioprojects.tsv`. That column holds one organism per BioProject and it
has been wrong twice in a single day's curation: PRJNA753546 registers seven plant species
plus aphids as `Zea mays`, and PRJNA630305 registers five tree species as `Fagus grandifolia`.
A per-host table built from it would be confidently wrong about roughly a thousand samples.

Instead `scientific_name` is fetched per run from the ENA portal, one request per BioProject
(see crypt-ena-portal-fetch), and cached. Where a study's own curation already resolved the
species per sample, that wins: the Lappe maize/teosinte survey carries a Species column in
its join table, taken from the paper's Supplementary Table 1.

The cohort is defined exactly as it is elsewhere:
    setting starts with `field`      the tissue was field-grown
    library_representative != no    one row per biological sample
    tissue is not root-like         aerial tissue only, and not the 3 insect runs

Prints a readable summary and writes `data/host_summary.tsv`.

Run:  python 01a_Literature/report/host_summary.py
      python 01a_Literature/report/host_summary.py --refetch
"""

from __future__ import annotations

import argparse, collections, csv, json, re, sys, urllib.request
from pathlib import Path

# Entry point in a subdirectory. 01a_Literature is not a valid package name (it starts with
# a digit) so `python -m` is unavailable; this puts the module root and the repo root on the
# path before any local import, which must therefore come after it.
_HERE = Path(__file__).resolve()
sys.path[:0] = [str(_HERE.parents[1]), str(_HERE.parents[2])]

from _layout import BIOSAMPLE, HOST_SUMMARY, RUNS, RUN_SPECIES, STUDIES

from _hostgroup import host_group
from _cohort import all_runs, cohort

ROOT = Path(__file__).resolve().parents[2]
CACHE = RUN_SPECIES
OUT = HOST_SUMMARY

PORTAL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
          "?accession={}&result=read_run&format=tsv&fields=run_accession,scientific_name")

# NON_AERIAL and the rest of the cohort predicate moved to cohort.py (2026-10-06):
# it was duplicated here and in every ad-hoc query, and the copies drifted.

# Tidy the long tail of synonyms and infraspecific names into one label per host, so the table
# groups by the plant a reader would name. Subspecies are kept where they matter biologically:
# teosinte is not maize for this study's purposes, and wild grape is not cultivated grape.
TIDY = [
    (r"^Oryza sativa.*",                     "Oryza sativa (rice)"),
    (r"^Oryza rufipogon",                    "Oryza rufipogon (wild rice)"),
    (r"^Zea mays subsp\.? mays|^Zea mays$",  "Zea mays (maize)"),
    (r"^Zea mays (ssp|subsp)\.? (mexicana|parviglumis)", "Zea mays teosinte subspp."),
    (r"^Zea (nicaraguensis|diploperennis|luxurians|perennis)", "Zea spp. (wild teosinte)"),
    (r"^Triticum aestivum",                  "Triticum aestivum (bread wheat)"),
    # A bare "Triticum" comes from BioSample host fields that record only the genus.
    (r"^Triticum$|^Triticum sp",             "Triticum aestivum (bread wheat)"),
    (r"^Triticum turgidum|^Triticum durum",  "Triticum durum/turgidum"),
    (r"^Vitis vinifera subsp\.? sylvestris", "Vitis vinifera subsp. sylvestris (wild grape)"),
    (r"^Vitis vinifera",                     "Vitis vinifera (grapevine)"),
    (r"^Vitis (riparia|amurensis|davidii)",  "Vitis spp. (wild grapevine)"),
    (r"^Arabidopsis halleri",                "Arabidopsis halleri"),
    (r"^Arabidopsis thaliana",               "Arabidopsis thaliana"),
    (r"^Sorghum bicolor",                    "Sorghum bicolor"),
    (r"^Brassica napus",                     "Brassica napus (oilseed rape)"),
    (r"^Populus nigra",                      "Populus nigra (black poplar)"),
    (r"^Fagus grandifolia",                  "Fagus grandifolia (American beech)"),
    (r"^Helianthus annuus",                  "Helianthus annuus (sunflower)"),
    (r"^Secale cereale",                     "Secale cereale (rye)"),
    (r"^Triticale|^x?Triticosecale",         "Triticale"),
    (r"^Panicum virgatum|^Switchgrass",      "Panicum virgatum (switchgrass)"),
]


def tidy(name: str) -> str:
    n = (name or "").strip()
    if not n:
        return "(unknown)"
    for pat, label in TIDY:
        if re.match(pat, n, re.I):
            return label
    return n


def fetch_species(bps: list[str], refetch: bool) -> dict:
    cached = {}
    if CACHE.exists() and not refetch:
        cached = {r["Run"]: r["scientific_name"]
                  for r in csv.DictReader(CACHE.open(), delimiter="\t")}
    need = [b for b in bps if not any(True for _ in ())] if refetch else None
    out = dict(cached)
    for i, bp in enumerate(bps, 1):
        if not refetch and any(True for r in cached):
            pass
        # Only fetch a project we have no rows for at all.
        if not refetch and cached and bp in getattr(fetch_species, "_seen", set()):
            continue
        try:
            with urllib.request.urlopen(PORTAL.format(bp), timeout=90) as fh:
                lines = fh.read().decode().splitlines()
        except Exception as e:                      # noqa: BLE001
            print(f"  {bp}: portal failed ({e})", file=sys.stderr)
            continue
        n = 0
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
                n += 1
        print(f"  [{i}/{len(bps)}] {bp}: {n} runs", file=sys.stderr)
    cols = ["Run", "scientific_name"]
    with CACHE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for run, sp in sorted(out.items()):
            w.writerow({"Run": run, "scientific_name": sp})
    return out


# Projects where the archive's scientific_name is not the HOST, or is absent, and the paper
# settles it. Both cases below were found by running this script and reading its output.
WHOLE_PROJECT_HOST = {
    # The Sato Arabidopsis field cluster: the ENA portal does not mirror these NCBI projects,
    # so every run came back blank. The paper is explicit - "199 A. thaliana accessions".
    "PRJNA1055060": "Arabidopsis thaliana", "PRJNA1055734": "Arabidopsis thaliana",
    "PRJNA1055104": "Arabidopsis thaliana", "PRJNA1056126": "Arabidopsis thaliana",
    "PRJNA1055317": "Arabidopsis thaliana", "PRJNA1055736": "Arabidopsis thaliana",
    "PRJNA1055424": "Arabidopsis thaliana", "PRJNA1055939": "Arabidopsis thaliana",
    # Rust surveys whose archive name is the pathogen and whose host is stated elsewhere.
    # PRJEB65589 from Lewis et al. 2024 Notes S1 ("Pgt-infected bread wheat (Triticum
    # aestivum) plants collected in the UK..."); PRJNA486288 from its BioSample host field
    # ("Soft red winter wheat: IL11-28222"); PRJNA1231252 from its cultivar field (Morocco,
    # Amboise, Benchmark, Kalmar are all wheat lines).
    "PRJEB65589": "Triticum aestivum", "PRJNA486288": "Triticum aestivum",
    "PRJNA1231252": "Triticum aestivum",
    # Resolved from the papers on 2026-10-08, after leon asked the obvious question: if we
    # hold the paper, how can the crop be unknown? It could not. Each of these states the
    # host outright and the count matches the deposit.
    #   PRJEB84390  "A total of 100 Pt-infected wheat field samples were subjected to RNA
    #               extraction"                                (100 runs) 10.1186/s12864-025-12230-4
    #   PRJEB36485  "RNA-seq analysis of Pst-infected wheat tissue ... across wheat-growing
    #               regions within South Africa"                (49 runs) 10.1111/ppa.13468
    #   PRJEB47693  "bread wheat (Triticum aestivum), durum wheat (T. durum) or barley
    #               (Hordeum vulgare) stem or leaf samples"     10.1111/ppa.13532 — per-sample
    #               hosts already come from its BioSample attributes; this covers the one run
    #               whose attribute names the fungus instead.
    "PRJEB84390": "Triticum aestivum", "PRJEB36485": "Triticum aestivum",
    "PRJEB47693": "Triticum aestivum",
    #   PRJNA1181034  no readable paper, but its own BioProject Description says it twice:
    #                 "Reduces Fusarium Head Blight by activating wheat resistance ... its
    #                 impact on the wheat transcriptome."
    "PRJNA1181034": "Triticum aestivum",
}

# A name matching this is the PATHOGEN, not the host, so WHOLE_PROJECT_HOST must outrank it.
PATHOGEN_NAME = re.compile(r"Puccinia|Blumeria|Zymoseptoria|Fusarium|Pyrenophora|Rhizoctonia"
                           r"|Magnaporthe|Pyricularia|Ustilaginoidea|Colletotrichum|^PDA$", re.I)


def curated_species() -> dict:
    """Per-sample species a study's own curation already resolved. These override the archive.

    Two reasons the archive cannot be trusted here, both demonstrated by this cohort:

    `scientific_name` is not always the HOST. Ada21's rust surveillance samples are infected
    WHEAT LEAVES registered under `Puccinia striiformis f. sp. tritici`, because the submitter's
    subject was the pathogen. Run naively, this script reported 929 samples whose host was a
    rust fungus. Table S1 carries the real host per sample (T. aestivum 916, Triticale 10,
    T. durum 8, Secale cereale 3) and ada21_runs.csv preserves it.

    And it is sometimes simply absent: the ENA portal does not mirror the eight NCBI projects
    of the Sato Arabidopsis cluster, so all 2,398 came back blank.
    """
    out = {}
    for rel, col in (("doi_10.1186_s12864-022-09001-w/lappe_runs.csv", "Species"),
                     ("doi_10.1186_s12864-021-07488-3/ada21_runs.csv", "HostSpecies"),
                     # Lewis 2024 Dataset S2: Wheat 50, Rye 2. WHOLE_PROJECT_HOST had folded
                     # the whole of PRJEB65589 into wheat, which was right for 50 of 52.
                     ("doi_10.1111_nph.19864/lewis24_runs.csv", "HostSpecies")):
        p = STUDIES / rel
        if not p.exists():
            continue
        for r in csv.DictReader(p.open()):
            if r.get(col):
                out[r["Run"]] = r[col]
    return out


def biosample_host() -> dict:
    """Per-run host from the BioSample `host` attribute, via its BioSample accession.

    Added 2026-10-08. `scientific_name` names the PATHOGEN for every rust and FHB survey in
    this cohort, which put 301 cohort runs under a fungus. For four of those eight projects
    the submitters recorded the real host as a BioSample attribute and nobody was reading it:
    PRJEB47693 gives Triticum aestivum 28 / Secale cereale 4 / Hordeum vulgare 1 per sample,
    and PRJNA950118 gives Triticum aestivum 54 / Hordeum jubatum 4. Those two carry the only
    field BARLEY in the cohort, so skipping this attribute reported barley as absent.

    Ranks below curated_species(), which comes from a paper's own table, and above the
    archive's scientific_name.
    """
    out = {}
    for r in all_runs(RUNS):
        bs = r.get("BioSample")
        if not bs:
            continue
        a = _attrs_cache.setdefault(r["BioProject"], _load_attrs(r["BioProject"]))
        h = (a.get(bs) or {}).get("host") if isinstance(a, dict) else None
        if h and not re.match(r"Puccinia|Blumeria|Zymoseptoria|Fusarium|Pyrenophora|Rhizoctonia",
                              str(h), re.I):
            out[r["Run"]] = str(h)
    return out


_attrs_cache: dict = {}


def _load_attrs(bp: str) -> dict:
    p = BIOSAMPLE / f"{bp}.json"
    if not p.exists():
        return {}
    try:
        d = json.loads(p.read_text())
    except Exception:
        return {}
    return d.get("attrs") if isinstance(d.get("attrs"), dict) else d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refetch", action="store_true")
    ap.add_argument("--per-run", metavar="TSV",
                    help="also write one row per cohort run with its resolved host, for "
                         "figures that need the host but must not re-derive it: the archive "
                         "name is the PATHOGEN for the rust surveys, so any script reading "
                         "run_species.tsv directly reports 56 wheat runs instead of 959")
    args = ap.parse_args()

    # The cohort predicate lives in cohort.py and nowhere else; see its docstring for why.
    runs = cohort(all_runs(RUNS))
    bps = sorted({r["BioProject"] for r in runs})
    print(f"cohort: {len(runs):,} samples across {len(bps)} projects", file=sys.stderr)

    have = {}
    if CACHE.exists() and not args.refetch:
        have = {r["Run"]: r["scientific_name"]
                for r in csv.DictReader(CACHE.open(), delimiter="\t")}
    missing = sorted({r["BioProject"] for r in runs if r["Run"] not in have})
    if missing:
        print(f"fetching species for {len(missing)} project(s)", file=sys.stderr)
        have = fetch_species(missing, refetch=False) | have if have else fetch_species(missing, False)
    have.update(biosample_host())
    have.update(curated_species())      # a paper's own table outranks the archive attribute

    rows = []
    for r in runs:
        sp = have.get(r["Run"]) or ""
        # WHOLE_PROJECT_HOST was written as a fallback for a MISSING name, so it never fired
        # for the rust surveys: those runs do have a name, it is just the fungus's. Adding
        # PRJEB65589 to the table on 2026-10-08 therefore changed nothing until this check
        # existed. A stated host beats a pathogen name.
        if not sp or PATHOGEN_NAME.search(sp):
            sp = WHOLE_PROJECT_HOST.get(r["BioProject"], "") or sp
        rows.append(dict(r, host=tidy(sp)))
    unknown = sum(1 for r in rows if r["host"] == "(unknown)")

    if args.per_run:
        cols = ["Run", "BioProject", "host", "location", "collection_date", "setting",
                "sampling_selection", "tissue"]
        with open(args.per_run, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore", quoting=csv.QUOTE_NONE, escapechar=None)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.per_run}  ({len(rows):,} rows)", file=sys.stderr)

    agg = collections.defaultdict(lambda: dict(
        n=0, projects=set(), locs=set(), years=set(), settings=collections.Counter(),
        sel=collections.Counter(), tissues=collections.Counter()))
    for r in rows:
        a = agg[r["host"]]
        a["n"] += 1
        a["projects"].add(r["BioProject"])
        if r["location"]:
            a["locs"].add(r["location"])
        if r["collection_date"][:4].isdigit():
            a["years"].add(r["collection_date"][:4])
        a["settings"][r["setting"]] += 1
        a["sel"][r["sampling_selection"] or "(none)"] += 1
        a["tissues"][r["tissue"]] += 1

    # Two evidence standards, not one denominator and a discard pile. See triage.py's
    # sampling_selection comment: in samples with a KNOWN pathogen a finding is a pathogen
    # beyond it, and in samples with none a single unreported pathogen is already the finding.
    # Fungicide-treated samples count toward a numerator but their absences are uninformative.
    COINFECTION = {"disease-selected", "inoculated"}
    DISCOVERY = {"unselected", "symptom-avoided"}
    out = []
    from _hostgroup import ORDER
    for host, a in sorted(agg.items(),
                          key=lambda kv: (ORDER[host_group(kv[0])], -kv[1]["n"])):
        yrs = sorted(a["years"])
        out.append(dict(
            host=host, host_group=host_group(host), samples=a["n"], projects=len(a["projects"]),
            localities=len(a["locs"]),
            years=f"{yrs[0]}-{yrs[-1]}" if len(yrs) > 1 else (yrs[0] if yrs else ""),
            n_years=len(yrs),
            wild=sum(v for k, v in a["settings"].items() if "natural population" in k),
            coinfection_set=sum(v for k, v in a["sel"].items() if k in COINFECTION),
            discovery_set=sum(v for k, v in a["sel"].items() if k in DISCOVERY),
            positives_only=sum(v for k, v in a["sel"].items() if k == "fungicide-treated"),
            unflagged=sum(v for k, v in a["sel"].items() if k == "(none)"),
            main_tissue=a["tissues"].most_common(1)[0][0],
            settings="; ".join(f"{k}:{v}" for k, v in a["settings"].most_common()),
            selection="; ".join(f"{k}:{v}" for k, v in a["sel"].most_common())))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print(f"\n{'host':<40}{'samp':>6}{'proj':>5}{'loc':>5}{'years':>12}{'wild':>6}"
          f"{'coinf':>7}{'disc':>7}  tissue")
    print("-" * 108)
    for r in out:
        print(f"{r['host'][:39]:<40}{r['samples']:>6}{r['projects']:>5}{r['localities']:>5}"
              f"{r['years']:>12}{r['wild'] or '':>6}{r['coinfection_set'] or '':>7}"
              f"{r['discovery_set'] or '':>7}  {r['main_tissue'][:22]}")
    print("-" * 108)
    print(f"{'TOTAL':<40}{sum(r['samples'] for r in out):>6}"
          f"{len({p for a in agg.values() for p in a['projects']}):>5}"
          f"{len({l for a in agg.values() for l in a['locs']}):>5}")
    if unknown:
        print(f"\n{unknown} samples have no species from the archive", file=sys.stderr)
    print(f"\nwrote {OUT.relative_to(ROOT)}", file=sys.stderr)


if __name__ == "__main__":
    main()
