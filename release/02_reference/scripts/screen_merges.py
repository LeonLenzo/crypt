#!/usr/bin/env python3
"""Resolve NCBI provisional names onto the USDA binomial they belong to.

    python release/02_reference/scripts/screen_merges.py

Reads  data/curation/screen_species.tsv    the USDA master, resolved to NCBI taxids
       data/raw/raw_taxdmp.zip             names.dmp + nodes.dmp
       data/raw/raw_assembly_summary.txt   to count what each merge actually moves
Writes data/curation/screen_merges.tsv     one row per provisional node, with its action

NCBI and USDA can hold the same fungus under two names, and then the reference sees two
organisms. *Cercospora flagellaris* is the case that exposed it: USDA has the binomial with
33 host records and NO genome, while every genome is filed under the provisional node
*Cercospora* cf. *flagellaris*, which USDA never names. So the name arrives as an
observation with nothing to build, and the genome arrives as a congener.

The criterion is a token match on the binomial, never a substring: the NCBI name must parse
as `<Genus> cf. <epithet>` (a trailing strain code is allowed and ignored), and
`<Genus> <epithet>` must appear as a `usda_name` in the master. Genus and epithet are
compared as whole tokens, so `Cercospora cf. flagellaris` reaches `Cercospora flagellaris`
and nothing else.

    merge     the USDA binomial resolves to its own NCBI taxid -> move the provisional
              node's assemblies onto it. The two names become one taxon.
    reroute   USDA names the binomial but NCBI has no node for it, so the provisional node
              is the only node this fungus has -> keep its taxid, promote it to an
              observation.

`aff.` IS DELIBERATELY NOT MERGED. `cf.` ("compare with") says the isolate is probably that
species; `aff.` ("has affinity to") says it is probably a DIFFERENT, undescribed species
related to it. Merging an `aff.` node would assert the identity its name denies, and an
undescribed relative is exactly what the congener tier is for. They are written out with
action=`keep` so the decision is visible rather than silent.

Everything here is reversible: the table is the whole record, and dropping a row un-merges
that taxon on the next run.
"""

import csv
import io
import re
import zipfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
SPECIES = HERE / "data/curation/screen_species.tsv"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
SUMMARY = HERE / "data/raw/raw_assembly_summary.txt"
OUT = HERE / "data/curation/screen_merges.tsv"
csv.field_size_limit(10 ** 7)

# <Genus> <qualifier>. <epithet> [strain code...]  -- qualifier captured so aff. is visible
PROVISIONAL = re.compile(r"^([A-Z][a-z]+) (cf|aff|nr)\. ([a-z][a-z-]+)(?:\s+.*)?$")
MERGEABLE = {"cf", "nr"}


def load_names():
    """-> name[taxid], rank[taxid]. Scientific names only."""
    name, rank = {}, {}
    with zipfile.ZipFile(TAXDMP) as zf:
        with zf.open("nodes.dmp") as raw:
            for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
                p = [f.strip() for f in ln.split("|")]
                rank[int(p[0])] = p[2]
        with zf.open("names.dmp") as raw:
            for ln in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
                p = [f.strip() for f in ln.split("|")]
                if len(p) >= 4 and p[3] == "scientific name":
                    name[int(p[0])] = p[1]
    return name, rank


def main():
    name, rank = load_names()

    # the USDA master, keyed by its binomial
    usda = {}
    with open(SPECIES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = " ".join(r["usda_name"].split()[:2])
            # prefer a row that resolved to a taxid over one that did not
            if key not in usda or (not usda[key]["ncbi_taxid"] and r["ncbi_taxid"]):
                usda[key] = r

    n_asm = defaultdict(int)
    with open(SUMMARY, encoding="utf-8", errors="replace") as fh:
        hdr = None
        for ln in fh:
            if ln.startswith("#assembly_accession"):
                hdr = ln[1:].rstrip("\n").split("\t")
                continue
            if ln.startswith("#") or hdr is None:
                continue
            d = dict(zip(hdr, ln.rstrip("\n").split("\t")))
            if d["species_taxid"].isdigit():
                n_asm[int(d["species_taxid"])] += 1

    rows = []
    for tid, nm in name.items():
        m = PROVISIONAL.match(nm)
        if not m:
            continue
        genus, qual, epithet = m.groups()
        hit = usda.get(f"{genus} {epithet}")
        if not hit:
            continue
        if not n_asm[tid]:
            continue                       # nothing to move; not worth a curation row
        to_tid = hit["ncbi_taxid"]
        if qual not in MERGEABLE:
            action, to_tid, note = "keep", "", f"{qual}. denotes a related undescribed species"
        elif to_tid.isdigit() and int(to_tid) != tid:
            action, note = "merge", "USDA binomial has its own NCBI node"
        elif not to_tid:
            action, note = "reroute", "USDA names it; NCBI has no node for the binomial"
        else:
            continue                       # already the same node
        rows.append({
            "from_taxid": tid, "from_name": nm, "qualifier": f"{qual}.",
            "action": action, "to_taxid": to_tid,
            "to_name": name.get(int(to_tid), "") if str(to_tid).isdigit() else "",
            "usda_name": hit["usda_name"], "usda_records": hit["n_records"],
            "to_has_own_assemblies": ("yes" if str(to_tid).isdigit() and n_asm[int(to_tid)]
                                      else "no"),
            "assemblies_moved": n_asm[tid], "note": note,
        })

    rows.sort(key=lambda r: (r["action"], r["from_name"]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["from_taxid", "from_name", "qualifier", "action", "to_taxid", "to_name",
            "usda_name", "usda_records", "to_has_own_assemblies", "assemblies_moved", "note"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    by = defaultdict(list)
    for r in rows:
        by[r["action"]].append(r)
    for act in ("merge", "reroute", "keep"):
        rs = by.get(act, [])
        print(f"{act:8s} {len(rs):>3d} taxa, {sum(r['assemblies_moved'] for r in rs):>3d} assemblies")
        for r in rs:
            tgt = r["to_name"] or r["usda_name"]
            flag = "" if r["to_has_own_assemblies"] == "yes" else "   <- target had no genome"
            print(f"         {r['from_name']:40.40s} -> {tgt:32.32s} "
                  f"{r['assemblies_moved']:>3d} asm{flag}")
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
