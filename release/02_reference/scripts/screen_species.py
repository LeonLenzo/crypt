#!/usr/bin/env python3
"""Assign every USDA fungus name an NCBI taxid, or record why it has none.

    python release/02_reference/scripts/screen_species.py

Reads  data/inputs/scope_rule.tsv        which name classes count as a match
       data/raw/raw_usda_fungus_host.csv the master list, 99,059 distinct fungus names
       data/raw/raw_taxdmp.zip           names.dmp + nodes.dmp, read from the zip
Writes data/curation/screen_species.tsv  one row per USDA name, with its status

Every name gets a row, including the ~69% NCBI has no node for. A reference absent from the
database is invisible and unrecoverable at classification time, because the reads go silently
to the nearest taxon present, so an unresolvable name is recorded rather than dropped.

`status` carries how the name resolved, and the distinctions in it are load-bearing:

  scientific_name   the current NCBI name
  synonym           NCBI files this name as a synonym of the matched node. 8,470 names
                    arrive this way and 2,902 taxa are reachable by NO other route, so
                    dropping this class loses taxa, not just spellings
  equivalent_name   an alternative name for the same node
  includes          loose subsumption rather than equivalence, so accepted but flagged
  *_trimmed         matched only after cutting an infraspecific tail back to the binomial.
                    Kept separate because the tail sometimes matters: dropping it from
                    `Ravenelia spinulosa var. spinulosa` is harmless, dropping it from a
                    forma specialis is not, and the declared set is written at f. sp. rank
  not_fungi         resolves to a real NCBI node outside Fungi. 1,326 names do, being
                    homonyms shared with animals and plants
  no_ncbi_node      no node under any accepted name class. USDA is a herbarium records
                    database reaching back to 1960s literature and most of those fungi have
                    never been sequenced, so this is a property of NCBI, not a match failure

NOT APPLIED HERE: the host side. The criterion in scope_rule.tsv says "recorded on at least
one plant host", and screening that needs the 64,961 host names resolved too. This script
counts hosts per fungus but does not test whether they are plants, so what it implements is
"recorded in USDA at all". The dataset is a plant-pathology compilation so the two nearly
coincide, but they are not the same statement and the difference is not yet measured.
"""

import argparse
import csv
import io
import zipfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
RULE = HERE / "data/inputs/scope_rule.tsv"
USDA = HERE / "data/raw/raw_usda_fungus_host.csv"
TAXDMP = HERE / "data/raw/raw_taxdmp.zip"
OUT = HERE / "data/curation/screen_species.tsv"
FUNGI = 4751

norm = lambda s: " ".join(s.split())


def read_rule():
    out = defaultdict(dict)
    with open(RULE) as fh:
        for row in csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                  delimiter="\t"):
            out[row["source"]][row["field"]] = row["value"]
    return out


def usda_names(crop_genera, aliases):
    """{name: (n_records, n_distinct_hosts, cohort host genera it is recorded on)}.

    The cohort-host set is a LABEL and never a filter: see scope_rule.tsv. Matching is on the
    host's GENUS, because USDA writes the host at every rank from `Triticum sp.` to
    `Zea mays var. saccharata`, and on the alias table, because a few records give a common
    name instead of a binomial.
    """
    rec = defaultdict(int)
    hosts = defaultdict(set)
    crops = defaultdict(set)
    with open(USDA, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            n = norm(row.get("sciName") or "")
            if not n:
                continue
            rec[n] += 1
            h = norm(row.get("host") or "")
            if not h:
                continue
            hosts[n].add(h)
            g = aliases.get(h.lower()) or h.split()[0]
            if g in crop_genera:
                crops[n].add(g)
    return {n: (rec[n], len(hosts[n]), crops[n]) for n in rec}


def dmp_rows(zf, member):
    """names.dmp and nodes.dmp are pipe-delimited with tab padding around each field."""
    with zf.open(member) as raw:
        for line in io.TextIOWrapper(raw, encoding="utf-8", errors="replace"):
            yield [f.strip() for f in line.split("|")]


def load_taxonomy(keys, classes):
    """({key: {class: taxid}} restricted to `keys`, {taxid: parent}, {taxid: rank}).

    names.dmp is 318 MB and nodes.dmp 233 MB uncompressed, so names are filtered to the keys
    actually being looked up while streaming rather than indexed in full.
    """
    idx = defaultdict(dict)
    parent, rank = {}, {}
    with zipfile.ZipFile(TAXDMP) as zf:
        for p in dmp_rows(zf, "names.dmp"):
            if len(p) < 4 or p[3] not in classes:
                continue
            k = norm(p[1]).lower()
            if k in keys:
                idx[k].setdefault(p[3], int(p[0]))
        for p in dmp_rows(zf, "nodes.dmp"):
            if len(p) < 3:
                continue
            t = int(p[0])
            parent[t], rank[t] = int(p[1]), p[2]
    return idx, parent, rank


def under(tid, root, parent):
    """Ancestor walk with a cycle guard: NCBI's root is its own parent."""
    seen = set()
    while tid and tid not in seen:
        if tid == root:
            return True
        seen.add(tid)
        nxt = parent.get(tid)
        if nxt == tid:
            break
        tid = nxt
    return False


def main():
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    rule = read_rule()
    classes = [c.strip() for c in rule["scope"]["name_classes"].split(",")] + ["includes"]

    crop_genera = [g.strip() for g in rule["scope"]["cohort_host_genera"].split(",")]
    aliases = dict(a.split("=", 1) for a in rule["scope"]["cohort_host_aliases"].split(","))
    aliases = {k.strip().lower(): v.strip() for k, v in aliases.items()}
    counts = usda_names(crop_genera, aliases)
    print(f"{len(counts):,} distinct USDA fungus names")

    keys = set()
    for n in counts:
        keys.add(n.lower())
        toks = n.split()
        if len(toks) > 2:
            keys.add(" ".join(toks[:2]).lower())
    idx, parent, rank = load_taxonomy(keys, set(classes))
    print(f"{len(idx):,} of {len(keys):,} lookup keys appear in names.dmp")

    name_of = {}
    with zipfile.ZipFile(TAXDMP) as zf:
        wanted = {t for d in idx.values() for t in d.values()}
        for p in dmp_rows(zf, "names.dmp"):
            if len(p) >= 4 and p[3] == "scientific name":
                t = int(p[0])
                if t in wanted:
                    name_of[t] = p[1]

    def match(n):
        """(taxid, status) for one USDA name, trying exact then binomial-trimmed."""
        for key, suffix in ((n.lower(), ""), ):
            for cls in classes:
                if cls in idx.get(key, {}):
                    return idx[key][cls], cls.replace(" ", "_") + suffix
        toks = n.split()
        if len(toks) > 2:
            key = " ".join(toks[:2]).lower()
            for cls in classes:
                if cls in idx.get(key, {}):
                    return idx[key][cls], cls.replace(" ", "_") + "_trimmed"
        return None, "no_ncbi_node"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    tally = defaultdict(int)
    resolved_taxids = set()
    on_crop, crop_taxids = defaultdict(int), defaultdict(set)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["usda_name", "n_records", "n_hosts", "ncbi_taxid", "matched_name",
                    "rank", "status", "on_wheat"])
        for n in sorted(counts):
            nrec, nhost, crops = counts[n]
            tid, status = match(n)
            if tid is not None and not under(tid, FUNGI, parent):
                tid, status = None, "not_fungi"
            if tid is not None:
                resolved_taxids.add(tid)
            tally[status] += 1
            w.writerow([n, nrec, nhost, tid or "", name_of.get(tid, "") if tid else "",
                        rank.get(tid, "") if tid else "", status,
                        "yes" if crops else "no"])
            for g in crops:
                on_crop[g] += 1
                if tid:
                    crop_taxids[g].add(tid)

    print("\nSTATUS")
    for k, v in sorted(tally.items(), key=lambda kv: -kv[1]):
        print(f"  {k:<28}{v:>8,}  ({100 * v / len(counts):.1f}%)")
    n_res = sum(v for k, v in tally.items() if k not in ("no_ncbi_node", "not_fungi"))
    print(f"\n  {n_res:,} resolved names -> {len(resolved_taxids):,} distinct fungal taxids")
    print(f"  {n_res - len(resolved_taxids):,} names are NCBI synonyms of a name already "
          f"in the list")
    print("\nRECORDED ON WHEAT (a label, not a filter)")
    for g in crop_genera:
        print(f"  {g:<16}{on_crop[g]:>7,} names   {len(crop_taxids[g]):>6,} distinct taxids")
    allc = set().union(*crop_taxids.values()) if crop_taxids else set()
    print(f"  {'any of them':<16}{'':>7}   {len(allc):>6,} distinct taxids "
          f"of {len(resolved_taxids):,}")

    print(f"\n  -> {OUT}")
    print("  (pathogen labels are a separate step: screen_pathogens.py)")


if __name__ == "__main__":
    main()
