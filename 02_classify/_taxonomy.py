#!/usr/bin/env python3
"""_taxonomy.py — ONE taxonomy backbone, read once from holobase.

Every taxid in this module resolves here and nowhere else. That is a deliberate constraint,
not tidiness: a second backbone means two answers to "what is this taxon", and that is the
exact class of problem that has already cost this project real data.

  * **holobase's `taxon_id` is an internal surrogate, not an NCBI taxid.** The schema says so
    and the real key is `taxon.ncbi_taxid`. Joining db_v3's taxids against `taxon_id` returns a
    confident **0% overlap** rather than an error. Both spaces exist here and are never mixed:
    internal ids are called `tid`, NCBI taxids are called `ncbi`.
  * **Homonyms select nothing, silently.** The `datasets` CLI refuses an ambiguous NAME
    ("exact match for more than one taxid") and the clade then returns no assemblies with no
    sign anything went wrong. *Cryptococcus* is a fungus (5206) and a scale insect (79213);
    *Rhizophagus* a fungus (1129544) and a beetle (295919). Between them that cost 42
    assemblies and 17 species. **So resolve to a taxid here and query NCBI by taxid**, which
    is unambiguous by construction. That is what this module is for.
  * **Synonyms.** `taxon_name` holds 5.2 M NCBI names with a `name_class`, so
    *Magnaporthe oryzae* resolves to the same node as *Pyricularia oryzae* (internal 263712,
    ncbi 318829). A name-keyed membership list without this silently holds both as two taxa.
  * **Rank.** A seed list may name a forma specialis where the database holds the species node,
    which breaks a depth rule that counts per taxid. `species_of()` is the fix and callers
    should apply it before any comparison or count.
  * **`kingdom` spells plants `Viridiplantae`, not `Plants`.** Filtering on 'Plants' silently
    returns zero.
  * **An ancestor walk needs a cycle guard**, because NCBI's root is its own parent. An
    unbounded walk never terminates for a taxon with no ancestor at the requested rank.

`taxon_merge` is empty (0 rows), so there is no merge decision layer to honour yet.

One pass over the ~3 M-row `taxon` table fills everything, including the species and phylum
maps `_assembly.py` needs, which used to cost a full pass each.
"""

import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _layout import HOLOBASE  # noqa: E402

# Phylum -> BUSCO lineage, kept here with the taxonomy that supplies the phylum.
PHYLUM_LINEAGE = {
    "Basidiomycota": "basidiomycota_odb10",
    "Ascomycota": "ascomycota_odb10",
    "Chytridiomycota": "fungi_odb10",
    "Blastocladiomycota": "fungi_odb10",
    "Mucoromycota": "fungi_odb10",
    "Zoopagomycota": "fungi_odb10",
    "Oomycota": "stramenopiles_odb10",
}

_CACHE = None


class Taxonomy:
    """holobase's taxonomy, in both id spaces, with the walks the selectors need."""

    def __init__(self, db_path: Path = HOLOBASE):
        self.parent, self.rank, self.ncbi, self.name, self.kingdom = {}, {}, {}, {}, {}
        self.children = defaultdict(list)
        self.by_ncbi = {}
        con = sqlite3.connect(db_path)
        for tid, par, rk, nc, nm, kgd in con.execute(
                "SELECT taxon_id,parent_id,rank,ncbi_taxid,scientific_name,kingdom FROM taxon"):
            self.parent[tid], self.rank[tid], self.name[tid] = par, rk, nm
            self.kingdom[tid] = kgd
            if nc is not None:
                self.ncbi[tid] = nc
                self.by_ncbi[nc] = tid
            if par is not None and par != tid:
                self.children[par].append(tid)
        # Names, scientific first so a synonym can never shadow an accepted name.
        self.by_name = {}
        for tid, nm, cls in con.execute(
                "SELECT taxon_id,name,name_class FROM taxon_name"):
            key = nm.strip().lower()
            if cls == "scientific name" or key not in self.by_name:
                self.by_name.setdefault(key, tid)
                if cls == "scientific name":
                    self.by_name[key] = tid
        con.close()
        self._species_cache = {}

    # ── resolution ────────────────────────────────────────────────────────────

    def resolve(self, token) -> int:
        """A name, an NCBI taxid, or a taxid as a string -> internal id, or 0.

        Accepts a taxid because a membership list MUST be allowed to carry one: it is the only
        way to name a homonym unambiguously, which is why column 1 of `residents.tsv` may be a
        taxid. A bare integer is read as an NCBI taxid, never as an internal id, because
        internal ids are an implementation detail of holobase and never appear in a seeds file.
        """
        if token is None:
            return 0
        s = str(token).strip()
        if not s:
            return 0
        if s.isdigit():
            return self.by_ncbi.get(int(s), 0)
        return self.by_name.get(s.lower(), 0)

    def ncbi_of(self, tid: int) -> int:
        return self.ncbi.get(tid, 0)

    def label(self, tid: int) -> str:
        return self.name.get(tid, str(tid))

    # ── walks, all cycle-guarded ──────────────────────────────────────────────

    def ancestor_at(self, tid: int, rank: str) -> int:
        """Nearest ancestor-or-self at `rank`, or 0. Guards the self-parent root."""
        seen = set()
        while tid in self.rank and tid not in seen:
            seen.add(tid)
            if self.rank[tid] == rank:
                return tid
            nxt = self.parent.get(tid)
            if nxt is None or nxt == tid or nxt in seen:
                break
            tid = nxt
        return 0

    def species_of(self, tid: int) -> int:
        """Nearest ancestor-or-self at species rank, memoised over shared lineages."""
        if tid in self._species_cache:
            return self._species_cache[tid]
        chain, found, cur = [], 0, tid
        seen = set()
        while cur in self.rank and cur not in seen:
            seen.add(cur)
            if cur in self._species_cache:
                found = self._species_cache[cur]
                break
            if self.rank[cur] == "species":
                found = cur
                break
            chain.append(cur)
            nxt = self.parent.get(cur)
            if nxt is None or nxt == cur or nxt in seen:
                break
            cur = nxt
        for t in chain:
            self._species_cache[t] = found
        self._species_cache[tid] = found
        return found

    def is_ancestor(self, anc: int, tid: int) -> bool:
        seen = set()
        while tid in self.parent and tid not in seen:
            if tid == anc:
                return True
            seen.add(tid)
            nxt = self.parent.get(tid)
            if nxt is None or nxt == tid or nxt in seen:
                break
            tid = nxt
        return tid == anc

    def descendants(self, root: int):
        """Every node in the subtree at `root`, including root."""
        stack, seen = [root], set()
        while stack:
            t = stack.pop()
            if t in seen:
                continue
            seen.add(t)
            yield t
            stack.extend(self.children.get(t, ()))

    # ── the maps _assembly.py needs, from the same single pass ─────────────────

    def species_map(self) -> dict:
        """{ncbi taxid: species-level ncbi taxid}. Species map to themselves, so a caller may
        apply `to_species.get(tid, tid)` unconditionally."""
        out = {}
        for tid, nc in self.ncbi.items():
            sp = self.species_of(tid)
            if sp and sp in self.ncbi:
                out[nc] = self.ncbi[sp]
        return out

    def phylum_map(self) -> dict:
        """{ncbi taxid: phylum name}, for the BUSCO lineage call."""
        out = {}
        for tid, nc in self.ncbi.items():
            p = self.ancestor_at(tid, "phylum")
            if p:
                out[nc] = self.name.get(p)
        return out


def load(db_path: Path = HOLOBASE) -> Taxonomy:
    """The process-wide Taxonomy. Reading it twice is ~3 M rows twice, so it is cached."""
    global _CACHE
    if _CACHE is None:
        if not db_path.exists():
            raise FileNotFoundError(
                f"{db_path} not found. There is no fallback on purpose: every membership "
                "decision in this module resolves through holobase, and a silent fallback to "
                "name matching is how homonyms and synonyms get into the list.")
        _CACHE = Taxonomy(db_path)
    return _CACHE
