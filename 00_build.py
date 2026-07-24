#!/usr/bin/env python3
"""
00_build.py — build PHI-base + ICTV reference database with full NCBI taxonomy.

Must be run with system python3 (not miniconda) due to ete3/sqlite3
binary incompatibility:
    python3 00_build.py

Reads PHI-base CSV and ICTV VMR, resolves all pathogen/virus/host taxids
through the local NCBI taxonomy (ete3 NCBITaxa), expands each to include
all descendant taxa (strains, subspecies), and collects all known name
variants (scientific names + synonyms) for STAT name resolution.

PHI-base supplies fungi, bacteria, oomycetes, nematodes.
ICTV VMR supplies plant viruses (Host source: plants / invertebrates, plants).

Output: output/00_build/phibase_db.json

The database is consumed by all downstream pipeline steps.
"""

import argparse
import csv
import sys
import json
import shutil
import tempfile
import urllib.request
from datetime import date
from pathlib import Path

from ete3 import NCBITaxa

VIRIDIPLANTAE = 33090   # scope: plant hosts only

PHIBASE_URL = (
    "https://raw.githubusercontent.com/PHI-base/data/master/releases/"
    "phi-base_current.csv"
)
ICTV_VMR_URL = "https://ictv.global/vmr/current"

# Host source values in ICTV VMR that indicate plant viruses.
# "invertebrates, plants" = insect-vectored viruses (geminiviruses, tospoviruses etc.)
ICTV_PLANT_HOST_SOURCES = {"plants", "plants (s)", "invertebrates, plants"}

CONTAMINANT_TAXIDS = {
    10847,   # Enterobacteria phage phiX174 — Illumina spike-in control
    9606,    # Homo sapiens — common contamination
    83333,   # Escherichia coli K-12 — common lab strain
    1408252, # Cellulophaga phage phi38:1 — library prep contaminant
}


def _is_authority_string(name: str) -> bool:
    """Exclude author+year strings like 'Fusarium oxysporum Schltdl., 1824'."""
    return "," in name and any(c.isdigit() for c in name)


# ── PHI-base fetch ────────────────────────────────────────────────────────────

def fetch_phibase(dest: Path, url: str = PHIBASE_URL) -> None:
    """Download the latest PHI-base CSV from GitHub to dest."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Fetching PHI-base from {url} …", flush=True)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        tmp_path = Path(tmp.name)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "crypt/00_build"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            size = 0
            with open(tmp_path, "wb") as out:
                while chunk := resp.read(65536):
                    out.write(chunk)
                    size += len(chunk)

        with open(tmp_path, encoding="utf-8-sig", errors="replace") as f:
            header = f.readline()
        if "Pathogen_species" not in header and "Record ID" not in header:
            raise ValueError(f"Downloaded file does not look like PHI-base CSV "
                             f"(header: {header[:120]!r})")

        shutil.move(str(tmp_path), dest)
        print(f"  Saved {size / 1e6:.1f} MB → {dest}", flush=True)

    except Exception as e:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"PHI-base download failed: {e}") from e


# ── ICTV VMR fetch ────────────────────────────────────────────────────────────

def fetch_ictv_vmr(dest: Path, url: str = ICTV_VMR_URL) -> None:
    """Download the latest ICTV VMR Excel file. URL redirects to current release."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Fetching ICTV VMR from {url} …", flush=True)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
        tmp_path = Path(tmp.name)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "crypt/00_build"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            final_url = resp.url
            size = 0
            with open(tmp_path, "wb") as out:
                while chunk := resp.read(65536):
                    out.write(chunk)
                    size += len(chunk)

        if size < 10000:
            raise ValueError(f"Downloaded file too small ({size} bytes) — likely an error page")

        shutil.move(str(tmp_path), dest)
        print(f"  Saved {size / 1e6:.1f} MB → {dest}  (from {final_url})", flush=True)

    except Exception as e:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"ICTV VMR download failed: {e}") from e


# ── ICTV VMR parsing ──────────────────────────────────────────────────────────

def _parse_ictv_plant_viruses(vmr_path: Path) -> list[str]:
    """Return unique plant virus species names from the ICTV VMR Excel file."""
    try:
        import openpyxl
    except ImportError:
        raise ImportError(
            "openpyxl is required to parse the ICTV VMR.\n"
            "Install with:  pip install openpyxl"
        )

    wb = openpyxl.load_workbook(vmr_path, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))

    hdrs = {str(h).strip(): i for i, h in enumerate(rows[0]) if h}
    host_col    = hdrs.get("Host source")
    species_col = hdrs.get("Species")

    if host_col is None or species_col is None:
        raise ValueError(f"ICTV VMR missing expected columns. Found: {list(hdrs)}")

    seen: set[str] = set()
    species: list[str] = []
    for row in rows[1:]:
        host = str(row[host_col] or "").lower().strip()
        if host not in ICTV_PLANT_HOST_SOURCES:
            continue
        sp = str(row[species_col] or "").strip()
        if sp and sp not in seen:
            seen.add(sp)
            species.append(sp)

    wb.close()
    return species


def _load_ictv_viruses(vmr_path: Path,
                       ncbi: NCBITaxa) -> tuple[set[int], dict[str, int]]:
    """
    Parse ICTV VMR, filter to plant viruses, resolve species names to NCBI taxids.
    Returns:
      virus_seed_taxids: set of resolved NCBI taxids for plant virus species
      ictv_name_map:     {lowercase_ictv_name: taxid} for name_to_taxid supplement
    """
    species_names = _parse_ictv_plant_viruses(vmr_path)
    print(f"  ICTV plant virus species (all clades): {len(species_names):,}", flush=True)

    name_map = ncbi.get_name_translator(species_names)
    resolved   = {name: taxids[0] for name, taxids in name_map.items() if taxids}
    unresolved = len(species_names) - len(resolved)
    print(f"  Resolved to NCBI taxids: {len(resolved):,}  "
          f"(unresolved / not in NCBI: {unresolved:,})", flush=True)

    virus_seed_taxids = set(resolved.values())
    ictv_name_map     = {name.lower(): taxid for name, taxid in resolved.items()}
    return virus_seed_taxids, ictv_name_map


# ── PHI-base parsing ──────────────────────────────────────────────────────────

def _parse_phibase(csv_path: Path) -> tuple[dict[int, set[int]], set[int], set[int]]:
    """
    Parse PHI-base CSV.
    Returns pathogen_to_hosts, all_pathogen_taxids, all_host_taxids.
    """
    pathogen_to_hosts: dict[int, set[int]] = {}
    bad_rows = 0

    with open(csv_path, encoding="utf-8-sig", errors="replace") as f:
        reader = csv.DictReader(f)
        next(reader)  # skip second human-readable header row
        for row in reader:
            try:
                p_taxid = int(row.get("Pathogen_NCBI_species_Taxonomy ID", "").strip())
                h_taxid = int(row.get("Host_NCBI_Taxonomy_ID", "").strip())
            except ValueError:
                bad_rows += 1
                continue
            if p_taxid <= 0 or h_taxid <= 0:
                bad_rows += 1
                continue
            pathogen_to_hosts.setdefault(p_taxid, set()).add(h_taxid)

    all_pathogen_taxids = set(pathogen_to_hosts.keys())
    all_host_taxids     = {h for hosts in pathogen_to_hosts.values() for h in hosts}

    print(f"  PHI-base: {sum(len(v) for v in pathogen_to_hosts.values()):,} "
          f"interactions, {len(all_pathogen_taxids)} pathogen taxids, "
          f"{len(all_host_taxids)} host taxids "
          f"({bad_rows} rows skipped — missing/non-numeric taxid)")
    return pathogen_to_hosts, all_pathogen_taxids, all_host_taxids


def _filter_plant_hosts(host_taxids: set[int],
                        pathogen_to_hosts: dict[int, set[int]],
                        ncbi: NCBITaxa) -> tuple[set[int], dict[int, set[int]]]:
    """Restrict to plant-host entries (Viridiplantae = 33090)."""
    plant_hosts: set[int] = set()
    for taxid in host_taxids:
        try:
            lineage = ncbi.get_lineage(taxid)
            if lineage and VIRIDIPLANTAE in lineage:
                plant_hosts.add(taxid)
        except Exception:
            pass

    filtered: dict[int, set[int]] = {}
    for p_taxid, h_taxids in pathogen_to_hosts.items():
        plant_h = h_taxids & plant_hosts
        if plant_h:
            filtered[p_taxid] = plant_h

    plant_pathogens = set(filtered.keys())
    print(f"  After plant filter: {len(plant_pathogens)} pathogens, "
          f"{len(plant_hosts)} hosts")
    return plant_hosts, filtered


# ── Taxonomy expansion ────────────────────────────────────────────────────────

def _expand_taxids(seed_taxids: set[int], ncbi: NCBITaxa,
                   label: str) -> tuple[set[int], dict[int, int]]:
    """
    Expand seed taxids to include all descendants.
    Returns (expanded_set, {expanded_taxid → seed_taxid}).
    """
    expanded: set[int]        = set()
    taxid_to_seed: dict[int, int] = {}
    skipped = 0

    for seed in seed_taxids:
        try:
            descendants = ncbi.get_descendant_taxa(seed, collapse_subspecies=False)
            expanded.add(seed)
            taxid_to_seed[seed] = seed
            for d in descendants:
                expanded.add(d)
                taxid_to_seed[d] = seed
        except Exception:
            skipped += 1

    print(f"  {label}: {len(seed_taxids)} seed → "
          f"{len(expanded):,} total (incl. descendants)"
          + (f"; {skipped} taxids not found in local taxonomy" if skipped else ""))
    return expanded, taxid_to_seed


def _build_name_map(taxids: set[int],
                    ncbi: NCBITaxa) -> tuple[dict[int, str], dict[str, int]]:
    """
    Build taxid→canonical_name and lowercase name→taxid (incl. synonyms).
    Authority strings ('Author, Year') are excluded.
    """
    taxid_to_name: dict[int, str] = {}
    name_to_taxid: dict[str, int] = {}

    if not taxids:
        return taxid_to_name, name_to_taxid

    id_list = ",".join(str(t) for t in taxids)

    rows = ncbi.db.execute(
        f"SELECT taxid, spname FROM species WHERE taxid IN ({id_list})"
    )
    for taxid, name in rows:
        taxid_to_name[taxid] = name
        if not _is_authority_string(name):
            name_to_taxid[name.lower()] = taxid

    rows = ncbi.db.execute(
        f"SELECT taxid, spname FROM synonym WHERE taxid IN ({id_list})"
    )
    for taxid, name in rows:
        if not _is_authority_string(name):
            name_to_taxid.setdefault(name.lower(), taxid)

    return taxid_to_name, name_to_taxid


# ── Main build ────────────────────────────────────────────────────────────────

def build(phibase_csv: Path, vmr_path: Path | None = None,
          host_scope: str = "plant") -> dict:
    """Build the full reference database. Returns a dict ready for JSON serialisation."""
    ncbi = NCBITaxa()
    print(f"NCBITaxa loaded: {ncbi.dbfile}")

    print("\n[1/6] Parsing PHI-base …")
    pathogen_to_hosts, seed_pathogens, seed_hosts = _parse_phibase(phibase_csv)

    if host_scope == "plant":
        print("\n[2/6] Filtering to plant hosts …")
        seed_hosts, pathogen_to_hosts = _filter_plant_hosts(
            seed_hosts, pathogen_to_hosts, ncbi
        )
        seed_pathogens = set(pathogen_to_hosts.keys())
    else:
        print("\n[2/6] No host scope filter applied.")

    print("\n[3/6] Loading ICTV plant viruses …")
    if vmr_path and vmr_path.exists():
        virus_seed_taxids, ictv_name_map = _load_ictv_viruses(vmr_path, ncbi)
    else:
        virus_seed_taxids, ictv_name_map = set(), {}
        print("  No ICTV VMR provided — virus step skipped")

    print("\n[4/6] Expanding pathogen taxonomy …")
    all_pathogen_taxids, pathogen_to_seed = _expand_taxids(seed_pathogens, ncbi, "pathogens")

    if virus_seed_taxids:
        all_virus_taxids, virus_to_seed = _expand_taxids(virus_seed_taxids, ncbi, "viruses")
        all_pathogen_taxids |= all_virus_taxids
        pathogen_to_seed.update(virus_to_seed)
    else:
        all_virus_taxids, virus_to_seed = set(), {}

    print("\n[5/6] Expanding host taxonomy …")
    all_host_taxids, host_to_seed = _expand_taxids(seed_hosts, ncbi, "hosts")

    print("\n[6/6] Building name lookup tables …")
    all_taxids = all_pathogen_taxids | all_host_taxids
    taxid_to_name, name_to_taxid = _build_name_map(all_taxids, ncbi)
    # Supplement with ICTV species names (different naming convention to NCBI)
    for name, taxid in ictv_name_map.items():
        name_to_taxid.setdefault(name, taxid)
    print(f"  name→taxid entries: {len(name_to_taxid):,}")

    # Reverse PHI-base host/pathogen map
    host_to_pathogens: dict[int, list[int]] = {}
    for p_taxid, h_taxids in pathogen_to_hosts.items():
        for h_taxid in h_taxids:
            host_to_pathogens.setdefault(h_taxid, [])
            if p_taxid not in host_to_pathogens[h_taxid]:
                host_to_pathogens[h_taxid].append(p_taxid)

    # Viruses: known to infect all PHI-base plant host seeds
    if virus_seed_taxids:
        for v_seed in virus_seed_taxids:
            pathogen_to_hosts[v_seed] = set(seed_hosts)
        for h_seed in seed_hosts:
            host_to_pathogens.setdefault(h_seed, [])
            for v_seed in virus_seed_taxids:
                if v_seed not in host_to_pathogens[h_seed]:
                    host_to_pathogens[h_seed].append(v_seed)

    db = {
        "meta": {
            "built":                    str(date.today()),
            "host_scope":               host_scope,
            "phibase_path":             str(phibase_csv),
            "vmr_path":                 str(vmr_path) if vmr_path else None,
            "n_pathogen_species":       len(seed_pathogens),
            "n_virus_species":          len(virus_seed_taxids),
            "n_host_species":           len(seed_hosts),
            "n_pathogen_taxids_total":  len(all_pathogen_taxids),
            "n_virus_taxids_total":     len(all_virus_taxids),
            "n_host_taxids_total":      len(all_host_taxids),
            "n_name_entries":           len(name_to_taxid),
        },
        "pathogen_taxids":    sorted(all_pathogen_taxids),
        "virus_taxids":       sorted(all_virus_taxids),
        "host_taxids":        sorted(all_host_taxids),
        "contaminant_taxids": sorted(CONTAMINANT_TAXIDS),
        "taxid_to_name":  {str(k): v for k, v in taxid_to_name.items()},
        "name_to_taxid":  name_to_taxid,
        "pathogen_taxid_to_seed": {str(k): v for k, v in pathogen_to_seed.items()},
        "host_taxid_to_seed":     {str(k): v for k, v in host_to_seed.items()},
        "pathogen_to_hosts": {str(k): sorted(v)
                              for k, v in pathogen_to_hosts.items()},
        "host_to_pathogens": {str(k): sorted(v)
                              for k, v in host_to_pathogens.items()},
    }
    return db


# ── Save / load ───────────────────────────────────────────────────────────────

def save_db(db: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(db, f, indent=2)
    size_mb = path.stat().st_size / 1e6
    print(f"\nSaved → {path}  ({size_mb:.1f} MB)")


def load_db(path: Path) -> "PhibaseDB":
    with open(path) as f:
        raw = json.load(f)
    return PhibaseDB(raw)


# ── Convenience wrapper ───────────────────────────────────────────────────────

class PhibaseDB:
    """Thin wrapper around the raw JSON database dict. Provides O(1) lookups."""

    def __init__(self, raw: dict):
        self._meta              = raw["meta"]
        self._pathogen_taxids   = set(raw["pathogen_taxids"])
        self._virus_taxids      = set(raw.get("virus_taxids", []))
        self._host_taxids       = set(raw["host_taxids"])
        self._contaminants      = set(raw["contaminant_taxids"])
        self._taxid_to_name     = {int(k): v for k, v in raw["taxid_to_name"].items()}
        self._name_to_taxid     = raw["name_to_taxid"]
        self._pathogen_to_seed  = {int(k): v
                                   for k, v in raw["pathogen_taxid_to_seed"].items()}
        self._host_to_seed      = {int(k): v
                                   for k, v in raw["host_taxid_to_seed"].items()}
        self._pathogen_to_hosts = {int(k): set(v)
                                   for k, v in raw["pathogen_to_hosts"].items()}
        self._host_to_pathogens = {int(k): set(v)
                                   for k, v in raw["host_to_pathogens"].items()}

    def resolve(self, name: str) -> int | None:
        return self._name_to_taxid.get(name.lower())

    def is_pathogen(self, taxid: int) -> bool:
        return taxid in self._pathogen_taxids

    def is_virus(self, taxid: int) -> bool:
        return taxid in self._virus_taxids

    def is_host(self, taxid: int) -> bool:
        return taxid in self._host_taxids

    def is_contaminant(self, taxid: int) -> bool:
        return taxid in self._contaminants

    def name(self, taxid: int) -> str:
        return self._taxid_to_name.get(taxid, str(taxid))

    def seed_pathogen(self, taxid: int) -> int | None:
        return self._pathogen_to_seed.get(taxid)

    def seed_host(self, taxid: int) -> int | None:
        return self._host_to_seed.get(taxid)

    def known_interaction(self, pathogen_taxid: int, host_taxid: int) -> bool:
        """
        True if PHI-base / ICTV records this pathogen as infecting this host.
        For plant viruses, returns True for any plant host (all PHI-base host seeds).
        """
        p_seed = self._pathogen_to_seed.get(pathogen_taxid, pathogen_taxid)
        h_seed = self._host_to_seed.get(host_taxid, host_taxid)
        return h_seed in self._pathogen_to_hosts.get(p_seed, set())

    def hosts_of(self, pathogen_taxid: int) -> set[int]:
        seed = self._pathogen_to_seed.get(pathogen_taxid, pathogen_taxid)
        return self._pathogen_to_hosts.get(seed, set())

    def pathogens_of(self, host_taxid: int) -> set[int]:
        seed = self._host_to_seed.get(host_taxid, host_taxid)
        return self._host_to_pathogens.get(seed, set())

    @property
    def meta(self) -> dict:
        return self._meta

    def __repr__(self) -> str:
        m = self._meta
        v = m.get("n_virus_taxids_total", 0)
        return (f"PhibaseDB(built={m['built']}, "
                f"{m['n_pathogen_taxids_total']:,} pathogen taxids "
                f"[incl. {v:,} virus], "
                f"{m['n_host_taxids_total']:,} host taxids, "
                f"{m['n_name_entries']:,} names)")


# ── CLI ───────────────────────────────────────────────────────────────────────

OUT_DIR = Path("output/00_build")


# ── Logging ───────────────────────────────────────────────────────────────────

class _Tee:
    """Duplicate sys.stdout to a log file so all print() calls are captured."""
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._log    = open(path, "w", buffering=1)
        self._stdout = sys.__stdout__

    def write(self, text: str) -> None:
        self._stdout.write(text)
        self._log.write(text)

    def flush(self) -> None:
        self._stdout.flush()
        self._log.flush()

    def close(self) -> None:
        sys.stdout = self._stdout
        self._log.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phibase", default=str(OUT_DIR / "phi-base_current.csv"),
                    help="Path to PHI-base CSV")
    ap.add_argument("--ictv", default=str(OUT_DIR / "ictv_vmr.xlsx"),
                    help="Path to ICTV VMR Excel file")
    ap.add_argument("--out", default=str(OUT_DIR / "phibase_db.json"),
                    help="Output JSON path")
    ap.add_argument("--scope", default="plant", choices=["plant", "all"],
                    help="Host scope filter (default: plant)")
    ap.add_argument("--fetch", action="store_true",
                    help="Force re-download of PHI-base and ICTV VMR")
    ap.add_argument("--phibase-url", default=PHIBASE_URL,
                    help="PHI-base download URL")
    ap.add_argument("--ictv-url", default=ICTV_VMR_URL,
                    help="ICTV VMR download URL")
    args = ap.parse_args()

    log = _Tee(OUT_DIR / "build.log")
    sys.stdout = log

    try:
        phibase_path = Path(args.phibase)
        if args.fetch or not phibase_path.exists():
            fetch_phibase(phibase_path, url=args.phibase_url)

        vmr_path = Path(args.ictv)
        if args.fetch or not vmr_path.exists():
            fetch_ictv_vmr(vmr_path, url=args.ictv_url)

        print(f"Building PHI-base + ICTV database from {phibase_path} "
              f"(scope={args.scope})\n")
        db = build(phibase_path, vmr_path=vmr_path, host_scope=args.scope)

        out_path = Path(args.out)
        save_db(db, out_path)

        loaded = load_db(out_path)
        print(f"\nSmoke test: {loaded}")
        test_name = "Fusarium oxysporum"
        taxid = loaded.resolve(test_name)
        print(f"  resolve('{test_name}') → {taxid} ({loaded.name(taxid) if taxid else 'not found'})")
        if taxid:
            print(f"  is_pathogen({taxid}) → {loaded.is_pathogen(taxid)}")
            print(f"  is_virus({taxid})    → {loaded.is_virus(taxid)}")
            print(f"  hosts_of({taxid})    → {len(loaded.hosts_of(taxid))} plant hosts")
            tomato = 4081
            print(f"  known_interaction(F.oxysporum, tomato) → "
                  f"{loaded.known_interaction(taxid, tomato)}")
        tmv_name = "Tobacco mosaic virus"
        tmv_taxid = loaded.resolve(tmv_name)
        print(f"  resolve('{tmv_name}') → {tmv_taxid} "
              f"({'virus ✓' if tmv_taxid and loaded.is_virus(tmv_taxid) else 'not found'})")

        m = db["meta"]
        summary = (
            f"── 00_build summary ────────────────────────────\n"
            f"Built:                    {m['built']}\n"
            f"Host scope:               {m['host_scope']}\n"
            f"PHI-base source:          {phibase_path}\n"
            f"ICTV VMR source:          {vmr_path}\n"
            f"\n"
            f"Pathogen species (seeds):  {m['n_pathogen_species']:>7,}\n"
            f"Pathogen taxids (total):   {m['n_pathogen_taxids_total']:>7,}\n"
            f"  of which virus taxids:   {m['n_virus_taxids_total']:>7,}\n"
            f"Host species (seeds):      {m['n_host_species']:>7,}\n"
            f"Host taxids (total):       {m['n_host_taxids_total']:>7,}\n"
            f"Name lookup entries:       {m['n_name_entries']:>7,}\n"
            f"\n"
            f"DB:  {out_path}\n"
            f"Log: {OUT_DIR / 'build.log'}\n"
        )
        (OUT_DIR / "_summary.txt").write_text(summary)
        print(f"\n{summary}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
