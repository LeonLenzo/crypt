#!/usr/bin/env python3
"""Shared utilities for the crypt pipeline."""

import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path


class _Tee:
    """Duplicate sys.stdout to a log file so all print() calls are captured."""
    def __init__(self, path: Path) -> None:
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


def make_log_dir(base: Path) -> Path:
    """Create a timestamped subdirectory under base/history/."""
    ts      = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_dir = base / "history" / ts
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir


def link_latest(base: Path, target: Path) -> None:
    """Create/update a .latest symlink at base/{name}.latest pointing to target.
    Log files (.log) keep their extension: find.log.latest
    Text files (.txt) drop it: find_summary.latest
    """
    import os
    stem = target.stem if target.suffix == ".txt" else target.name
    link = base / f"{stem}.latest"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(os.path.relpath(target, base))


def http_get(url: str, headers: dict[str, str], retries: int = 5,
             no_retry_429: bool = False) -> bytes:
    """GET with exponential backoff on 429/5xx.
    no_retry_429=True: raise immediately on 429 (for external APIs with own rate limiters).
    """
    delay = 1.0
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429:
                if no_retry_429:
                    raise
                time.sleep(delay); delay *= 2
            elif e.code >= 500:
                time.sleep(delay); delay *= 2
            else:
                raise
        except Exception:
            time.sleep(delay); delay *= 2
    raise RuntimeError(f"Failed after {retries} retries: {url}")


def load_json(path: Path) -> dict:
    """Load JSON from path; return {} if absent or empty."""
    if not path.exists():
        return {}
    text = path.read_text().strip()
    return json.loads(text) if text else {}


def save_json(data: dict, path: Path) -> None:
    """Write JSON to path atomically (temp file + rename) to avoid corruption on kill."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(data, f)
    tmp.rename(path)


def _md5(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def build_manifest(data_dir: Path, out_path: Path = None,
                   checksum_suffixes: tuple = (".k2d",)) -> Path:
    """Walk a large, gitignored data directory and write a small manifest.tsv
    (relative path, size in bytes, mtime, +md5 for files matching
    checksum_suffixes) — so the repo shows what exists on remote/HPC scratch
    without holding the data itself. Skips its own output file if re-run in
    place. Returns the manifest path.

    checksum_suffixes: file extensions worth checksumming (e.g. Kraken2 .k2d
    DB files, where byte-level reproducibility matters). Everything else is
    sized/timestamped only — checksumming hundreds of GB of FASTQ/CDS for no
    real benefit isn't worth the time.
    """
    out_path = out_path or (data_dir / "manifest.tsv")
    rows = []
    n_files = 0
    total_bytes = 0
    for p in sorted(data_dir.rglob("*")):
        if p == out_path or not p.is_file():
            continue
        rel = p.relative_to(data_dir)
        st = p.stat()
        checksum = _md5(p) if p.suffix in checksum_suffixes else ""
        rows.append((str(rel), st.st_size, int(st.st_mtime), checksum))
        n_files += 1
        total_bytes += st.st_size

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("path\tsize_bytes\tmtime_epoch\tmd5\n")
        for rel, size, mtime, checksum in rows:
            f.write(f"{rel}\t{size}\t{mtime}\t{checksum}\n")

    print(f"{out_path}: {n_files} files, {total_bytes / 1e9:.2f} GB")
    return out_path


# ── taxon name -> NCBI taxid resolution ───────────────────────────────────────
#
# Used wherever an LLM extracts a species name as free text (scientific name,
# common name, cultivar-suffixed, "common (scientific)" pairs) and something
# downstream needs an actual NCBI taxid — e.g. meta_classify.py resolving
# named_hosts/named_pathogens right after extraction (deterministic lookup,
# not asking the LLM to recall taxids from memory, which invites confident-
# looking wrong numbers), and kraken_run_select.py consuming those taxids.
#
# Common-name aliases observed across the full 2,719-sample field/aerial cohort
# (see kraken_restructure_plan.md / session notes 2026-08-28) plus generically
# common PHI-base-relevant crop names. Extend as new unresolved names turn up.
HOST_NAME_ALIASES = {
    "wheat": "Triticum aestivum", "durum wheat": "Triticum durum",
    "diploid wheat": "Triticum monococcum",
    "rice": "Oryza sativa", "maize": "Zea mays", "corn": "Zea mays",
    "sweet corn": "Zea mays", "barley": "Hordeum vulgare",
    "potato": "Solanum tuberosum", "tomato": "Solanum lycopersicum",
    "soybean": "Glycine max", "cotton": "Gossypium hirsutum",
    "sorghum": "Sorghum bicolor", "millet": "Setaria italica",
    "oat": "Avena sativa", "rye": "Secale cereale",
    "grape": "Vitis vinifera", "grapevine": "Vitis vinifera",
    "vine": "Vitis vinifera", "apple": "Malus domestica",
    "banana": "Musa acuminata", "coffee": "Coffea arabica",
    "cassava": "Manihot esculenta", "sweet potato": "Ipomoea batatas",
    "ipomea batatas": "Ipomoea batatas",   # observed misspelling (Ipomea -> Ipomoea)
    "lentil": "Lens culinaris", "pea": "Pisum sativum",
    "field pea": "Pisum sativum",
    "canola": "Brassica napus", "rapeseed": "Brassica napus",
    "oilseed rape": "Brassica napus",
    "sugarcane": "Saccharum officinarum", "sugar beet": "Beta vulgaris",
    "chard": "Beta vulgaris", "swiss chard": "Beta vulgaris",
    "strawberry": "Fragaria x ananassa", "triticale": "Triticosecale",
    "bean": "Phaseolus vulgaris", "common bean": "Phaseolus vulgaris",
    "broad bean": "Vicia faba", "faba bean": "Vicia faba",
    "cowpea": "Vigna unguiculata",
    "cucumber": "Cucumis sativus", "melon": "Cucumis melo",
    "watermelon": "Citrullus lanatus", "squash": "Cucurbita pepo",
    "pumpkin": "Cucurbita pepo",
    "lettuce": "Lactuca sativa",
    "pepper": "Capsicum annuum", "chili": "Capsicum annuum",
    "spinach": "Spinacia oleracea",
    "sunflower": "Helianthus annuus",
    "pear": "Pyrus communis",
    "nectarine": "Prunus persica",   # cultivated peach variety, no separate species taxid
    "switchgrass": "Panicum virgatum",
    "teosinte": "Zea mays",   # wild ancestor of maize; closest cultivated relative
    "barberis vulgaris": "Berberis vulgaris",   # observed misspelling (Barberis -> Berberis)
    "subterranean clover": "Trifolium subterraneum",
}

_TAXON_TRAILING_AUTHORITY = re.compile(
    r"\s+(l\.?|mill\.?|dc\.?|lam\.?|linn\.?|willd\.?|thunb\.?)$"
)
_TAXON_PAREN = re.compile(r"\(([^)]+)\)")


def _clean_taxon_name(s: str) -> str:
    s = s.strip()
    s = re.sub(r"\s*'[^']*'", "", s)          # strip 'cultivar' quotes
    s = _TAXON_TRAILING_AUTHORITY.sub("", s.lower()).strip()
    return s


def _taxon_candidates(part: str) -> list:
    """Generate resolution candidate strings from one name part: the cleaned
    string itself, anything in parentheses, anything outside parentheses, and
    progressively shorter word-prefixes (drops subspecies/variety trailing
    words, and — critically — still includes single-word candidates so bare
    common names like "wheat" reach the alias table)."""
    part = part.strip()
    cands = []
    for m in _TAXON_PAREN.finditer(part):
        cands.append(_clean_taxon_name(m.group(1)))
    outside = _TAXON_PAREN.sub("", part).strip()
    if outside:
        cands.append(_clean_taxon_name(outside))
    cands.append(_clean_taxon_name(part))
    expanded = []
    for c in cands:
        words = c.split()
        for n in range(len(words), 0, -1):
            expanded.append(" ".join(words[:n]))
    seen, uniq = set(), []
    for c in expanded:
        if c and c not in seen:
            seen.add(c)
            uniq.append(c)
    return uniq


def resolve_taxon_name(raw: str, name_to_taxid: dict, aliases: dict = None) -> tuple:
    """Resolve a free-text taxon name (host or pathogen, possibly ';'-separated
    multi-value, possibly a common name) to an NCBI taxid. Tries, in order:
    phibase_db.json's name_to_taxid (scientific names) -> live ete3 NCBITaxa
    (best-effort, optional dependency — skipped if ete3 isn't installed) —
    against BOTH the raw candidate strings AND their alias-table substitutions
    (e.g. HOST_NAME_ALIASES maps "spinach" -> "Spinacia oleracea", but that
    scientific name isn't itself in phibase_db either — it only resolves via
    ete3 — so the alias substitution needs the same two-tier treatment as a
    raw candidate, not just a single phibase_db lookup).
    Returns (taxid_or_None, resolved_name_or_None, method_str) where method_str
    is one of: phibase_db | alias_table | ete3_live | no_name_stated | unresolved."""
    aliases = aliases or {}
    if not raw or not raw.strip():
        return None, None, "no_name_stated"

    # Build the full candidate list once: raw name-parts plus their alias
    # substitutions, each tagged with where it came from (for the returned
    # method_str) so both get an equal shot at phibase_db AND ete3.
    candidates = []   # [(candidate_str, tag), ...]
    seen = set()
    for part in raw.split(";"):
        for cand in _taxon_candidates(part):
            if cand not in seen:
                seen.add(cand)
                candidates.append((cand, "phibase_db"))
            if cand in aliases:
                sci = aliases[cand].lower()
                if sci not in seen:
                    seen.add(sci)
                    candidates.append((sci, "alias_table"))

    for cand, tag in candidates:
        if cand in name_to_taxid:
            return name_to_taxid[cand], cand, tag

    try:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        for cand, _tag in candidates:
            trans = ncbi.get_name_translator([cand.capitalize()])
            if trans:
                tid = list(trans.values())[0][0]
                return tid, cand, "ete3_live"
    except Exception:
        pass
    return None, None, "unresolved"


def upload_to_acacia(local_dir: Path, s3_prefix: str, bucket: str,
                     endpoint: str = "https://projects.pawsey.org.au",
                     profile: str = "acacia") -> bool:
    """Sync a local directory to Pawsey's Acacia S3 object store via `aws s3 sync`.
    Returns True on success. Used for archiving large Setonix-scratch data
    (survives scratch wipes) — currently kraken/ DB build assemblies; any module
    with large regeneratable-but-expensive data on scratch can reuse this."""
    s3_uri = f"s3://{bucket}/{s3_prefix}/"
    print(f"\nUploading {local_dir} -> {s3_uri} …", flush=True)
    cmd = ["aws", "s3", "sync", str(local_dir), s3_uri,
           "--profile", profile, "--endpoint-url", endpoint]
    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        print(f"WARNING: upload to Acacia failed (exit {result.returncode})", flush=True)
        return False
    print(f"Upload complete -> {s3_uri}", flush=True)
    return True
