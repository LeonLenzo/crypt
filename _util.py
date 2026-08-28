#!/usr/bin/env python3
"""Shared utilities for the crypt pipeline."""

import hashlib
import json
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
