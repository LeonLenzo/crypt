#!/usr/bin/env python3
"""Shared utilities for the crypt pipeline."""

import json
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
