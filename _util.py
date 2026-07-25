#!/usr/bin/env python3
"""Shared utilities for the crypt pipeline."""

import json
import sys
import time
import urllib.error
import urllib.request
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


def http_get(url: str, headers: dict[str, str], retries: int = 5) -> bytes:
    """GET with exponential backoff on 429/5xx."""
    delay = 1.0
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(delay); delay *= 2
            else:
                raise
        except Exception:
            time.sleep(delay); delay *= 2
    raise RuntimeError(f"Failed after {retries} retries: {url}")


def load_json(path: Path) -> dict:
    """Load JSON from path; return {} if absent."""
    return json.loads(path.read_text()) if path.exists() else {}


def save_json(data: dict, path: Path) -> None:
    """Write JSON to path, creating parent dirs as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
