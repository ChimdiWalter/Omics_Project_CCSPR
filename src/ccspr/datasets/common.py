from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

from ccspr.utils.io import ensure_dir


def download_first_available(urls: Iterable[str], out_path: str | Path, timeout: int = 120) -> str:
    out_path = Path(out_path)
    ensure_dir(out_path.parent)
    last_err = None
    for url in urls:
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200 and len(r.content) > 0:
                with open(out_path, "wb") as f:
                    f.write(r.content)
                return url
            last_err = RuntimeError(f"HTTP {r.status_code} for {url}")
        except Exception as e:
            last_err = e
    raise RuntimeError(f"Failed to download any URL to {out_path}: {last_err}")


def read_tsv_auto(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    with open(path, "rb") as f:
        head = f.read(2)
    if head == b"\x1f\x8b":
        with gzip.open(path, "rb") as f:
            raw = f.read()
        return pd.read_csv(io.BytesIO(raw), sep="\t", low_memory=False)
    # Some sources return plain text even when filename ends in .gz.
    return pd.read_csv(path, sep="\t", low_memory=False, comment="#", compression=None)


def normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).replace(".", "-").upper()
    if sid.startswith("TCGA-") and len(sid) >= 12:
        sid = sid[:12]
    return sid


def standardize_index(ids):
    return [normalize_sample_id(x) for x in ids]
