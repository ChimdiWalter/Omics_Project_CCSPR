from __future__ import annotations

import csv
import gzip
import io
import re
from pathlib import Path

import pandas as pd
import requests

from ccspr.utils.io import ensure_dir


def download_url(url: str, out_path: str | Path, timeout: int = 300, chunk_size: int = 1024 * 1024) -> Path:
    out_path = Path(out_path)
    ensure_dir(out_path.parent)
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    return out_path


def _to_snake(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "unknown"


def _parse_characteristic(values: list[str]) -> tuple[str, list[str]]:
    key = None
    parsed: list[str] = []
    for v in values:
        v = str(v).strip().strip('"')
        if ":" in v:
            k, val = v.split(":", 1)
            if key is None and k.strip():
                key = _to_snake(k)
            parsed.append(val.strip())
        else:
            parsed.append(v)
    if key is None:
        key = "characteristic"
    return key, parsed


def parse_geo_series_matrix(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip("\n") for ln in f]

    meta_rows: dict[str, list[str]] = {}
    char_blocks: list[list[str]] = []

    for ln in lines:
        if not ln.startswith("!Sample_"):
            continue
        parts = next(csv.reader([ln], delimiter="\t", quotechar='"'))
        field = parts[0][len("!Sample_") :]
        values = [p.strip().strip('"') for p in parts[1:]]
        if field.startswith("characteristics_ch1"):
            char_blocks.append(values)
        else:
            meta_rows[field] = values

    if "title" not in meta_rows:
        raise ValueError(f"Could not parse sample titles from series matrix: {path}")

    n = len(meta_rows["title"])
    df = pd.DataFrame({"sample_title": meta_rows["title"]})
    for field in ["geo_accession", "source_name_ch1", "description"]:
        vals = meta_rows.get(field)
        if vals and len(vals) == n:
            df[field] = vals

    seen: dict[str, int] = {}
    for block in char_blocks:
        if len(block) != n:
            continue
        key, parsed_vals = _parse_characteristic(block)
        idx = seen.get(key, 0)
        seen[key] = idx + 1
        col = key if idx == 0 else f"{key}_{idx + 1}"
        df[col] = parsed_vals

    return df


def read_geo_table(path: str | Path, sep: str = "\t") -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            raw = f.read()
        return pd.read_csv(io.BytesIO(raw), sep=sep, low_memory=False)
    return pd.read_csv(path, sep=sep, low_memory=False)
