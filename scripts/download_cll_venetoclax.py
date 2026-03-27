#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ccspr.datasets.geo_utils import download_url, parse_geo_series_matrix, read_geo_table
from ccspr.utils.io import ensure_dir

GSE = "GSE161711"
COUNTS_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/{GSE}/suppl/{GSE}_log2CPM_counts.txt.gz"
SERIES_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/{GSE}/matrix/{GSE}_series_matrix.txt.gz"


def _dedupe(values: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for v in values:
        vv = str(v)
        if not vv:
            vv = "unknown_feature"
        n = seen.get(vv, 0)
        seen[vv] = n + 1
        out.append(vv if n == 0 else f"{vv}__dup{n}")
    return out


def _sample_title(sample_id: str) -> str:
    return re.sub(r"^S\d+_", "", sample_id)


def _infer_label(sample_id: str, tissue: str) -> str:
    sid = str(sample_id).upper()
    tissue_l = str(tissue).lower()
    if "_PB_" in sid:
        return "PB_CD19"
    if "_LN_CD19" in sid:
        return "LN_CD19"
    if "_LN_MIX" in sid:
        if "NLN-" in sid or "normal" in tissue_l:
            return "NLN_LN_mix"
        return "CLL_LN_mix"
    if "peripheral blood" in tissue_l:
        return "PB"
    if "lymph node" in tissue_l:
        return "LN"
    return "unknown"


def build_processed_matrix(root: Path, force: bool = False) -> dict[str, str]:
    raw_dir = ensure_dir(root / "raw")
    out_matrix = root / "processed_matrix.tsv"
    out_meta = root / "sample_metadata.tsv"
    out_readme = root / "README.txt"
    readme_text = (
        "CLL bulk public substitute dataset used by this repo:\n"
        "  GEO accession: GSE161711\n"
        "  Files downloaded:\n"
        "    raw/GSE161711_log2CPM_counts.txt.gz\n"
        "    raw/GSE161711_series_matrix.txt.gz\n"
        "  Generated files:\n"
        "    processed_matrix.tsv\n"
        "    sample_metadata.tsv\n"
        "\n"
        "Label mapping in processed_matrix.tsv is inferred from sample naming/metadata:\n"
        "  PB_CD19, LN_CD19, CLL_LN_mix, NLN_LN_mix\n"
    )

    counts_gz = raw_dir / f"{GSE}_log2CPM_counts.txt.gz"
    series_gz = raw_dir / f"{GSE}_series_matrix.txt.gz"
    download_url(COUNTS_URL, counts_gz)
    download_url(SERIES_URL, series_gz)

    if out_matrix.exists() and out_meta.exists() and not force:
        out_readme.write_text(readme_text)
        return {
            "processed_matrix_path": str(out_matrix),
            "metadata_path": str(out_meta),
            "readme_path": str(out_readme),
            "counts_path": str(counts_gz),
            "series_matrix_path": str(series_gz),
        }

    expr = read_geo_table(counts_gz, sep="\t")
    meta = parse_geo_series_matrix(series_gz)

    gene_id = expr["gene_id"].astype(str) if "gene_id" in expr.columns else pd.Series(np.arange(expr.shape[0], dtype=int))
    gene_name = expr["gene_name"].astype(str) if "gene_name" in expr.columns else gene_id.astype(str)
    feature_names = _dedupe((gene_name.where(gene_name.str.lower() != "nan", gene_id)).tolist())

    non_feature = {
        "gene_id",
        "gene_type",
        "gene_name",
        "chr",
        "source",
        "start",
        "stop",
        "strand",
        "gene_status",
    }
    sample_cols = [c for c in expr.columns if c not in non_feature]
    if not sample_cols:
        raise ValueError("No sample columns detected in GSE161711 log2CPM matrix.")

    x = expr[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).T

    meta = meta.copy()
    meta["sample_title_clean"] = meta["sample_title"].astype(str)
    meta = meta.drop_duplicates(subset=["sample_title_clean"])
    meta_idx = meta.set_index("sample_title_clean", drop=False)

    rows = []
    labels = []
    for sid in sample_cols:
        title = _sample_title(str(sid))
        tissue = ""
        geo = ""
        source_name = ""
        if title in meta_idx.index:
            m = meta_idx.loc[title]
            tissue = str(m.get("tissue", ""))
            geo = str(m.get("geo_accession", ""))
            source_name = str(m.get("source_name_ch1", ""))
        label = _infer_label(str(sid), tissue)
        labels.append(label)
        rows.append(
            {
                "sample_id": str(sid),
                "sample_title": title,
                "geo_accession": geo,
                "source_name": source_name,
                "tissue": tissue,
                "label": label,
            }
        )

    proc = pd.DataFrame(x, index=sample_cols, columns=feature_names).reset_index(drop=True)
    proc.insert(0, "label", labels)
    proc.insert(0, "sample_id", [str(s) for s in sample_cols])
    proc.to_csv(out_matrix, sep="\t", index=False)

    meta_df = pd.DataFrame(rows)
    meta_df.to_csv(out_meta, sep="\t", index=False)

    out_readme.write_text(readme_text)

    return {
        "processed_matrix_path": str(out_matrix),
        "metadata_path": str(out_meta),
        "readme_path": str(out_readme),
        "counts_path": str(counts_gz),
        "series_matrix_path": str(series_gz),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Download and convert GSE161711 into CLL bulk substitute matrix.")
    ap.add_argument("--root", default="data/cll_venetoclax")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = ensure_dir(Path(args.root))
    out = build_processed_matrix(root=root, force=bool(args.force))
    print(out)


if __name__ == "__main__":
    main()
