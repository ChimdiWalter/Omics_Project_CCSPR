#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

from ccspr.datasets.geo_utils import download_url, parse_geo_series_matrix
from ccspr.utils.io import ensure_dir

GSE = "GSE165087"
RAW_TAR_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/{GSE}/suppl/{GSE}_RAW.tar"
SERIES_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/{GSE}/matrix/{GSE}_series_matrix.txt.gz"
FILELIST_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/{GSE}/suppl/filelist.txt"


def _dedupe(values: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for v in values:
        vv = str(v).strip()
        if not vv:
            vv = "unknown_gene"
        n = seen.get(vv, 0)
        seen[vv] = n + 1
        out.append(vv if n == 0 else f"{vv}__dup{n}")
    return out


def _sample_time_group(time_value: str) -> str:
    t = str(time_value).strip().lower()
    if "pretreatment" in t:
        return "pretreatment"
    if "watchful" in t:
        return "watchful_waiting"
    if "relapse" in t:
        return "relapse"
    if "ibrutinib" in t:
        return "on_ibrutinib"
    return "other"


def _stratified_sample_obs(obs: pd.DataFrame, key: str, n_max: int, seed: int) -> np.ndarray:
    if len(obs) <= int(n_max):
        return np.arange(len(obs), dtype=int)
    rng = np.random.default_rng(int(seed))
    keep = []
    counts = obs[key].astype(str).value_counts()
    for c, n in counts.items():
        idx = np.where(obs[key].astype(str).to_numpy() == c)[0]
        take = max(100, int(round((n / len(obs)) * int(n_max))))
        take = min(take, len(idx))
        keep.extend(rng.choice(idx, size=take, replace=False).tolist())
    keep = np.array(sorted(set(keep)), dtype=int)
    if len(keep) > int(n_max):
        keep = rng.choice(keep, size=int(n_max), replace=False)
    return np.sort(keep)


def _load_one_sample(prefix: str, extracted_dir: Path, meta_idx: pd.DataFrame) -> ad.AnnData:
    mtx_path = extracted_dir / f"{prefix}.matrix.mtx.gz"
    barcodes_path = extracted_dir / f"{prefix}.barcodes.tsv.gz"
    features_path = extracted_dir / f"{prefix}.features.tsv.gz"

    with gzip.open(mtx_path, "rb") as f:
        mat = mmread(f).tocsr()
    with gzip.open(barcodes_path, "rt", encoding="utf-8", errors="ignore") as f:
        barcodes = [ln.strip() for ln in f if ln.strip()]
    with gzip.open(features_path, "rt", encoding="utf-8", errors="ignore") as f:
        feats = [ln.rstrip("\n").split("\t") for ln in f]

    if mat.shape[1] != len(barcodes):
        raise ValueError(f"Barcode count mismatch for {prefix}: matrix has {mat.shape[1]}, barcodes {len(barcodes)}")
    if mat.shape[0] != len(feats):
        raise ValueError(f"Feature count mismatch for {prefix}: matrix has {mat.shape[0]}, features {len(feats)}")

    gene_id = [x[0] if len(x) > 0 else "" for x in feats]
    gene_name = [x[1] if len(x) > 1 else (x[0] if len(x) > 0 else "") for x in feats]
    feature_type = [x[2] if len(x) > 2 else "" for x in feats]

    var_names = _dedupe(gene_name)
    sample_title = prefix.split("_", 1)[1] if "_" in prefix else prefix
    gsm = prefix.split("_", 1)[0]
    md = meta_idx.loc[sample_title] if sample_title in meta_idx.index else pd.Series(dtype=object)
    source_name = str(md.get("source_name_ch1", ""))
    time_value = str(md.get("time", ""))
    disease_state = str(md.get("disease_state", ""))

    adata = ad.AnnData(X=csr_matrix(mat.T))
    adata.var_names = var_names
    adata.var["gene_id"] = gene_id
    adata.var["gene_name"] = gene_name
    adata.var["feature_type"] = feature_type

    obs_names = [f"{prefix}:{bc}" for bc in barcodes]
    adata.obs_names = obs_names
    adata.obs["cell_barcode"] = barcodes
    adata.obs["sample_id"] = prefix
    adata.obs["geo_accession"] = gsm
    adata.obs["sample_title"] = sample_title
    adata.obs["source_name"] = source_name
    adata.obs["time"] = time_value
    adata.obs["time_group"] = _sample_time_group(time_value)
    adata.obs["disease_state"] = disease_state
    adata.obs["label"] = adata.obs["time_group"].astype(str)
    return adata


def _preprocess_with_scanpy(
    adata: ad.AnnData,
    min_genes: int,
    min_cells: int,
    n_hvg: int,
    n_pcs: int,
) -> ad.AnnData:
    os.environ.setdefault("NUMBA_DISABLE_COVERAGE", "1")
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
    import scanpy as sc

    a = adata.copy()
    sc.pp.filter_cells(a, min_genes=int(min_genes))
    sc.pp.filter_genes(a, min_cells=int(min_cells))
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    n_hvg_eff = min(int(n_hvg), max(200, a.n_vars))
    sc.pp.highly_variable_genes(a, n_top_genes=n_hvg_eff, flavor="seurat")
    if "highly_variable" in a.var.columns and a.var["highly_variable"].any():
        a = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(a, max_value=10)
    n_comp = min(int(n_pcs), max(2, a.n_obs - 1), max(2, a.n_vars - 1))
    sc.tl.pca(a, n_comps=n_comp, svd_solver="arpack")
    return a


def build_h5ad_from_gse165087(
    root: Path,
    force_download: bool = False,
    force_convert: bool = False,
    max_cells_raw: int | None = 40000,
    max_cells_processed: int | None = 12000,
    min_genes: int = 200,
    min_cells: int = 3,
    n_hvg: int = 2000,
    n_pcs: int = 50,
    seed: int = 42,
) -> dict[str, str]:
    raw_dir = ensure_dir(root / "raw")
    extracted_dir = ensure_dir(raw_dir / "extracted")
    out_raw_h5ad = raw_dir / "raw_combined.h5ad"
    out_processed_h5ad = root / "data.h5ad"
    out_sample_meta = root / "sample_metadata.tsv"
    out_readme = root / "README.txt"
    readme_text = (
        "CLL/RS scRNA public substitute dataset used by this repo:\n"
        "  GEO accession: GSE165087\n"
        "  Files downloaded:\n"
        "    raw/GSE165087_RAW.tar\n"
        "    raw/GSE165087_series_matrix.txt.gz\n"
        "    raw/filelist.txt\n"
        "  Generated files:\n"
        "    raw/raw_combined.h5ad\n"
        "    data.h5ad\n"
        "    sample_metadata.tsv\n"
        "\n"
        "data.h5ad includes metadata-derived labels:\n"
        "  label=time_group (pretreatment, relapse, on_ibrutinib, watchful_waiting, other)\n"
    )

    tar_path = raw_dir / f"{GSE}_RAW.tar"
    series_path = raw_dir / f"{GSE}_series_matrix.txt.gz"
    filelist_path = raw_dir / "filelist.txt"

    if force_download:
        for p in [tar_path, series_path, filelist_path]:
            if p.exists():
                p.unlink()

    download_url(RAW_TAR_URL, tar_path, timeout=900)
    download_url(SERIES_URL, series_path, timeout=300)
    download_url(FILELIST_URL, filelist_path, timeout=300)

    if force_convert:
        for p in [out_raw_h5ad, out_processed_h5ad, out_sample_meta]:
            if p.exists():
                p.unlink()

    has_extracted = any(extracted_dir.glob("*.matrix.mtx.gz"))
    if (not has_extracted) or force_convert:
        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(path=extracted_dir)

    if out_processed_h5ad.exists() and out_raw_h5ad.exists() and out_sample_meta.exists() and (not force_convert):
        out_readme.write_text(readme_text)
        return {
            "raw_tar": str(tar_path),
            "series_matrix": str(series_path),
            "raw_h5ad": str(out_raw_h5ad),
            "processed_h5ad": str(out_processed_h5ad),
            "sample_metadata": str(out_sample_meta),
            "readme_path": str(out_readme),
        }

    meta = parse_geo_series_matrix(series_path)
    meta = meta.drop_duplicates(subset=["sample_title"]).copy()
    meta_idx = meta.set_index("sample_title", drop=False)

    matrix_paths = sorted(extracted_dir.glob("*.matrix.mtx.gz"))
    if not matrix_paths:
        raise ValueError(f"No 10X matrix files found in {extracted_dir}")

    sample_prefixes = [p.name[: -len(".matrix.mtx.gz")] for p in matrix_paths]
    adatas = []
    sample_rows = []
    for prefix in sample_prefixes:
        adata = _load_one_sample(prefix, extracted_dir, meta_idx)
        adatas.append(adata)
        sample_rows.append(
            {
                "sample_id": prefix,
                "geo_accession": str(adata.obs["geo_accession"].iloc[0]),
                "sample_title": str(adata.obs["sample_title"].iloc[0]),
                "source_name": str(adata.obs["source_name"].iloc[0]),
                "time": str(adata.obs["time"].iloc[0]),
                "time_group": str(adata.obs["time_group"].iloc[0]),
                "disease_state": str(adata.obs["disease_state"].iloc[0]),
                "n_cells_raw": int(adata.n_obs),
            }
        )

    combined = ad.concat(adatas, join="inner", merge="same", label="batch", keys=sample_prefixes, index_unique=None)
    if not isinstance(combined.X, csr_matrix):
        combined.X = csr_matrix(combined.X)

    if max_cells_raw is not None and combined.n_obs > int(max_cells_raw):
        keep = _stratified_sample_obs(combined.obs, "time_group", int(max_cells_raw), int(seed))
        combined = combined[keep].copy()

    combined.write_h5ad(out_raw_h5ad)

    proc_in = combined
    if max_cells_processed is not None and proc_in.n_obs > int(max_cells_processed):
        keep = _stratified_sample_obs(proc_in.obs, "time_group", int(max_cells_processed), int(seed) + 13)
        proc_in = proc_in[keep].copy()

    processed = _preprocess_with_scanpy(
        proc_in,
        min_genes=int(min_genes),
        min_cells=int(min_cells),
        n_hvg=int(n_hvg),
        n_pcs=int(n_pcs),
    )
    processed.uns["source_dataset"] = GSE
    processed.uns["preprocessing"] = {
        "min_genes": int(min_genes),
        "min_cells": int(min_cells),
        "n_hvg": int(n_hvg),
        "n_pcs": int(min(int(n_pcs), processed.obsm["X_pca"].shape[1])),
        "max_cells_raw": None if max_cells_raw is None else int(max_cells_raw),
        "max_cells_processed": None if max_cells_processed is None else int(max_cells_processed),
        "seed": int(seed),
    }
    processed.write_h5ad(out_processed_h5ad)

    sample_meta = pd.DataFrame(sample_rows)
    sample_meta.to_csv(out_sample_meta, sep="\t", index=False)

    out_readme.write_text(readme_text)

    return {
        "raw_tar": str(tar_path),
        "series_matrix": str(series_path),
        "raw_h5ad": str(out_raw_h5ad),
        "processed_h5ad": str(out_processed_h5ad),
        "sample_metadata": str(out_sample_meta),
        "readme_path": str(out_readme),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Download and convert GSE165087 into CLL/RS scRNA h5ad substitute.")
    ap.add_argument("--root", default="data/cll_rs_scrna")
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--force-convert", action="store_true")
    ap.add_argument("--max-cells-raw", type=int, default=40000)
    ap.add_argument("--max-cells-processed", type=int, default=12000)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--min-cells", type=int, default=3)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    max_cells_raw = None if int(args.max_cells_raw) <= 0 else int(args.max_cells_raw)
    max_cells_processed = None if int(args.max_cells_processed) <= 0 else int(args.max_cells_processed)

    out = build_h5ad_from_gse165087(
        root=ensure_dir(Path(args.root)),
        force_download=bool(args.force_download),
        force_convert=bool(args.force_convert),
        max_cells_raw=max_cells_raw,
        max_cells_processed=max_cells_processed,
        min_genes=int(args.min_genes),
        min_cells=int(args.min_cells),
        n_hvg=int(args.n_hvg),
        n_pcs=int(args.n_pcs),
        seed=int(args.seed),
    )
    print(out)


if __name__ == "__main__":
    main()
