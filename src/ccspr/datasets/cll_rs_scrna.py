from __future__ import annotations

from pathlib import Path

import numpy as np


def load_cll_rs_scrna(
    data_path: str | Path,
    label_key: str = "label",
    min_genes: int = 200,
    min_cells: int = 3,
    n_hvg: int = 2000,
    n_pcs: int = 50,
    max_cells: int | None = 5000,
    seed: int = 0,
    use_preprocessed_if_available: bool = True,
) -> dict:
    try:
        import scanpy as sc
    except Exception as e:
        raise ImportError("scanpy is required for CLL/RS scRNA loading. Install scanpy and anndata.") from e

    data_path = Path(data_path)
    if not data_path.exists():
        raise FileNotFoundError(
            f"scRNA path not found: {data_path}. "
            "For the RS/CLL study data, request controlled access first "
            "(dbGaP phs002458.v2.p1; EGA EGAS00001005495 / EGAD00001007922)."
        )

    if data_path.suffix == ".h5ad":
        adata = sc.read_h5ad(data_path)
    else:
        adata = sc.read_10x_mtx(data_path, var_names="gene_symbols", make_unique=True)

    if (
        bool(use_preprocessed_if_available)
        and data_path.suffix == ".h5ad"
        and "X_pca" in adata.obsm
        and label_key in adata.obs.columns
    ):
        if max_cells is not None and adata.n_obs > int(max_cells):
            rng = np.random.default_rng(int(seed))
            idx = rng.choice(adata.n_obs, int(max_cells), replace=False)
            adata = adata[idx].copy()

        x_pca = np.asarray(adata.obsm["X_pca"], dtype=float)
        n_use = min(int(n_pcs), x_pca.shape[1]) if int(n_pcs) > 0 else x_pca.shape[1]
        x = x_pca[:, :n_use]
        y = adata.obs[label_key].astype(str).to_numpy()
        return {
            "X": np.nan_to_num(x),
            "y": y,
            "sample_ids": adata.obs_names.astype(str).to_numpy(),
            "feature_names": np.array([f"PC{i+1}" for i in range(x.shape[1])], dtype=object),
            "meta": {
                "n_cells": int(adata.n_obs),
                "n_features": int(x.shape[1]),
                "label_key": label_key,
                "source": str(data_path),
                "preprocessed": True,
            },
        }

    if max_cells is not None and adata.n_obs > int(max_cells):
        rng = np.random.default_rng(int(seed))
        idx = rng.choice(adata.n_obs, int(max_cells), replace=False)
        adata = adata[idx].copy()

    sc.pp.filter_cells(adata, min_genes=int(min_genes))
    sc.pp.filter_genes(adata, min_cells=int(min_cells))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=int(n_hvg), flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(int(n_pcs), adata.n_vars - 1))

    if label_key not in adata.obs.columns:
        raise ValueError(
            f"label_key '{label_key}' not found in adata.obs. Available: {list(adata.obs.columns)}"
        )

    x = adata.obsm["X_pca"].astype(float)
    y = adata.obs[label_key].astype(str).to_numpy()

    return {
        "X": np.nan_to_num(x),
        "y": y,
        "sample_ids": adata.obs_names.astype(str).to_numpy(),
        "feature_names": np.array([f"PC{i+1}" for i in range(x.shape[1])], dtype=object),
        "meta": {
            "n_cells": int(adata.n_obs),
            "n_features": int(x.shape[1]),
            "label_key": label_key,
            "source": str(data_path),
        },
    }
