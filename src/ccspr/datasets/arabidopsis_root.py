"""Arabidopsis root atlas scRNA-seq dataset loader.

Source: Denyer et al. 2019, Developmental Cell 48(6):840-852.e5
GEO: GSE123013

Returns the same dict interface as other ccspr dataset loaders:
  {X, y, sample_ids, feature_names, meta}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def load_arabidopsis_root(
    data_path: str | Path = "data/arabidopsis_root/data.h5ad",
    label_key: str = "leiden",
    n_hvg: int = 2000,
    n_pcs: int = 50,
    max_cells: int = 3000,
    min_cluster_size: int = 20,
    seed: int = 42,
    use_pca: bool = True,
) -> dict:
    """Load preprocessed Arabidopsis root atlas data.

    Parameters
    ----------
    data_path : path to .h5ad file produced by download_arabidopsis_root.py
    label_key : obs column for labels (default: 'leiden' clusters)
    n_hvg : number of highly variable genes to select (if re-selecting)
    n_pcs : number of PCs for feature matrix
    max_cells : maximum cells to retain (stratified subsampling)
    min_cluster_size : drop clusters smaller than this
    seed : random seed for subsampling
    use_pca : if True, use PCA embeddings as X; if False, use HVG expression

    Returns
    -------
    dict with keys: X, y, sample_ids, feature_names, meta
    """
    import anndata as ad

    data_path = Path(data_path)
    if not data_path.exists():
        raise FileNotFoundError(
            f"Arabidopsis data not found at {data_path}. "
            "Run: python3 scripts/download_arabidopsis_root.py"
        )

    adata = ad.read_h5ad(data_path)

    # Use specified label column
    if label_key not in adata.obs.columns:
        raise KeyError(
            f"Label column '{label_key}' not in obs. "
            f"Available: {adata.obs.columns.tolist()}"
        )

    labels = adata.obs[label_key].astype(str).values

    # Filter small clusters
    from collections import Counter
    counts = Counter(labels)
    keep_labels = {l for l, c in counts.items() if c >= min_cluster_size}
    mask = np.array([l in keep_labels for l in labels])
    adata = adata[mask].copy()
    labels = labels[mask]
    print(f"  After filtering clusters < {min_cluster_size}: {adata.n_obs} cells, "
          f"{len(keep_labels)} clusters")

    # Stratified subsampling if needed
    if max_cells > 0 and adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        keep_idx = []
        unique_labels = np.unique(labels)
        for lab in unique_labels:
            lab_idx = np.where(labels == lab)[0]
            n_take = max(2, int(round(len(lab_idx) * (max_cells / adata.n_obs))))
            n_take = min(n_take, len(lab_idx))
            keep_idx.extend(rng.choice(lab_idx, n_take, replace=False).tolist())
        keep_idx = np.array(sorted(set(keep_idx)), dtype=int)
        if len(keep_idx) > max_cells:
            keep_idx = rng.choice(keep_idx, max_cells, replace=False)
            keep_idx = np.sort(keep_idx)
        adata = adata[keep_idx].copy()
        labels = labels[keep_idx]
        print(f"  After subsampling: {adata.n_obs} cells")

    # Build feature matrix
    if use_pca and "X_pca" in adata.obsm:
        X = adata.obsm["X_pca"][:, :n_pcs].astype(np.float64)
        feature_names = np.array([f"PC{i+1}" for i in range(X.shape[1])])
    else:
        # Use HVG expression
        if "highly_variable" in adata.var.columns:
            hvg_mask = adata.var["highly_variable"].values
            X = adata[:, hvg_mask].X
            feature_names = adata.var_names[hvg_mask].values.astype(str)
        else:
            X = adata.X
            feature_names = adata.var_names.values.astype(str)

        if hasattr(X, "toarray"):
            X = X.toarray()
        X = X.astype(np.float64)

        # Limit to n_hvg if needed
        if X.shape[1] > n_hvg:
            var = np.var(X, axis=0)
            top_idx = np.argsort(-var)[:n_hvg]
            X = X[:, top_idx]
            feature_names = feature_names[top_idx]

    # Encode labels as strings
    y = labels.astype(str)

    sample_ids = adata.obs_names.values.astype(str)

    meta = {
        "source": "GSE123013",
        "organism": "Arabidopsis thaliana",
        "tissue": "root",
        "paper": "Denyer et al. 2019",
        "n_cells": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "n_classes": int(len(np.unique(y))),
        "label_key": label_key,
        "use_pca": use_pca,
    }

    print(f"  Arabidopsis root: {X.shape[0]} cells, {X.shape[1]} features, "
          f"{meta['n_classes']} classes")

    return {
        "X": X,
        "y": y,
        "sample_ids": sample_ids,
        "feature_names": feature_names,
        "meta": meta,
    }
