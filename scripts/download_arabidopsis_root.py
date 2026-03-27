#!/usr/bin/env python3
"""Download Arabidopsis root atlas scRNA-seq data from GEO (GSE123013).

Denyer et al. 2019: Spatiotemporal Developmental Trajectories in the
Arabidopsis Root Revealed Using High-Throughput Single-Cell RNA Sequencing.
Developmental Cell 48(6):840-852.e5.

Downloads the 5-way merged raw count matrix and produces a processed
AnnData object with cell type annotations.
"""
from __future__ import annotations

import argparse
import gzip
import os
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


MERGE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123013/suppl/"
    "GSE123013%5F5way%5Fmerge%5Fraw%2Etsv%2Egz"
)

FALLBACK_URL = (
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123013/suppl/"
    "GSE123013%5F5way%5Fmerge%5Fraw%2Etsv%2Egz"
)


def download_file(url: str, out_path: Path, timeout: int = 300) -> None:
    """Download a file with progress reporting."""
    if out_path.exists():
        print(f"  Already exists: {out_path}")
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Downloading: {url}")
    print(f"  Destination: {out_path}")
    try:
        urllib.request.urlretrieve(url, str(out_path))
    except Exception as e:
        print(f"  Primary URL failed: {e}")
        print(f"  Trying fallback...")
        urllib.request.urlretrieve(FALLBACK_URL, str(out_path))
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  Downloaded: {size_mb:.1f} MB")


def build_anndata(raw_tsv_gz: Path, root: Path, wt_only: bool = True) -> Path:
    """Parse the 5-way merged count matrix into AnnData.

    The TSV has genes as rows and cells as columns. Cell barcodes encode
    the sample of origin.
    """
    import anndata as ad
    import scanpy as sc

    out_path = root / "data.h5ad"
    if out_path.exists():
        print(f"  AnnData already exists: {out_path}")
        return out_path

    print("  Reading count matrix (this may take a minute)...")
    df = pd.read_csv(raw_tsv_gz, sep="\t", index_col=0)
    print(f"  Raw matrix shape: {df.shape[0]} genes x {df.shape[1]} cells")

    # Clean quoted names
    df.index = df.index.astype(str).str.strip('"')
    df.columns = df.columns.astype(str).str.strip('"')

    # Transpose: cells as rows, genes as columns
    adata = ad.AnnData(X=df.T.values.astype(np.float32))
    adata.obs_names = pd.Index(df.columns.astype(str))
    adata.var_names = pd.Index(df.index.astype(str))

    # Infer sample of origin from barcode
    # Format: "SAMPLE_SAMPLE_BARCODE" e.g. "WT1_WT1_AAACCTGAGACAGACC"
    # or quoted: '"WT1_WT1_AAACCTGAGACAGACC"'
    barcodes = adata.obs_names.tolist()

    sample_ids = []
    for bc in barcodes:
        bc_clean = bc.strip('"')
        parts = bc_clean.split("_")
        if len(parts) >= 2:
            sample_ids.append(parts[0])  # e.g. "WT1", "WT2", "WT3", "rhd6", "gl2"
        else:
            sample_ids.append("unknown")

    adata.obs["sample_id"] = sample_ids
    unique_samples = sorted(adata.obs["sample_id"].unique())
    counts_per_sample = adata.obs["sample_id"].value_counts()
    print(f"  Unique sample IDs: {unique_samples}")
    for s in unique_samples:
        print(f"    {s}: {counts_per_sample[s]} cells")

    if wt_only:
        # Keep only WT replicates
        wt_ids = [s for s in unique_samples if s.upper().startswith("WT")]
        if wt_ids:
            print(f"  Filtering to WT samples: {wt_ids}")
            adata = adata[adata.obs["sample_id"].isin(wt_ids)].copy()
            print(f"  WT cells: {adata.n_obs}")
        else:
            print(f"  Warning: Could not identify WT samples from IDs {unique_samples}")
            print(f"  Keeping all cells")

    # Basic QC filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After QC filtering: {adata.n_obs} cells, {adata.n_vars} genes")

    # Normalize and log-transform
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars))
    print(f"  HVGs: {adata.var['highly_variable'].sum()}")

    # PCA
    sc.tl.pca(adata, n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1))

    # Neighbors + Leiden clustering for cell type assignment
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, adata.obsm["X_pca"].shape[1]))
    sc.tl.leiden(adata, resolution=0.8)
    print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

    # UMAP for visualization (optional but useful)
    sc.tl.umap(adata)

    # Save
    adata.write(out_path)
    print(f"  Saved: {out_path} ({out_path.stat().st_size / (1024*1024):.1f} MB)")

    # Save metadata
    meta = adata.obs[["sample_id", "leiden"]].copy()
    meta["barcode"] = adata.obs_names
    meta.to_csv(root / "cell_metadata.csv", index=False)
    print(f"  Saved: {root / 'cell_metadata.csv'}")

    return out_path


def main():
    ap = argparse.ArgumentParser(description="Download Arabidopsis root atlas (GSE123013)")
    ap.add_argument("--root", default="data/arabidopsis_root", help="Output directory")
    ap.add_argument("--all-samples", action="store_true", help="Keep all 5 samples (default: WT only)")
    args = ap.parse_args()

    root = Path(args.root)
    raw_dir = root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    print("=== Downloading Arabidopsis Root Atlas (GSE123013) ===")
    tsv_gz = raw_dir / "GSE123013_5way_merge_raw.tsv.gz"
    download_file(MERGE_URL, tsv_gz)

    print("\n=== Building AnnData ===")
    build_anndata(tsv_gz, root, wt_only=not args.all_samples)

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
