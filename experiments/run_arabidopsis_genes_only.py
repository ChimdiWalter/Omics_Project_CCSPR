#!/usr/bin/env python3
"""Extract TIP genes + generate figures for Arabidopsis.

Lightweight companion to run_arabidopsis.py — skips benchmark and
diagnostics (already saved) and only does the missing gene extraction
and figure generation.
"""
from __future__ import annotations

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

_src = str(Path(__file__).resolve().parent.parent / "src")
sys.path.insert(0, _src)

from ccspr.datasets.arabidopsis_root import load_arabidopsis_root
from ccspr.stability.tip import tip_bootstrap_topk

BACKENDS = ["euclidean", "dtm"]
SEED = 42
CACHE_DIR = "data/cache"


def log(msg: str):
    print(f"[{datetime.utcnow().isoformat()}Z] {msg}", flush=True)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-path", default="data/arabidopsis_root/data.h5ad")
    ap.add_argument("--backends", nargs="+", default=BACKENDS)
    ap.add_argument("--max-cells", type=int, default=500)
    ap.add_argument("--n-boot", type=int, default=3)
    ap.add_argument("--top-k", type=int, default=20)
    ap.add_argument("--var-genes", type=int, default=2000)
    ap.add_argument("--out-dir", default="results/arabidopsis")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Load HVG matrix for gene name extraction
    log("Loading Arabidopsis HVG matrix...")
    ds_hvg = load_arabidopsis_root(
        data_path=args.data_path,
        max_cells=args.max_cells,
        seed=SEED,
        use_pca=False,
        n_hvg=args.var_genes,
    )
    x_hvg = ds_hvg["X"]
    gene_names = ds_hvg["feature_names"]
    log(f"  HVG matrix: {x_hvg.shape}, gene names: {len(gene_names)}")

    # TIP gene extraction
    log("=== TIP Gene Extraction ===")
    all_gene_rows = []
    for backend in args.backends:
        log(f"  {backend}: extracting top genes...")
        t0 = time.time()
        iters_use = 2 if backend in ("ricci", "diffusion", "phate_like") else 0
        tip_out = tip_bootstrap_topk(
            x_hvg, n_boot=args.n_boot, subsample_frac=0.8,
            top_k=args.top_k, mode=backend, k=10, iters=iters_use,
            alpha=0.5, solver="elastic", lambda_=1e-6, normalize="std",
            rescale_median=1.0, seed=SEED + 7000, cache_dir=CACHE_DIR,
        )
        tip = tip_out["tip"]
        idx = np.argsort(-tip)[:args.top_k]
        genes = gene_names[idx].astype(str).tolist()
        scores = tip[idx].tolist()
        for rank, (g, s) in enumerate(zip(genes, scores), 1):
            all_gene_rows.append({
                "backend": backend, "rank": rank,
                "gene": g, "tip_score": round(s, 4),
            })
        elapsed = time.time() - t0
        log(f"  {backend}: top genes = {genes[:5]} ({elapsed:.0f}s)")

    gene_df = pd.DataFrame(all_gene_rows)
    gene_df.to_csv(tables_dir / "arabidopsis_top_features_by_backend.csv", index=False)
    log(f"Saved: {tables_dir / 'arabidopsis_top_features_by_backend.csv'}")

    # Figures
    log("=== Generating Figures ===")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # F1 bar chart from existing benchmark
    bench_path = tables_dir / "arabidopsis_benchmark.csv"
    if bench_path.exists():
        bench_df = pd.read_csv(bench_path)
        fig, ax = plt.subplots(figsize=(7, 4))
        methods = bench_df["backend"].tolist()
        f1s = bench_df["mean_f1"].tolist()
        stds = bench_df["std_f1"].tolist()
        colors = ["#4C72B0" if m != "standard" else "#DD8452" for m in methods]
        ax.bar(range(len(methods)), f1s, yerr=stds, capsize=4, color=colors)
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(methods, rotation=30, ha="right")
        ax.set_ylabel("Weighted F1")
        ax.set_title("Arabidopsis Root: Backend Comparison")
        ax.set_ylim(0, 1.05)
        fig.tight_layout()
        fig.savefig(figures_dir / "arabidopsis_f1_bars.png", dpi=220)
        plt.close(fig)
        log(f"Saved: {figures_dir / 'arabidopsis_f1_bars.png'}")

    # TIP profiles from existing diagnostics
    diag_path = tables_dir / "arabidopsis_diagnostics.csv"
    if diag_path.exists():
        diag_df = pd.read_csv(diag_path)
        fig, ax = plt.subplots(figsize=(6, 4))
        x_pos = np.arange(len(diag_df))
        ax.bar(x_pos - 0.15, diag_df["mean_lifetime"], 0.3, label="Lifetime", color="#4C72B0")
        ax.bar(x_pos + 0.15, diag_df["mean_prominence"], 0.3, label="Prominence", color="#DD8452")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(diag_df["backend"], rotation=30, ha="right")
        ax.set_ylabel("Mean value")
        ax.set_title("Arabidopsis Root: TIP Diagnostics")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures_dir / "arabidopsis_diagnostics.png", dpi=220)
        plt.close(fig)
        log(f"Saved: {figures_dir / 'arabidopsis_diagnostics.png'}")

    # Gene ranking comparison
    if all_gene_rows:
        backends_found = gene_df["backend"].unique()
        if len(backends_found) >= 2:
            sets = {b: set(gene_df[gene_df["backend"] == b]["gene"]) for b in backends_found}
            shared = set.intersection(*sets.values())
            fig, ax = plt.subplots(figsize=(6, 4))
            for i, b in enumerate(backends_found):
                unique = sets[b] - set.union(*[s for k, s in sets.items() if k != b])
                ax.bar(i, len(unique), color="#4C72B0", label="Unique" if i == 0 else "")
                ax.bar(i, len(shared), bottom=len(unique), color="#DD8452",
                       label="Shared" if i == 0 else "")
            ax.set_xticks(range(len(backends_found)))
            ax.set_xticklabels(backends_found, rotation=30, ha="right")
            ax.set_ylabel("Number of genes")
            ax.set_title(f"Arabidopsis Root: Gene Overlap (top {args.top_k})")
            ax.legend()
            fig.tight_layout()
            fig.savefig(figures_dir / "arabidopsis_gene_overlap.png", dpi=220)
            plt.close(fig)
            log(f"Saved: {figures_dir / 'arabidopsis_gene_overlap.png'}")

    log("=== DONE ===")


if __name__ == "__main__":
    main()
