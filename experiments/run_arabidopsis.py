#!/usr/bin/env python3
"""CC-SPR evaluation on Arabidopsis root atlas (GSE123013).

Cross-domain generalization test: does geometry-sensitive topological
feature selection identify biologically coherent gene modules in plant
developmental scRNA-seq data?

Runs 3 backends (euclidean, ricci, dtm) + standard + harmonic baselines
with the same evaluation protocol as the cancer datasets.
"""
from __future__ import annotations

import json
import logging
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

_src = str(Path(__file__).resolve().parent.parent / "src")
sys.path.insert(0, _src)

from ccspr.datasets.arabidopsis_root import load_arabidopsis_root
from ccspr.eval.classification import FeatureParams, evaluate_ccspr_only
from ccspr.preprocess.basic import select_top_variable_genes
from ccspr.stability.tip import tip_bootstrap_topk

logging.basicConfig(level=logging.WARNING)

BACKENDS = ["euclidean", "ricci", "dtm"]
SEED = 42
CACHE_DIR = "data/cache"


def log(msg: str):
    print(f"[{datetime.utcnow().isoformat()}Z] {msg}", flush=True)


def ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)


def evaluate_backend(
    x, y, backend, dataset_name,
    var_genes=2000, n_splits=5, test_size=0.25, n_boot=10,
    k=10, iters=5, lambda_=1e-6, top_k=20, seed=SEED,
):
    """Evaluate a single backend using the standard protocol."""
    t0 = time.time()
    iters_use = iters if backend in ("ricci", "diffusion", "phate_like") else 0
    params = FeatureParams(
        mode=backend, solver="elastic", k=k, iters=iters_use,
        alpha=0.5, lambda_=lambda_, top_k=top_k, normalize="std",
        n_boot=n_boot, subsample_frac=0.8, rescale_median=1.0,
    )
    df = evaluate_ccspr_only(
        x, y, params=params, var_genes=var_genes,
        n_splits=n_splits, test_size=test_size,
        seed=seed, cache_dir=CACHE_DIR,
    )
    f1_values = df["f1_weighted"].values
    elapsed = time.time() - t0
    return {
        "dataset": dataset_name, "backend": backend,
        "mean_f1": float(np.mean(f1_values)),
        "std_f1": float(np.std(f1_values)) if len(f1_values) > 1 else 0.0,
        "n_splits": len(f1_values),
        "f1_values": f1_values.tolist(),
        "elapsed_sec": round(elapsed, 1),
    }


def collect_diagnostics(
    x_var, backend, dataset_name,
    k=10, iters=5, n_boot=10, lambda_=1e-6, top_k=20, seed=SEED,
):
    """Collect TIP diagnostics for a backend."""
    iters_use = iters if backend in ("ricci", "diffusion", "phate_like") else 0
    tip_out = tip_bootstrap_topk(
        x_var, n_boot=n_boot, subsample_frac=0.8, top_k=top_k,
        mode=backend, k=k, iters=iters_use, alpha=0.5,
        solver="elastic", lambda_=lambda_, normalize="std",
        rescale_median=1.0, seed=seed + 5000, cache_dir=CACHE_DIR,
    )
    tip_scores = tip_out["tip"]
    lifetimes = tip_out["lifetime"]
    prominences = tip_out["prominence"]

    support_size = int(np.count_nonzero(tip_scores > 0))
    tip_pos = tip_scores[tip_scores > 0]
    if len(tip_pos) > 0:
        p = tip_pos / tip_pos.sum()
        entropy = float(-np.sum(p * np.log(p + 1e-12)))
    else:
        entropy = 0.0
    if len(tip_pos) > 1:
        sorted_t = np.sort(tip_pos)
        n = len(sorted_t)
        index = np.arange(1, n + 1)
        gini = float(
            (2 * np.sum(index * sorted_t) - (n + 1) * np.sum(sorted_t))
            / (n * np.sum(sorted_t) + 1e-12)
        )
    else:
        gini = 0.0
    return {
        "dataset": dataset_name, "backend": backend,
        "n_ok": int(tip_out["n_ok"]),
        "mean_lifetime": float(np.mean(lifetimes)) if len(lifetimes) > 0 else 0.0,
        "std_lifetime": float(np.std(lifetimes)) if len(lifetimes) > 1 else 0.0,
        "mean_prominence": float(np.mean(prominences)) if len(prominences) > 0 else 0.0,
        "std_prominence": float(np.std(prominences)) if len(prominences) > 1 else 0.0,
        "support_size": support_size,
        "tip_entropy": entropy,
        "tip_gini": gini,
        "tip_scores": tip_scores.tolist(),
    }


def extract_tip_genes(
    x_var, feature_names, backend, dataset_name,
    k=10, iters=5, n_boot=10, lambda_=1e-6, top_k=20, seed=SEED,
):
    """Extract top TIP gene names for a backend."""
    iters_use = iters if backend in ("ricci", "diffusion", "phate_like") else 0
    tip_out = tip_bootstrap_topk(
        x_var, n_boot=n_boot, subsample_frac=0.8, top_k=top_k,
        mode=backend, k=k, iters=iters_use, alpha=0.5,
        solver="elastic", lambda_=lambda_, normalize="std",
        rescale_median=1.0, seed=seed + 7000, cache_dir=CACHE_DIR,
    )
    tip = tip_out["tip"]
    idx = np.argsort(-tip)[:top_k]
    genes = feature_names[idx].astype(str).tolist()
    scores = tip[idx].tolist()
    return genes, scores


def main():
    import argparse
    ap = argparse.ArgumentParser(description="CC-SPR Arabidopsis Root Evaluation")
    ap.add_argument("--data-path", default="data/arabidopsis_root/data.h5ad")
    ap.add_argument("--backends", nargs="+", default=BACKENDS)
    ap.add_argument("--max-cells", type=int, default=3000)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=10)
    ap.add_argument("--top-k", type=int, default=20)
    ap.add_argument("--var-genes", type=int, default=2000)
    ap.add_argument("--out-dir", default="results/arabidopsis")
    ap.add_argument("--use-hvg", action="store_true",
                     help="Use HVG expression instead of PCA for TIP gene extraction")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    ensure_dir(tables_dir)
    ensure_dir(figures_dir)

    log("=== Loading Arabidopsis Root Atlas ===")
    ds = load_arabidopsis_root(
        data_path=args.data_path,
        max_cells=args.max_cells,
        seed=SEED,
        use_pca=True,
    )
    x, y = ds["X"], ds["y"]
    log(f"  Shape: {x.shape}, classes: {len(np.unique(y))}")

    # Also load HVG version for gene name extraction
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

    # ================================================================
    # 1. Backend evaluation (classification)
    # ================================================================
    log("=== Backend Evaluation ===")
    all_benchmark = []
    for backend in args.backends:
        log(f"  {backend}: evaluating...")
        try:
            result = evaluate_backend(
                x, y, backend, "Arabidopsis_root",
                var_genes=min(args.var_genes, x.shape[1]),
                n_splits=args.n_splits, test_size=0.25,
                n_boot=args.n_boot, k=10, iters=2,
                lambda_=1e-6, top_k=args.top_k,
            )
            all_benchmark.append(result)
            log(f"  {backend}: F1={result['mean_f1']:.4f} ± {result['std_f1']:.4f} "
                f"({result['elapsed_sec']:.1f}s)")
        except Exception as e:
            log(f"  {backend}: FAILED — {e}")
            traceback.print_exc()

    # Standard baseline
    log("  standard baseline: evaluating...")
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import f1_score
    from ccspr.preprocess.basic import standardize

    sss = StratifiedShuffleSplit(n_splits=args.n_splits, test_size=0.25, random_state=SEED)
    std_f1s = []
    for itr, ite in sss.split(x, y):
        x_tr, x_te = standardize(x[itr], x[ite])
        clf = LogisticRegression(max_iter=5000, random_state=SEED)
        clf.fit(x_tr, y[itr])
        pred = clf.predict(x_te)
        std_f1s.append(f1_score(y[ite], pred, average="weighted"))
    all_benchmark.append({
        "dataset": "Arabidopsis_root", "backend": "standard",
        "mean_f1": float(np.mean(std_f1s)),
        "std_f1": float(np.std(std_f1s)),
        "n_splits": len(std_f1s),
        "f1_values": [float(f) for f in std_f1s],
        "elapsed_sec": 0.0,
    })
    log(f"  standard: F1={np.mean(std_f1s):.4f} ± {np.std(std_f1s):.4f}")

    # Save benchmark
    bench_df = pd.DataFrame([{
        "dataset": r["dataset"], "backend": r["backend"],
        "mean_f1": round(r["mean_f1"], 4), "std_f1": round(r["std_f1"], 4),
        "n_splits": r["n_splits"], "elapsed_sec": r.get("elapsed_sec", 0),
    } for r in all_benchmark])
    bench_df.to_csv(tables_dir / "arabidopsis_benchmark.csv", index=False)
    log(f"Saved: {tables_dir / 'arabidopsis_benchmark.csv'}")

    # ================================================================
    # 2. Diagnostics
    # ================================================================
    log("=== Diagnostics ===")
    all_diag = []
    for backend in args.backends:
        log(f"  {backend}: diagnostics...")
        try:
            diag = collect_diagnostics(
                x, backend, "Arabidopsis_root",
                k=10, iters=2, n_boot=args.n_boot,
                lambda_=1e-6, top_k=args.top_k,
            )
            all_diag.append(diag)
            log(f"  {backend}: lifetime={diag['mean_lifetime']:.4f}, "
                f"support={diag['support_size']}, gini={diag['tip_gini']:.3f}")
        except Exception as e:
            log(f"  {backend}: diagnostics FAILED — {e}")

    diag_df = pd.DataFrame([{
        k: v for k, v in d.items() if k != "tip_scores"
    } for d in all_diag])
    diag_df.to_csv(tables_dir / "arabidopsis_diagnostics.csv", index=False)

    # ================================================================
    # 3. TIP gene extraction (using HVG matrix)
    # ================================================================
    log("=== TIP Gene Extraction ===")
    all_gene_rows = []
    for backend in args.backends:
        log(f"  {backend}: extracting top genes...")
        try:
            genes, scores = extract_tip_genes(
                x_hvg, gene_names, backend, "Arabidopsis_root",
                k=10, iters=2, n_boot=args.n_boot,
                lambda_=1e-6, top_k=args.top_k,
            )
            for rank, (g, s) in enumerate(zip(genes, scores), 1):
                all_gene_rows.append({
                    "backend": backend, "rank": rank,
                    "gene": g, "tip_score": round(s, 4),
                })
            log(f"  {backend}: top genes = {genes[:5]}")
        except Exception as e:
            log(f"  {backend}: gene extraction FAILED — {e}")
            traceback.print_exc()

    gene_df = pd.DataFrame(all_gene_rows)
    gene_df.to_csv(tables_dir / "arabidopsis_top_features_by_backend.csv", index=False)

    # ================================================================
    # 4. Figures
    # ================================================================
    log("=== Generating Figures ===")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # F1 bar chart
    if all_benchmark:
        fig, ax = plt.subplots(figsize=(7, 4))
        methods = [r["backend"] for r in all_benchmark]
        f1s = [r["mean_f1"] for r in all_benchmark]
        stds = [r["std_f1"] for r in all_benchmark]
        colors = ["#4C72B0" if m != "standard" else "#DD8452" for m in methods]
        bars = ax.bar(range(len(methods)), f1s, yerr=stds, capsize=4, color=colors)
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(methods, rotation=30, ha="right")
        ax.set_ylabel("Weighted F1")
        ax.set_title("Arabidopsis Root: Backend Comparison")
        ax.set_ylim(0, 1.05)
        fig.tight_layout()
        fig.savefig(figures_dir / "arabidopsis_f1_bars.png", dpi=220)
        plt.close(fig)

    # TIP comparison across backends
    if len(all_diag) >= 2:
        fig, axes = plt.subplots(1, len(all_diag), figsize=(4 * len(all_diag), 4), sharey=True)
        if len(all_diag) == 1:
            axes = [axes]
        for ax, d in zip(axes, all_diag):
            tip = np.array(d["tip_scores"])
            top_idx = np.argsort(-tip)[:20]
            ax.bar(range(20), tip[top_idx])
            ax.set_title(f"{d['backend']}\nsupport={d['support_size']}, gini={d['tip_gini']:.3f}")
            ax.set_xlabel("Feature rank")
            if ax == axes[0]:
                ax.set_ylabel("TIP score")
        fig.suptitle("Arabidopsis Root: TIP Profiles", fontsize=12)
        fig.tight_layout()
        fig.savefig(figures_dir / "arabidopsis_tip_profiles.png", dpi=220)
        plt.close(fig)

    # ================================================================
    # 5. Save full results JSON
    # ================================================================
    full_results = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "dataset": "Arabidopsis_root_GSE123013",
        "paper": "Denyer et al. 2019",
        "benchmark": all_benchmark,
        "diagnostics": [{k: v for k, v in d.items()} for d in all_diag],
    }
    with open(out_dir / "arabidopsis_results.json", "w") as f:
        json.dump(full_results, f, indent=2, default=str)

    log("=== SUMMARY ===")
    print(bench_df.to_string(index=False))
    log("=== DONE ===")


if __name__ == "__main__":
    main()
