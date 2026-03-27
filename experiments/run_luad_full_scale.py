#!/usr/bin/env python3
"""Full-scale LUAD evaluation with all available samples.

Removes the max_samples=60 cap from the geometry_family run to test
whether backend rankings hold at full TCGA-LUAD scale (~500 samples).

Runs 3 key backends: euclidean, ricci, dtm (the ones that showed
differentiation at n=60).
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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ccspr.datasets.tcga_luad import download_tcga_luad, load_tcga_luad
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
    x, y, backend, var_genes=5000, n_splits=5, test_size=0.25,
    n_boot=10, k=10, iters=5, lambda_=1e-6, top_k=20, seed=SEED,
):
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
        "backend": backend,
        "mean_f1": float(np.mean(f1_values)),
        "std_f1": float(np.std(f1_values)) if len(f1_values) > 1 else 0.0,
        "n_splits": len(f1_values),
        "f1_values": f1_values.tolist(),
        "elapsed_sec": round(elapsed, 1),
    }


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Full-scale LUAD evaluation")
    ap.add_argument("--backends", nargs="+", default=BACKENDS)
    ap.add_argument("--max-samples", type=int, default=0,
                     help="0 = no cap (use all available)")
    ap.add_argument("--var-genes", type=int, default=5000)
    ap.add_argument("--n-splits", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=10)
    ap.add_argument("--top-k", type=int, default=20)
    ap.add_argument("--out-dir", default="results/luad_full")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    tables_dir = out_dir / "tables"
    ensure_dir(tables_dir)

    log("=== Loading Full TCGA-LUAD ===")
    download_tcga_luad("data/tcga_luad")
    ds = load_tcga_luad(root="data/tcga_luad", min_class_size=10)
    x, y = ds["X"], ds["y"]
    feat = ds["feature_names"]
    log(f"  Full LUAD: {x.shape[0]} samples, {x.shape[1]} features, "
        f"{len(np.unique(y))} classes")

    if args.max_samples > 0 and x.shape[0] > args.max_samples:
        rng = np.random.default_rng(SEED)
        keep = []
        for c in np.unique(y):
            idx = np.where(y == c)[0]
            take = max(2, int(round(len(idx) * (args.max_samples / len(y)))))
            take = min(take, len(idx))
            keep.extend(rng.choice(idx, take, replace=False).tolist())
        keep = np.array(sorted(set(keep)))[:args.max_samples]
        x, y = x[keep], y[keep]
        log(f"  Subsampled to {x.shape[0]} samples")

    # Standard baseline
    log("=== Standard Baseline ===")
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import f1_score
    from ccspr.preprocess.basic import standardize

    x_var, idx_var = select_top_variable_genes(x, top_n=args.var_genes)
    sss = StratifiedShuffleSplit(n_splits=args.n_splits, test_size=0.25, random_state=SEED)
    std_f1s = []
    for itr, ite in sss.split(x_var, y):
        xtr, xte = standardize(x_var[itr], x_var[ite])
        clf = LogisticRegression(max_iter=5000, random_state=SEED)
        clf.fit(xtr, y[itr])
        pred = clf.predict(xte)
        std_f1s.append(f1_score(y[ite], pred, average="weighted"))
    log(f"  Standard: F1={np.mean(std_f1s):.4f} ± {np.std(std_f1s):.4f}")

    all_results = [{
        "backend": "standard", "mean_f1": float(np.mean(std_f1s)),
        "std_f1": float(np.std(std_f1s)), "n_splits": len(std_f1s),
        "f1_values": [float(f) for f in std_f1s], "elapsed_sec": 0.0,
    }]

    # Backend evaluations
    for backend in args.backends:
        log(f"=== {backend} ===")
        try:
            result = evaluate_backend(
                x, y, backend, var_genes=args.var_genes,
                n_splits=args.n_splits, n_boot=args.n_boot,
                top_k=args.top_k,
            )
            all_results.append(result)
            log(f"  {backend}: F1={result['mean_f1']:.4f} ± {result['std_f1']:.4f} "
                f"({result['elapsed_sec']:.1f}s)")
        except Exception as e:
            log(f"  {backend}: FAILED — {e}")
            traceback.print_exc()

    # TIP gene extraction (for enrichment)
    log("=== TIP Gene Extraction ===")
    gene_rows = []
    for backend in args.backends:
        log(f"  {backend}: extracting genes...")
        try:
            iters_use = 5 if backend in ("ricci", "diffusion", "phate_like") else 0
            tip_out = tip_bootstrap_topk(
                x_var, n_boot=args.n_boot, subsample_frac=0.8,
                top_k=args.top_k, mode=backend, k=10, iters=iters_use,
                alpha=0.5, solver="elastic", lambda_=1e-6, normalize="std",
                rescale_median=1.0, seed=SEED + 7000, cache_dir=CACHE_DIR,
            )
            tip = tip_out["tip"]
            idx = np.argsort(-tip)[:args.top_k]
            genes_var = feat[idx_var]
            for rank, i in enumerate(idx, 1):
                gene_rows.append({
                    "method": f"{backend}_elastic", "rank": rank,
                    "gene": str(genes_var[i]), "tip": float(tip[i]),
                })
            log(f"  {backend}: top genes = {[str(genes_var[i]) for i in idx[:5]]}")
        except Exception as e:
            log(f"  {backend}: gene extraction FAILED — {e}")

    # Save results
    bench_df = pd.DataFrame(all_results)
    bench_df.to_csv(tables_dir / "luad_full_benchmark.csv", index=False)
    log(f"Saved: {tables_dir / 'luad_full_benchmark.csv'}")

    if gene_rows:
        gene_df = pd.DataFrame(gene_rows)
        gene_df.to_csv(tables_dir / "luad_full_top_features_by_method.csv", index=False)

    with open(out_dir / "luad_full_results.json", "w") as f:
        json.dump({
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "n_samples": int(x.shape[0]),
            "n_features": int(x.shape[1]),
            "n_classes": int(len(np.unique(y))),
            "benchmark": all_results,
        }, f, indent=2, default=str)

    log("=== SUMMARY ===")
    print(bench_df[["backend", "mean_f1", "std_f1", "n_splits"]].to_string(index=False))
    log("=== DONE ===")


if __name__ == "__main__":
    main()
