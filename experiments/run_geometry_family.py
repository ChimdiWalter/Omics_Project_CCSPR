#!/usr/bin/env python3
"""CC-SPR Geometry-Family Empirical Study.

Runs all 5 geometry backends (euclidean, ricci, diffusion, phate_like, dtm)
across LUAD, GSE161711, and bounded GSE165087. Collects F1, H1 diagnostics,
TIP scores, and representative support diagnostics.

Resource profile:
  LUAD:      max_samples=60, n_boot=5, 2-split (matches original runs)
  GSE161711: n=96, n_boot=5, 5-fold CV
  GSE165087: bounded 320 cells, n_boot=1, 1-split
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

# Ensure src/ is on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ccspr.datasets.tcga_luad import download_tcga_luad, load_tcga_luad
from ccspr.datasets.cll_venetoclax import load_cll_venetoclax
from ccspr.eval.classification import FeatureParams, evaluate_ccspr_only
from ccspr.geometry.distance import VALID_MODES
from ccspr.preprocess.basic import select_top_variable_genes
from ccspr.stability.tip import tip_bootstrap_topk
from ccspr.utils.io import ensure_dir

logging.basicConfig(level=logging.WARNING)

# ---------------------------------------------------------------------------
BACKENDS = ["euclidean", "ricci", "diffusion", "phate_like", "dtm"]
SEED = 42
CACHE_DIR = "data/cache"

RESULTS_DIR = Path("results/geometry_family")
TABLES_DIR = RESULTS_DIR / "tables"


def log(msg: str):
    print(f"[{datetime.utcnow().isoformat()}Z] {msg}", flush=True)


def _stratified_subsample(x, y, n_max, seed):
    if x.shape[0] <= n_max:
        return x, y
    rng = np.random.default_rng(seed)
    keep = []
    classes, counts = np.unique(y, return_counts=True)
    for c, cnt in zip(classes, counts):
        idx = np.where(y == c)[0]
        take = max(2, int(round(len(idx) * (n_max / len(y)))))
        take = min(take, len(idx))
        keep.extend(rng.choice(idx, take, replace=False).tolist())
    keep = np.array(sorted(set(keep)), dtype=int)
    if len(keep) > n_max:
        keep = rng.choice(keep, n_max, replace=False)
    return x[keep], y[keep]


def evaluate_backend(
    x, y, backend, dataset_name,
    var_genes=5000, n_splits=2, test_size=0.25, n_boot=5,
    k=10, iters=5, lambda_=1e-6, top_k=20, seed=SEED,
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
        "dataset": dataset_name, "backend": backend,
        "mean_f1": float(np.mean(f1_values)),
        "std_f1": float(np.std(f1_values)) if len(f1_values) > 1 else 0.0,
        "n_splits": len(f1_values),
        "f1_values": f1_values.tolist(),
        "elapsed_sec": round(elapsed, 1),
    }


def collect_diagnostics(
    x_var, backend, dataset_name,
    k=10, iters=5, n_boot=5, lambda_=1e-6, top_k=20, seed=SEED,
):
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
        gini = float((2 * np.sum(index * sorted_t) - (n + 1) * np.sum(sorted_t)) / (n * np.sum(sorted_t) + 1e-12))
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
        "lifetimes": lifetimes.tolist(),
        "prominences": prominences.tolist(),
    }


def main():
    ensure_dir(RESULTS_DIR)
    ensure_dir(TABLES_DIR)
    ensure_dir(Path(CACHE_DIR))

    all_benchmark = []
    all_diagnostics = []
    all_ablation = []
    timing = {}

    def save_checkpoint(phase):
        with open(RESULTS_DIR / "checkpoint.json", "w") as f:
            json.dump({
                "phase": phase,
                "timestamp": datetime.utcnow().isoformat() + "Z",
                "benchmark": all_benchmark,
                "diagnostics": all_diagnostics,
                "ablation": all_ablation,
                "timing": timing,
            }, f, indent=2, default=str)
        log(f"Checkpoint saved: phase={phase}")

    # ===================================================================
    # 1. LUAD — max_samples=60, n_boot=5, 2-split
    # ===================================================================
    log("=== PHASE 1: LUAD Geometry-Family Benchmark ===")
    download_tcga_luad("data/tcga_luad")
    ds = load_tcga_luad(root="data/tcga_luad", min_class_size=10)
    x_luad, y_luad = ds["X"], ds["y"]
    x_luad, y_luad = _stratified_subsample(x_luad, y_luad, 60, SEED)
    log(f"LUAD loaded+subsampled: {x_luad.shape[0]} samples, {x_luad.shape[1]} features")

    x_luad_var, _ = select_top_variable_genes(x_luad, top_n=5000)

    for backend in BACKENDS:
        log(f"  LUAD / {backend}: evaluating...")
        try:
            result = evaluate_backend(
                x_luad, y_luad, backend, "LUAD",
                var_genes=5000, n_splits=2, test_size=0.25, n_boot=5,
                k=10, iters=5, lambda_=1e-6, top_k=20,
            )
            all_benchmark.append(result)
            timing[f"LUAD_{backend}"] = result["elapsed_sec"]
            log(f"  LUAD / {backend}: F1={result['mean_f1']:.4f} ± {result['std_f1']:.4f} ({result['elapsed_sec']:.1f}s)")
        except Exception as e:
            log(f"  LUAD / {backend}: FAILED — {e}")
            traceback.print_exc()

    for backend in BACKENDS:
        log(f"  LUAD / {backend}: diagnostics...")
        try:
            diag = collect_diagnostics(x_luad_var, backend, "LUAD", k=10, iters=5, n_boot=5)
            all_diagnostics.append(diag)
            log(f"  LUAD / {backend}: lifetime={diag['mean_lifetime']:.4f}, support={diag['support_size']}, gini={diag['tip_gini']:.3f}")
        except Exception as e:
            log(f"  LUAD / {backend}: diagnostics FAILED — {e}")

    # LUAD ablations
    log("  LUAD: backend sensitivity ablations...")
    for backend in BACKENDS:
        for lam in [1e-6, 1e-5]:
            iters_use = 5 if backend in ("ricci", "diffusion", "phate_like") else 0
            try:
                params = FeatureParams(
                    mode=backend, solver="elastic", k=10, iters=iters_use,
                    alpha=0.5, lambda_=lam, top_k=20, normalize="std",
                    n_boot=1, subsample_frac=0.8, rescale_median=1.0,
                )
                df = evaluate_ccspr_only(x_luad, y_luad, params=params,
                                          var_genes=5000, n_splits=2, test_size=0.25,
                                          seed=SEED + 100, cache_dir=CACHE_DIR)
                all_ablation.append({
                    "dataset": "LUAD", "backend": backend,
                    "ablation": "lambda", "value": lam,
                    "mean_f1": float(df["f1_weighted"].mean()),
                    "std_f1": float(df["f1_weighted"].std()) if len(df) > 1 else 0.0,
                    "n": len(df),
                })
            except Exception as e:
                log(f"  LUAD / {backend} / lambda={lam}: FAILED — {e}")

    save_checkpoint("luad_done")

    # ===================================================================
    # 2. GSE161711 — n=96, n_boot=5, 5-fold CV
    # ===================================================================
    log("=== PHASE 2: GSE161711 Geometry-Family Benchmark ===")
    ds_gse = load_cll_venetoclax()
    x_gse, y_gse = ds_gse["X"], ds_gse["y"]
    log(f"GSE161711 loaded: {x_gse.shape[0]} samples, {x_gse.shape[1]} features")

    x_gse_var, _ = select_top_variable_genes(x_gse, top_n=5000)

    for backend in BACKENDS:
        log(f"  GSE161711 / {backend}: evaluating...")
        try:
            result = evaluate_backend(
                x_gse, y_gse, backend, "GSE161711",
                var_genes=5000, n_splits=5, test_size=0.25, n_boot=5,
                k=10, iters=5, lambda_=1e-6, top_k=20,
            )
            all_benchmark.append(result)
            timing[f"GSE161711_{backend}"] = result["elapsed_sec"]
            log(f"  GSE161711 / {backend}: F1={result['mean_f1']:.4f} ± {result['std_f1']:.4f} ({result['elapsed_sec']:.1f}s)")
        except Exception as e:
            log(f"  GSE161711 / {backend}: FAILED — {e}")
            traceback.print_exc()

    for backend in BACKENDS:
        log(f"  GSE161711 / {backend}: diagnostics...")
        try:
            diag = collect_diagnostics(x_gse_var, backend, "GSE161711", k=10, iters=5, n_boot=5)
            all_diagnostics.append(diag)
            log(f"  GSE161711 / {backend}: lifetime={diag['mean_lifetime']:.4f}, support={diag['support_size']}")
        except Exception as e:
            log(f"  GSE161711 / {backend}: diagnostics FAILED — {e}")

    # GSE161711 ablations
    log("  GSE161711: backend sensitivity ablations...")
    for backend in BACKENDS:
        for lam in [1e-6, 1e-5]:
            iters_use = 5 if backend in ("ricci", "diffusion", "phate_like") else 0
            try:
                params = FeatureParams(
                    mode=backend, solver="elastic", k=10, iters=iters_use,
                    alpha=0.5, lambda_=lam, top_k=20, normalize="std",
                    n_boot=1, subsample_frac=0.8, rescale_median=1.0,
                )
                df = evaluate_ccspr_only(x_gse, y_gse, params=params,
                                          var_genes=5000, n_splits=5, test_size=0.25,
                                          seed=SEED + 100, cache_dir=CACHE_DIR)
                all_ablation.append({
                    "dataset": "GSE161711", "backend": backend,
                    "ablation": "lambda", "value": lam,
                    "mean_f1": float(df["f1_weighted"].mean()),
                    "std_f1": float(df["f1_weighted"].std()) if len(df) > 1 else 0.0,
                    "n": len(df),
                })
            except Exception as e:
                log(f"  GSE161711 / {backend} / lambda={lam}: FAILED — {e}")

    save_checkpoint("gse161711_done")

    # ===================================================================
    # 3. GSE165087 — bounded 320 cells, n_boot=1
    # ===================================================================
    log("=== PHASE 3: GSE165087 Bounded ===")
    scrna_available = False
    try:
        from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
        ds_sc = load_cll_rs_scrna(
            data_path="data/cll_rs_scrna/data.h5ad",
            n_hvg=800, n_pcs=20, max_cells=320, seed=SEED,
        )
        x_sc, y_sc = ds_sc["X"], ds_sc["y"]
        log(f"GSE165087 loaded: {x_sc.shape[0]} samples, {x_sc.shape[1]} features")

        for backend in BACKENDS:
            log(f"  GSE165087 / {backend}: evaluating...")
            try:
                result = evaluate_backend(
                    x_sc, y_sc, backend, "GSE165087",
                    var_genes=min(800, x_sc.shape[1]),
                    n_splits=1, test_size=0.25, n_boot=1,
                    k=10, iters=2, lambda_=1e-6, top_k=20,
                )
                all_benchmark.append(result)
                timing[f"GSE165087_{backend}"] = result["elapsed_sec"]
                log(f"  GSE165087 / {backend}: F1={result['mean_f1']:.4f} ({result['elapsed_sec']:.1f}s)")
            except Exception as e:
                log(f"  GSE165087 / {backend}: FAILED — {e}")
                traceback.print_exc()

        # Tiny lambda sensitivity on GSE165087
        for backend in BACKENDS:
            for lam in [1e-6, 1e-5]:
                iters_use = 2 if backend in ("ricci", "diffusion", "phate_like") else 0
                try:
                    params = FeatureParams(
                        mode=backend, solver="elastic", k=10, iters=iters_use,
                        alpha=0.5, lambda_=lam, top_k=20, normalize="std",
                        n_boot=1, subsample_frac=0.8, rescale_median=1.0,
                    )
                    df = evaluate_ccspr_only(x_sc, y_sc, params=params,
                                              var_genes=min(800, x_sc.shape[1]),
                                              n_splits=1, test_size=0.25,
                                              seed=SEED + 400, cache_dir=CACHE_DIR)
                    all_ablation.append({
                        "dataset": "GSE165087", "backend": backend,
                        "ablation": "lambda", "value": lam,
                        "mean_f1": float(df["f1_weighted"].mean()),
                        "std_f1": 0.0, "n": 1,
                    })
                except Exception as e:
                    log(f"  GSE165087 / {backend} / lambda={lam}: FAILED — {e}")

        scrna_available = True
    except Exception as e:
        log(f"  GSE165087: SKIPPED — {e}")

    save_checkpoint("scrna_done")

    # ===================================================================
    # 4. Save all results
    # ===================================================================
    log("=== PHASE 4: Saving Results ===")

    bench_rows = [{
        "dataset": r["dataset"], "backend": r["backend"],
        "mean_f1": round(r["mean_f1"], 4), "std_f1": round(r["std_f1"], 4),
        "n_splits": r["n_splits"], "elapsed_sec": r["elapsed_sec"],
    } for r in all_benchmark]
    df_bench = pd.DataFrame(bench_rows)
    df_bench.to_csv(TABLES_DIR / "geometry_family_benchmark.csv", index=False)
    log(f"Saved: {TABLES_DIR / 'geometry_family_benchmark.csv'}")

    diag_rows = [{
        "dataset": d["dataset"], "backend": d["backend"],
        "n_ok": d["n_ok"],
        "mean_lifetime": round(d["mean_lifetime"], 4),
        "std_lifetime": round(d["std_lifetime"], 4),
        "mean_prominence": round(d["mean_prominence"], 4),
        "std_prominence": round(d["std_prominence"], 4),
        "support_size": d["support_size"],
        "tip_entropy": round(d["tip_entropy"], 4),
        "tip_gini": round(d["tip_gini"], 4),
    } for d in all_diagnostics]
    df_diag = pd.DataFrame(diag_rows)
    df_diag.to_csv(TABLES_DIR / "geometry_family_diagnostics.csv", index=False)
    log(f"Saved: {TABLES_DIR / 'geometry_family_diagnostics.csv'}")

    df_abl = pd.DataFrame(all_ablation)
    if not df_abl.empty:
        df_abl.to_csv(TABLES_DIR / "geometry_family_ablation.csv", index=False)
    log(f"Saved ablation table")

    pd.DataFrame([{"key": k, "seconds": v} for k, v in timing.items()]).to_csv(
        TABLES_DIR / "geometry_family_timing.csv", index=False)

    full_results = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "benchmark": all_benchmark,
        "diagnostics": all_diagnostics,
        "ablation": all_ablation,
        "timing": timing,
        "scrna_available": scrna_available,
    }
    with open(RESULTS_DIR / "geometry_family_results.json", "w") as f:
        json.dump(full_results, f, indent=2, default=str)

    log("=== SUMMARY ===")
    print("\nGeometry-Family Benchmark:")
    print(df_bench.to_string(index=False))
    if not df_diag.empty:
        print("\nDiagnostics:")
        print(df_diag.to_string(index=False))
    log("=== DONE ===")


if __name__ == "__main__":
    main()
