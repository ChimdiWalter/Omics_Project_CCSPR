#!/usr/bin/env python3
"""Lightweight experiments to strengthen the CC-SPR manuscript.

Priority order:
  1. LUAD multiseed paired analysis (reuses cached outputs)
  2. GSE161711 bounded ablation sweep (small safe grid)
  3. Diagnostics from existing outputs (support size, entropy, Gini, TIP, H1)
  4. Optional tiny scRNA sensitivity (k ablation, bounded)

No heavy scRNA configs are launched. All outputs skip if results already exist.
"""
from __future__ import annotations

import argparse
import resource
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import entropy as scipy_entropy
from scipy.stats import ttest_rel, wilcoxon
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedShuffleSplit

from ccspr.datasets.cll_venetoclax import load_cll_venetoclax
from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.eval.ablations import run_luad_ablations
from ccspr.eval.classification import (
    FeatureParams,
    _evaluate_feature_model,
    _fit_predict_logreg,
    evaluate_ccspr_only,
    evaluate_luad_protocol,
)
from ccspr.plots.figures import (
    plot_f1_bars_with_ci,
    plot_lambda_robustness,
    plot_lifetime_prominence,
    plot_tip_eu_vs_ricci,
)
from ccspr.preprocess.basic import select_top_variable_genes, standardize
from ccspr.utils.io import ensure_dir, read_yaml, save_json


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ts() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _rss_gb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / (1024.0 * 1024.0)


def _log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", flush=True)


def _paired_stats(a: np.ndarray, b: np.ndarray, label_a: str, label_b: str) -> dict:
    delta = a - b
    out = {
        "comparison": f"{label_a}_vs_{label_b}",
        "n_pairs": int(delta.size),
        "mean_delta": float(np.mean(delta)),
        "std_delta": float(np.std(delta, ddof=1)) if delta.size > 1 else 0.0,
        "median_delta": float(np.median(delta)),
    }
    if delta.size >= 2:
        t_stat, t_p = ttest_rel(a, b)
        out["paired_t_stat"] = float(t_stat)
        out["paired_t_p"] = float(t_p)
    else:
        out["paired_t_stat"] = np.nan
        out["paired_t_p"] = np.nan

    if delta.size >= 1 and np.any(np.abs(delta) > 0):
        try:
            w_stat, w_p = wilcoxon(a, b, zero_method="wilcox", correction=False)
            out["wilcoxon_stat"] = float(w_stat)
            out["wilcoxon_p"] = float(w_p)
        except Exception:
            out["wilcoxon_stat"] = np.nan
            out["wilcoxon_p"] = np.nan
    else:
        out["wilcoxon_stat"] = np.nan
        out["wilcoxon_p"] = np.nan
    return out


def _stratified_subsample(x: np.ndarray, y: np.ndarray, n_max: int, seed: int):
    if n_max <= 0 or x.shape[0] <= n_max:
        return x, y
    rng = np.random.default_rng(seed)
    keep = []
    classes, _ = np.unique(y, return_counts=True)
    for c in classes:
        idx = np.where(y == c)[0]
        take = max(2, int(round(len(idx) * (n_max / len(y)))))
        take = min(take, len(idx))
        keep.extend(rng.choice(idx, take, replace=False).tolist())
    keep = np.array(sorted(set(keep)), dtype=int)
    if len(keep) > n_max:
        keep = rng.choice(keep, n_max, replace=False)
    return x[keep], y[keep]


def _summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["model", "mean", "std", "n", "ci95"])
    out = (
        df.groupby("model")["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    out["ci95"] = 1.96 * out["std"].fillna(0.0) / np.sqrt(out["n"].clip(lower=1))
    return out


def _ablation_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["ablation", "ablation_value", "mean", "std", "n"])
    return (
        df.groupby(["ablation", "ablation_value"])["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )


def _write_checkpoint(path: Path, phase: str, meta: dict) -> None:
    payload = {"timestamp": _ts(), "phase": phase, "max_rss_gb": round(_rss_gb(), 3)}
    payload.update(meta)
    save_json(payload, path)


# ---------------------------------------------------------------------------
# 1. LUAD multiseed paired analysis
# ---------------------------------------------------------------------------

def run_luad_multiseed(cfg_path: str, out_dir: Path, tables_dir: Path, figures_dir: Path) -> bool:
    """Returns True if work was done, False if skipped."""
    tag = "luad_multiseed_v2"
    final_metrics = out_dir / f"metrics_{tag}.csv"
    final_stats = out_dir / f"paired_stats_{tag}.csv"
    if final_metrics.exists() and final_stats.exists():
        _log(f"LUAD multiseed outputs already exist ({tag}); skipping")
        return False

    cfg = read_yaml(cfg_path)
    seeds = [42, 52, 62, 72, 82, 92]
    n_splits = 3

    ds = load_tcga_luad(
        root=cfg["dataset"].get("root", "data/tcga_luad"),
        min_class_size=int(cfg["dataset"].get("min_class_size", 10)),
    )
    x, y = ds["X"], ds["y"]
    x, y = _stratified_subsample(
        x, y,
        n_max=int(cfg["dataset"].get("max_samples", 0)),
        seed=int(cfg.get("seed", 42)),
    )
    _log(f"LUAD dataset: samples={x.shape[0]} features={x.shape[1]}")

    harmonic_params = FeatureParams(**cfg["protocol"]["harmonic_params"])
    ccspr_params = FeatureParams(**cfg["protocol"]["ccspr_params"])

    from dataclasses import replace
    eu_sparse_params = replace(ccspr_params, mode="euclidean", iters=0)

    main_rows = []
    paired_sparse_rows = []

    for sd in seeds:
        _log(f"LUAD multiseed: seed={sd}")
        m_main, _ = evaluate_luad_protocol(
            x, y,
            var_genes=int(cfg["protocol"].get("var_genes", 5000)),
            n_splits=n_splits,
            test_size=float(cfg["protocol"].get("test_size", 0.25)),
            seed=sd,
            cache_dir=cfg.get("cache_dir", "data/cache"),
            harmonic_params=harmonic_params,
            ccspr_params=ccspr_params,
        )
        m_main["seed"] = sd
        main_rows.append(m_main)

        m_ricci = evaluate_ccspr_only(
            x, y,
            params=ccspr_params,
            var_genes=int(cfg["protocol"].get("var_genes", 5000)),
            n_splits=n_splits,
            test_size=float(cfg["protocol"].get("test_size", 0.25)),
            seed=sd + 7000,
            cache_dir=cfg.get("cache_dir", "data/cache"),
        )
        m_eu = evaluate_ccspr_only(
            x, y,
            params=eu_sparse_params,
            var_genes=int(cfg["protocol"].get("var_genes", 5000)),
            n_splits=n_splits,
            test_size=float(cfg["protocol"].get("test_size", 0.25)),
            seed=sd + 7000,
            cache_dir=cfg.get("cache_dir", "data/cache"),
        )
        merged = (
            m_ricci[["split", "f1_weighted"]]
            .rename(columns={"f1_weighted": "f1_ricci_sparse"})
            .merge(
                m_eu[["split", "f1_weighted"]].rename(columns={"f1_weighted": "f1_eu_sparse"}),
                on="split",
                how="inner",
            )
        )
        merged["seed"] = sd
        paired_sparse_rows.append(merged)
        _log(f"LUAD seed={sd} done; max_rss_gb={_rss_gb():.2f}")

    df_main = pd.concat(main_rows, ignore_index=True)
    df_main.to_csv(final_metrics, index=False)

    df_paired = pd.concat(paired_sparse_rows, ignore_index=True)
    df_paired["delta_ricci_minus_eu_sparse"] = (
        df_paired["f1_ricci_sparse"] - df_paired["f1_eu_sparse"]
    )
    df_paired.to_csv(out_dir / f"paired_{tag}_ricci_vs_eu_sparse.csv", index=False)

    pivot = (
        df_main.pivot_table(index=["seed", "split"], columns="model", values="f1_weighted", aggfunc="first")
        .reset_index()
    )

    stats_rows = []
    for la, lb in [("ccspr", "harmonic"), ("ccspr", "standard")]:
        if {la, lb}.issubset(set(pivot.columns)):
            stats_rows.append(_paired_stats(
                pivot[la].to_numpy(dtype=float),
                pivot[lb].to_numpy(dtype=float),
                la, lb,
            ))
    stats_rows.append(_paired_stats(
        df_paired["f1_ricci_sparse"].to_numpy(dtype=float),
        df_paired["f1_eu_sparse"].to_numpy(dtype=float),
        "ricci_sparse", "euclidean_sparse",
    ))

    df_stats = pd.DataFrame(stats_rows)
    df_stats.to_csv(final_stats, index=False)
    df_stats.to_csv(tables_dir / f"paired_stats_{tag}.csv", index=False)

    # Paired deltas boxplot
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    vals, labels = [], []
    if {"ccspr", "harmonic"}.issubset(set(pivot.columns)):
        vals.append(pivot["ccspr"].to_numpy(float) - pivot["harmonic"].to_numpy(float))
        labels.append("ccspr-harmonic")
    vals.append(df_paired["delta_ricci_minus_eu_sparse"].to_numpy(float))
    labels.append("ricci_sparse-eu_sparse")
    ax.boxplot([v for v in vals if v.size > 0], labels=[labels[i] for i, v in enumerate(vals) if v.size > 0])
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)
    ax.set_ylabel("Paired F1 difference")
    ax.set_title("LUAD multiseed paired deltas (v2)")
    fig.tight_layout()
    fig.savefig(figures_dir / f"paired_deltas_{tag}.png", dpi=180)
    plt.close(fig)

    _log(f"LUAD multiseed complete: {final_metrics}")
    return True


# ---------------------------------------------------------------------------
# 2. GSE161711 bounded ablation sweep
# ---------------------------------------------------------------------------

def run_gse161711_bounded_ablation(cfg_path: str, out_dir: Path, tables_dir: Path, figures_dir: Path) -> bool:
    """Returns True if work was done."""
    cfg = read_yaml(cfg_path)
    tag = str(cfg["output"].get("tag", "cll_venetoclax_gse161711_bounded"))
    final_metrics = out_dir / f"metrics_{tag}.csv"
    final_ablation = out_dir / f"ablation_summary_{tag}.csv"

    if bool(cfg.get("skip_if_complete", True)) and final_metrics.exists() and final_ablation.exists():
        _log(f"GSE161711 bounded ablation outputs already exist ({tag}); skipping")
        return False

    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)

    ds = load_cll_venetoclax(
        processed_matrix_path=cfg["dataset"].get("processed_matrix_path", "data/cll_venetoclax/processed_matrix.tsv"),
        label_col=cfg["dataset"].get("label_col", "label"),
    )
    _log(f"GSE161711 dataset: samples={ds['X'].shape[0]} features={ds['X'].shape[1]}")

    ccspr_params = FeatureParams(**cfg["protocol"]["ccspr_params"])

    ablation = run_luad_ablations(
        ds["X"],
        ds["y"],
        base_params=ccspr_params,
        var_genes=min(int(cfg["ablation"]["var_genes"]), ds["X"].shape[1]),
        n_splits=int(cfg["ablation"]["n_splits"]),
        test_size=float(cfg["ablation"]["test_size"]),
        seed=seed + 99,
        grids=cfg["ablation"].get("grids"),
        cache_dir=cfg.get("cache_dir", "data/cache"),
    )

    if not ablation.empty:
        ablation["phase"] = "ablation"
        ablation["dataset"] = tag

    ablation.to_csv(final_metrics, index=False)
    ablation_summary = _ablation_summary(ablation)
    ablation_summary.to_csv(final_ablation, index=False)
    ablation_summary.to_csv(tables_dir / f"ablation_summary_{tag}.csv", index=False)

    if not ablation.empty and "lambda" in ablation["ablation"].values:
        plot_lambda_robustness(ablation, str(figures_dir / f"lambda_robustness_{tag}.png"))

    _log(f"GSE161711 bounded ablation complete: {final_metrics}")
    return True


# ---------------------------------------------------------------------------
# 3. Diagnostics from existing outputs
# ---------------------------------------------------------------------------

def _gini(x: np.ndarray) -> float:
    """Gini coefficient of a 1-D array."""
    x = np.abs(x.astype(float))
    if x.sum() == 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    idx = np.arange(1, n + 1)
    return float((2 * np.sum(idx * x) - (n + 1) * np.sum(x)) / (n * np.sum(x)))


def run_diagnostics(out_dir: Path, tables_dir: Path, figures_dir: Path) -> bool:
    """Compute support size, entropy, Gini, TIP concentration, H1 lifetime/prominence from existing arrays."""
    tag = "lightweight_diagnostics"
    final_csv = tables_dir / f"{tag}.csv"

    # Collect all figure_arrays npz files
    npz_files = sorted(out_dir.glob("figure_arrays_*.npz"))
    if not npz_files:
        _log("No figure_arrays npz files found; skipping diagnostics")
        return False

    rows = []
    for npz_path in npz_files:
        dataset_tag = npz_path.stem.replace("figure_arrays_", "")
        _log(f"Computing diagnostics for {dataset_tag}")
        data = np.load(npz_path, allow_pickle=True)

        for backend, tip_key, life_key, prom_key in [
            ("euclidean", "tip_eu", "lifetime_eu", "prominence_eu"),
            ("ricci", "tip_ricci", "lifetime_ricci", "prominence_ricci"),
        ]:
            tip = data.get(tip_key, np.array([], dtype=float))
            lifetime = data.get(life_key, np.array([], dtype=float))
            prominence = data.get(prom_key, np.array([], dtype=float))

            if tip.size == 0:
                continue

            support = int(np.count_nonzero(tip > 0))
            tip_norm = tip / tip.sum() if tip.sum() > 0 else tip
            tip_entropy = float(scipy_entropy(tip_norm + 1e-12))
            tip_gini = _gini(tip)
            tip_top5 = float(np.sort(tip)[-5:].sum() / max(tip.sum(), 1e-12))
            tip_top10 = float(np.sort(tip)[-10:].sum() / max(tip.sum(), 1e-12))

            row = {
                "dataset": dataset_tag,
                "backend": backend,
                "support_size": support,
                "total_features": int(tip.size),
                "support_fraction": round(support / max(tip.size, 1), 4),
                "tip_entropy": round(tip_entropy, 4),
                "tip_gini": round(tip_gini, 4),
                "tip_top5_concentration": round(tip_top5, 4),
                "tip_top10_concentration": round(tip_top10, 4),
            }

            if lifetime.size > 0:
                row["h1_lifetime_mean"] = round(float(np.mean(lifetime)), 4)
                row["h1_lifetime_std"] = round(float(np.std(lifetime)), 4)
                row["h1_lifetime_median"] = round(float(np.median(lifetime)), 4)
            if prominence.size > 0:
                row["h1_prominence_mean"] = round(float(np.mean(prominence)), 4)
                row["h1_prominence_std"] = round(float(np.std(prominence)), 4)
                row["h1_prominence_median"] = round(float(np.median(prominence)), 4)

            rows.append(row)

    # Also pull diagnostics from LUAD/GSE161711 figure arrays if they exist in results/figures
    for npz_path in sorted(out_dir.glob("*.npz")):
        if npz_path.stem.startswith("figure_arrays_"):
            continue  # already processed above
        # skip non-figure arrays
        pass

    df = pd.DataFrame(rows)
    df.to_csv(final_csv, index=False)
    df.to_csv(out_dir / f"{tag}.csv", index=False)

    # Runtime and memory summary from checkpoint files
    runtime_rows = []
    for ckpt_path in sorted(out_dir.glob("checkpoint_*.json")):
        import json
        with open(ckpt_path) as f:
            ckpt = json.load(f)
        runtime_rows.append({
            "checkpoint_file": ckpt_path.name,
            "phase": ckpt.get("phase", "unknown"),
            "max_rss_gb": ckpt.get("max_rss_gb", np.nan),
            "timestamp": ckpt.get("timestamp", ""),
        })
    if runtime_rows:
        df_runtime = pd.DataFrame(runtime_rows)
        df_runtime.to_csv(tables_dir / "runtime_memory_summary.csv", index=False)
        _log(f"Runtime/memory summary: {len(runtime_rows)} checkpoints")

    # Diagnostics figure: support size + Gini bar chart
    if not df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

        # Support fraction
        ax = axes[0]
        labels = [f"{r['dataset']}\n({r['backend']})" for _, r in df.iterrows()]
        ax.barh(labels, df["support_fraction"])
        ax.set_xlabel("Support fraction")
        ax.set_title("TIP support fraction")

        # Gini
        ax = axes[1]
        ax.barh(labels, df["tip_gini"])
        ax.set_xlabel("Gini coefficient")
        ax.set_title("TIP Gini coefficient")

        # Top-10 concentration
        ax = axes[2]
        ax.barh(labels, df["tip_top10_concentration"])
        ax.set_xlabel("Top-10 concentration")
        ax.set_title("TIP top-10 concentration")

        fig.tight_layout()
        fig.savefig(figures_dir / "diagnostics_tip_summary.png", dpi=180)
        plt.close(fig)

    _log(f"Diagnostics complete: {final_csv}")
    return True


# ---------------------------------------------------------------------------
# 4. Optional tiny scRNA sensitivity
# ---------------------------------------------------------------------------

def run_scrna_tiny_sensitivity(cfg_path: str, out_dir: Path, tables_dir: Path, figures_dir: Path) -> bool:
    """Tiny k-ablation on bounded scRNA. Skips if outputs exist."""
    if not Path(cfg_path).exists():
        _log(f"scRNA tiny sensitivity config not found: {cfg_path}; skipping")
        return False

    cfg = read_yaml(cfg_path)
    tag = str(cfg["output"].get("tag", "cll_rs_scrna_gse165087_tiny_sensitivity"))

    final_metrics = out_dir / f"metrics_{tag}.csv"
    final_main = out_dir / f"main_results_{tag}.csv"
    final_ablation = out_dir / f"ablation_summary_{tag}.csv"

    if bool(cfg.get("skip_if_complete", True)) and final_metrics.exists() and final_main.exists() and final_ablation.exists():
        _log(f"scRNA tiny sensitivity outputs already exist ({tag}); skipping")
        return False

    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)
    checkpoint_json = out_dir / f"checkpoint_{tag}.json"

    _log("Loading bounded scRNA dataset for tiny sensitivity")
    from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
    ds = load_cll_rs_scrna(
        data_path=cfg["dataset"]["data_path"],
        label_key=cfg["dataset"].get("label_key", "label"),
        min_genes=int(cfg["dataset"].get("min_genes", 200)),
        min_cells=int(cfg["dataset"].get("min_cells", 3)),
        n_hvg=int(cfg["dataset"].get("n_hvg", 1000)),
        n_pcs=int(cfg["dataset"].get("n_pcs", 25)),
        max_cells=cfg["dataset"].get("max_cells", 400),
        seed=seed,
        use_preprocessed_if_available=bool(cfg["dataset"].get("use_preprocessed_if_available", True)),
    )
    _write_checkpoint(checkpoint_json, "dataset_loaded", {"n_cells": int(ds["X"].shape[0]), "n_features": int(ds["X"].shape[1])})
    _log(f"scRNA dataset ready: cells={ds['X'].shape[0]} features={ds['X'].shape[1]} max_rss_gb={_rss_gb():.2f}")

    ccspr_params = FeatureParams(**cfg["protocol"]["ccspr_params"])
    harmonic_params = FeatureParams(**cfg["protocol"]["harmonic_params"])
    var_genes = min(int(cfg["protocol"].get("var_genes", ds["X"].shape[1])), ds["X"].shape[1])

    splitter = StratifiedShuffleSplit(
        n_splits=int(cfg["protocol"].get("n_splits", 1)),
        test_size=float(cfg["protocol"].get("test_size", 0.25)),
        random_state=seed,
    )
    itr, ite = next(splitter.split(ds["X"], ds["y"]))
    x_train, x_test = ds["X"][itr], ds["X"][ite]
    y_train, y_test = ds["y"][itr], ds["y"][ite]
    x_train_var, idx_var = select_top_variable_genes(x_train, top_n=var_genes)
    x_test_var = x_test[:, idx_var]

    metrics_rows = []

    # Standard baseline
    _log("scRNA tiny: standard baseline")
    xtr_std, xts_std = standardize(x_train_var, x_test_var)
    y_pred_std = _fit_predict_logreg(xtr_std, y_train, xts_std, seed)
    metrics_rows.append({"split": 0, "model": "standard", "f1_weighted": float(f1_score(y_test, y_pred_std, average="weighted")), "distance_mode": "euclidean", "solver": "none"})

    # Harmonic
    _log("scRNA tiny: harmonic baseline")
    f1_h, _ = _evaluate_feature_model(x_train_var, y_train, x_test_var, y_test, harmonic_params, seed + 1000, cfg.get("cache_dir", "data/cache"))
    metrics_rows.append({"split": 0, "model": "harmonic", "f1_weighted": float(f1_h), "distance_mode": harmonic_params.mode, "solver": harmonic_params.solver, "k": harmonic_params.k})

    # CC-SPR main
    _log("scRNA tiny: CC-SPR main")
    f1_c, _ = _evaluate_feature_model(x_train_var, y_train, x_test_var, y_test, ccspr_params, seed + 2000, cfg.get("cache_dir", "data/cache"))
    metrics_rows.append({"split": 0, "model": "ccspr", "f1_weighted": float(f1_c), "distance_mode": ccspr_params.mode, "solver": ccspr_params.solver, "k": ccspr_params.k})

    metrics_main = pd.DataFrame(metrics_rows)

    # Ablation
    ablation_rows = []
    if bool(cfg.get("ablation", {}).get("enabled", False)):
        kind = str(cfg["ablation"].get("kind", "k"))
        values = list(cfg["ablation"].get("values", []))
        for idx, value in enumerate(values):
            if kind == "lambda":
                params = FeatureParams(**{**cfg["protocol"]["ccspr_params"], "lambda_": float(value)})
            else:
                params = FeatureParams(**{**cfg["protocol"]["ccspr_params"], "k": int(value)})
            _log(f"scRNA tiny ablation {kind}={value}")
            f1_a, _ = _evaluate_feature_model(x_train_var, y_train, x_test_var, y_test, params, seed + 3000 + idx, cfg.get("cache_dir", "data/cache"))
            ablation_rows.append({
                "split": 0, "model": "ccspr", "f1_weighted": float(f1_a),
                "distance_mode": params.mode, "solver": params.solver,
                "k": params.k, "iters": params.iters,
                "lambda": params.lambda_, "top_k": params.top_k,
                "normalize": params.normalize,
                "ablation": kind, "ablation_value": value,
                "phase": "ablation", "dataset": tag,
            })
            _write_checkpoint(checkpoint_json, f"ablation_{kind}_{value}", {"rows_written": len(ablation_rows)})
            _log(f"scRNA tiny ablation {kind}={value} done; max_rss_gb={_rss_gb():.2f}")

    metrics_main["phase"] = "main"
    metrics_main["dataset"] = tag
    metrics_main.to_csv(final_metrics, index=False)

    summary = _summary(metrics_main)
    summary.to_csv(final_main, index=False)
    summary.to_csv(tables_dir / f"main_results_{tag}.csv", index=False)

    ablation_df = pd.DataFrame(ablation_rows)
    ablation_summary = _ablation_summary(ablation_df)
    ablation_summary.to_csv(final_ablation, index=False)
    ablation_summary.to_csv(tables_dir / f"ablation_summary_{tag}.csv", index=False)

    _write_checkpoint(checkpoint_json, "complete", {"final_metrics": str(final_metrics)})
    _log(f"scRNA tiny sensitivity complete: {final_metrics}")
    return True


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Lightweight CC-SPR experiments")
    ap.add_argument("--luad-config", default="configs/tcga_luad.yaml")
    ap.add_argument("--gse161711-config", default="configs/cll_venetoclax_bounded_ablation.yaml")
    ap.add_argument("--scrna-config", default="configs/cll_rs_scrna_tiny_sensitivity.yaml")
    ap.add_argument("--skip-luad", action="store_true", help="Skip LUAD multiseed")
    ap.add_argument("--skip-gse161711", action="store_true", help="Skip GSE161711 bounded ablation")
    ap.add_argument("--skip-diagnostics", action="store_true", help="Skip diagnostics")
    ap.add_argument("--skip-scrna", action="store_true", help="Skip optional scRNA sensitivity")
    args = ap.parse_args()

    out_dir = Path("results")
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    ensure_dir(out_dir)
    ensure_dir(tables_dir)
    ensure_dir(figures_dir)
    ensure_dir(Path("data/cache"))

    master_checkpoint = out_dir / "checkpoint_lightweight_strengthen.json"
    completed = []

    _log("=" * 60)
    _log("CC-SPR Lightweight Manuscript Strengthening")
    _log("=" * 60)

    # Priority 1: LUAD multiseed
    if not args.skip_luad:
        _log("--- Priority 1: LUAD multiseed paired analysis ---")
        try:
            did_work = run_luad_multiseed(args.luad_config, out_dir, tables_dir, figures_dir)
            completed.append({"task": "luad_multiseed", "status": "done" if did_work else "skipped_existing"})
        except Exception as e:
            _log(f"LUAD multiseed FAILED: {e}")
            completed.append({"task": "luad_multiseed", "status": f"failed: {e}"})

    # Priority 2: GSE161711 bounded ablation
    if not args.skip_gse161711:
        _log("--- Priority 2: GSE161711 bounded ablation sweep ---")
        try:
            did_work = run_gse161711_bounded_ablation(args.gse161711_config, out_dir, tables_dir, figures_dir)
            completed.append({"task": "gse161711_bounded_ablation", "status": "done" if did_work else "skipped_existing"})
        except Exception as e:
            _log(f"GSE161711 bounded ablation FAILED: {e}")
            completed.append({"task": "gse161711_bounded_ablation", "status": f"failed: {e}"})

    # Priority 3: Diagnostics
    if not args.skip_diagnostics:
        _log("--- Priority 3: Diagnostics from existing outputs ---")
        try:
            did_work = run_diagnostics(out_dir, tables_dir, figures_dir)
            completed.append({"task": "diagnostics", "status": "done" if did_work else "skipped_no_data"})
        except Exception as e:
            _log(f"Diagnostics FAILED: {e}")
            completed.append({"task": "diagnostics", "status": f"failed: {e}"})

    # Priority 4: Optional scRNA tiny sensitivity
    if not args.skip_scrna:
        _log("--- Priority 4: Optional scRNA tiny sensitivity ---")
        try:
            did_work = run_scrna_tiny_sensitivity(args.scrna_config, out_dir, tables_dir, figures_dir)
            completed.append({"task": "scrna_tiny_sensitivity", "status": "done" if did_work else "skipped_existing"})
        except Exception as e:
            _log(f"scRNA tiny sensitivity FAILED: {e}")
            completed.append({"task": "scrna_tiny_sensitivity", "status": f"failed: {e}"})

    _write_checkpoint(master_checkpoint, "all_complete", {"tasks": completed})

    _log("=" * 60)
    _log("All lightweight experiments finished")
    for c in completed:
        _log(f"  {c['task']}: {c['status']}")
    _log(f"  max_rss_gb: {_rss_gb():.2f}")
    _log("=" * 60)


if __name__ == "__main__":
    main()
