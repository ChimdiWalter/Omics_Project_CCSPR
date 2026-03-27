#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, wilcoxon

from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.eval.classification import FeatureParams, evaluate_ccspr_only, evaluate_luad_protocol
from ccspr.utils.io import ensure_dir, read_yaml


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/tcga_luad.yaml")
    ap.add_argument("--seeds", default="42,52,62,72,82")
    ap.add_argument("--n-splits", type=int, default=3)
    ap.add_argument("--out-prefix", default="luad_multiseed")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]

    out_dir = Path(cfg["output"]["results_dir"])
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    ensure_dir(out_dir)
    ensure_dir(fig_dir)
    ensure_dir(tab_dir)

    ds = load_tcga_luad(
        root=cfg["dataset"].get("root", "data/tcga_luad"),
        min_class_size=int(cfg["dataset"].get("min_class_size", 10)),
    )
    x, y = ds["X"], ds["y"]
    x, y = _stratified_subsample(x, y, n_max=int(cfg["dataset"].get("max_samples", 0)), seed=int(cfg.get("seed", 42)))

    harmonic_params = FeatureParams(**cfg["protocol"]["harmonic_params"])
    ricci_sparse_params = FeatureParams(**cfg["protocol"]["ccspr_params"])
    eu_sparse_params = replace(ricci_sparse_params, mode="euclidean", iters=0)

    main_rows = []
    paired_sparse_rows = []

    for sd in seeds:
        m_main, _ = evaluate_luad_protocol(
            x,
            y,
            var_genes=int(cfg["protocol"].get("var_genes", 5000)),
            n_splits=int(args.n_splits),
            test_size=float(cfg["protocol"].get("test_size", 0.25)),
            seed=sd,
            cache_dir=cfg.get("cache_dir", "data/cache"),
            harmonic_params=harmonic_params,
            ccspr_params=ricci_sparse_params,
        )
        m_main["seed"] = sd
        main_rows.append(m_main)

        m_ricci = evaluate_ccspr_only(
            x,
            y,
            params=ricci_sparse_params,
            var_genes=int(cfg["protocol"].get("var_genes", 5000)),
            n_splits=int(args.n_splits),
            test_size=float(cfg["protocol"].get("test_size", 0.25)),
            seed=sd + 7000,
            cache_dir=cfg.get("cache_dir", "data/cache"),
        )
        m_eu = evaluate_ccspr_only(
            x,
            y,
            params=eu_sparse_params,
            var_genes=int(cfg["protocol"].get("var_genes", 5000)),
            n_splits=int(args.n_splits),
            test_size=float(cfg["protocol"].get("test_size", 0.25)),
            seed=sd + 7000,
            cache_dir=cfg.get("cache_dir", "data/cache"),
        )
        merged = m_ricci[["split", "f1_weighted"]].rename(columns={"f1_weighted": "f1_ricci_sparse"}).merge(
            m_eu[["split", "f1_weighted"]].rename(columns={"f1_weighted": "f1_eu_sparse"}),
            on="split",
            how="inner",
        )
        merged["seed"] = sd
        paired_sparse_rows.append(merged)

    df_main = pd.concat(main_rows, ignore_index=True)
    df_main.to_csv(out_dir / f"metrics_{args.out_prefix}.csv", index=False)

    df_paired_sparse = pd.concat(paired_sparse_rows, ignore_index=True)
    df_paired_sparse["delta_ricci_minus_eu_sparse"] = (
        df_paired_sparse["f1_ricci_sparse"] - df_paired_sparse["f1_eu_sparse"]
    )
    df_paired_sparse.to_csv(out_dir / f"paired_{args.out_prefix}_ricci_vs_eu_sparse.csv", index=False)

    pivot = df_main.pivot_table(index=["seed", "split"], columns="model", values="f1_weighted", aggfunc="first").reset_index()

    stats_rows = []
    if {"ccspr", "harmonic"}.issubset(set(pivot.columns)):
        a = pivot["ccspr"].to_numpy(dtype=float)
        b = pivot["harmonic"].to_numpy(dtype=float)
        stats_rows.append(_paired_stats(a, b, "ccspr", "harmonic"))

    if {"ccspr", "standard"}.issubset(set(pivot.columns)):
        a = pivot["ccspr"].to_numpy(dtype=float)
        b = pivot["standard"].to_numpy(dtype=float)
        stats_rows.append(_paired_stats(a, b, "ccspr", "standard"))

    stats_rows.append(
        _paired_stats(
            df_paired_sparse["f1_ricci_sparse"].to_numpy(dtype=float),
            df_paired_sparse["f1_eu_sparse"].to_numpy(dtype=float),
            "ricci_sparse",
            "euclidean_sparse",
        )
    )

    df_stats = pd.DataFrame(stats_rows)
    df_stats.to_csv(out_dir / f"paired_stats_{args.out_prefix}.csv", index=False)
    df_stats.to_csv(tab_dir / f"paired_stats_{args.out_prefix}.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    vals = [
        pivot["ccspr"].to_numpy(dtype=float) - pivot["harmonic"].to_numpy(dtype=float)
        if {"ccspr", "harmonic"}.issubset(set(pivot.columns))
        else np.array([]),
        df_paired_sparse["delta_ricci_minus_eu_sparse"].to_numpy(dtype=float),
    ]
    labels = ["ccspr-harmonic", "ricci_sparse-eu_sparse"]
    ax.boxplot([v for v in vals if v.size > 0], labels=[labels[i] for i, v in enumerate(vals) if v.size > 0])
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)
    ax.set_ylabel("Paired F1 difference")
    ax.set_title("LUAD multiseed paired deltas")
    fig.tight_layout()
    fig.savefig(fig_dir / f"paired_deltas_{args.out_prefix}.png", dpi=180)
    plt.close(fig)

    print("Saved:", out_dir / f"metrics_{args.out_prefix}.csv")
    print("Saved:", out_dir / f"paired_stats_{args.out_prefix}.csv")


if __name__ == "__main__":
    main()
