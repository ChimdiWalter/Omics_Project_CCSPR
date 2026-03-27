#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ccspr.datasets.tcga_luad import download_tcga_luad, load_tcga_luad
from ccspr.eval.ablations import run_luad_ablations
from ccspr.eval.classification import FeatureParams, evaluate_luad_protocol
from ccspr.plots.figures import (
    plot_f1_bars_with_ci,
    plot_lambda_robustness,
    plot_lifetime_prominence,
    plot_tip_eu_vs_ricci,
)
from ccspr.utils.io import ensure_dir, read_yaml, save_json


def _stratified_subsample(x: np.ndarray, y: np.ndarray, n_max: int, seed: int):
    if x.shape[0] <= n_max:
        return x, y
    rng = np.random.default_rng(seed)
    keep = []
    classes, counts = np.unique(y, return_counts=True)
    for c, n in zip(classes, counts):
        idx = np.where(y == c)[0]
        take = max(2, int(round(len(idx) * (n_max / len(y)))))
        take = min(take, len(idx))
        keep.extend(rng.choice(idx, take, replace=False).tolist())
    keep = np.array(sorted(set(keep)), dtype=int)
    if len(keep) > n_max:
        keep = rng.choice(keep, n_max, replace=False)
    return x[keep], y[keep]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/tcga_luad.yaml")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    seed = int(cfg["seed"])
    np.random.seed(seed)

    results_dir = Path(cfg["output"]["results_dir"])
    figures_dir = results_dir / "figures"
    tables_dir = results_dir / "tables"
    ensure_dir(results_dir)
    ensure_dir(figures_dir)
    ensure_dir(tables_dir)
    ensure_dir(Path(cfg["cache_dir"]))

    save_json(
        {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "config": cfg,
        },
        results_dir / "run.json",
    )

    if cfg["dataset"].get("download", True):
        download_tcga_luad(cfg["dataset"].get("root", "data/tcga_luad"))

    ds = load_tcga_luad(
        root=cfg["dataset"].get("root", "data/tcga_luad"),
        min_class_size=int(cfg["dataset"].get("min_class_size", 10)),
    )

    x = ds["X"]
    y = ds["y"]

    n_max = int(cfg["dataset"].get("max_samples", 0))
    if n_max > 0:
        x, y = _stratified_subsample(x, y, n_max=n_max, seed=seed)

    harmonic_params = FeatureParams(**cfg["protocol"]["harmonic_params"])
    ccspr_params = FeatureParams(**cfg["protocol"]["ccspr_params"])

    metrics_main, aux = evaluate_luad_protocol(
        x,
        y,
        var_genes=int(cfg["protocol"]["var_genes"]),
        n_splits=int(cfg["protocol"]["n_splits"]),
        test_size=float(cfg["protocol"]["test_size"]),
        seed=seed,
        cache_dir=cfg["cache_dir"],
        harmonic_params=harmonic_params,
        ccspr_params=ccspr_params,
    )
    metrics_main["phase"] = "main"

    ablation = run_luad_ablations(
        x,
        y,
        base_params=ccspr_params,
        var_genes=int(cfg["ablation"]["var_genes"]),
        n_splits=int(cfg["ablation"]["n_splits"]),
        test_size=float(cfg["ablation"]["test_size"]),
        seed=seed + 99,
        grids=cfg["ablation"].get("grids"),
        cache_dir=cfg["cache_dir"],
    )
    if not ablation.empty:
        ablation["phase"] = "ablation"

    metrics = pd.concat([metrics_main, ablation], ignore_index=True)
    metrics.to_csv(results_dir / "metrics.csv", index=False)

    main_summary = (
        metrics_main.groupby("model")["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    main_summary["ci95"] = 1.96 * main_summary["std"].fillna(0.0) / np.sqrt(main_summary["n"].clip(lower=1))
    main_summary.to_csv(results_dir / "main_results.csv", index=False)
    main_summary.to_csv(tables_dir / "main_results.csv", index=False)

    if not ablation.empty:
        ablation_summary = (
            ablation.groupby(["ablation", "ablation_value"])["f1_weighted"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={"count": "n"})
        )
        ablation_summary.to_csv(results_dir / "ablation_summary.csv", index=False)
        ablation_summary.to_csv(tables_dir / "ablation_summary.csv", index=False)
    else:
        pd.DataFrame(columns=["ablation", "ablation_value", "mean", "std", "n"]).to_csv(
            results_dir / "ablation_summary.csv", index=False
        )
        pd.DataFrame(columns=["ablation", "ablation_value", "mean", "std", "n"]).to_csv(
            tables_dir / "ablation_summary.csv", index=False
        )

    plot_tip_eu_vs_ricci(aux["tip_eu"], aux["tip_ricci"], str(figures_dir / "tip_eu_vs_ricci.png"))
    plot_lifetime_prominence(
        aux["lifetime_eu"],
        aux["lifetime_ricci"],
        aux["prominence_eu"],
        aux["prominence_ricci"],
        str(figures_dir / "h1_lifetime_eu_vs_ricci.png"),
        str(figures_dir / "h1_prominence_eu_vs_ricci.png"),
    )
    plot_f1_bars_with_ci(metrics_main, str(figures_dir / "f1_bars_ci.png"))
    if not ablation.empty:
        plot_lambda_robustness(ablation, str(figures_dir / "lambda_robustness.png"))

    print("Saved metrics:", results_dir / "metrics.csv")
    print("Saved figures dir:", figures_dir)


if __name__ == "__main__":
    main()
