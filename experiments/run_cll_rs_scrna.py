#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
from ccspr.eval.ablations import run_luad_ablations
from ccspr.eval.classification import FeatureParams, evaluate_luad_protocol
from ccspr.plots.figures import (
    plot_f1_bars_with_ci,
    plot_lambda_robustness,
    plot_lifetime_prominence,
    plot_tip_eu_vs_ricci,
)
from ccspr.utils.io import ensure_dir, read_yaml, save_json


def _save_table_both(df: pd.DataFrame, out_path: Path, global_path: Path) -> None:
    df.to_csv(out_path, index=False)
    df.to_csv(global_path, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/cll_rs_scrna.yaml")
    args = ap.parse_args()
    cfg = read_yaml(args.config)

    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)

    out = Path(cfg["output"]["results_dir"])
    figures_dir = out / "figures"
    tables_dir = out / "tables"
    global_fig_dir = Path(cfg["output"].get("global_figures_dir", "results/figures"))
    global_tables_dir = Path(cfg["output"].get("global_tables_dir", "results/tables"))
    ensure_dir(out)
    ensure_dir(figures_dir)
    ensure_dir(tables_dir)
    ensure_dir(global_fig_dir)
    ensure_dir(global_tables_dir)

    save_json({"timestamp": datetime.utcnow().isoformat() + "Z", "config": cfg}, out / "run.json")

    ds = load_cll_rs_scrna(
        data_path=cfg["dataset"]["data_path"],
        label_key=cfg["dataset"].get("label_key", "label"),
        min_genes=int(cfg["dataset"].get("min_genes", 200)),
        min_cells=int(cfg["dataset"].get("min_cells", 3)),
        n_hvg=int(cfg["dataset"].get("n_hvg", 2000)),
        n_pcs=int(cfg["dataset"].get("n_pcs", 50)),
        max_cells=cfg["dataset"].get("max_cells", 5000),
        use_precomputed_pca=bool(cfg["dataset"].get("use_precomputed_pca", True)),
        seed=seed,
    )
    pd.DataFrame([ds["meta"]]).to_csv(out / "dataset_meta.csv", index=False)

    harmonic_params = FeatureParams(**cfg["protocol"]["harmonic_params"])
    eu_sparse_params = FeatureParams(**cfg["protocol"]["eu_sparse_params"])
    ccspr_params = FeatureParams(**cfg["protocol"]["ccspr_params"])

    metrics_main, aux = evaluate_luad_protocol(
        ds["X"],
        ds["y"],
        var_genes=min(int(cfg["protocol"]["var_genes"]), ds["X"].shape[1]),
        n_splits=int(cfg["protocol"]["n_splits"]),
        test_size=float(cfg["protocol"]["test_size"]),
        seed=seed,
        cache_dir=cfg["cache_dir"],
        harmonic_params=harmonic_params,
        eu_sparse_params=eu_sparse_params,
        ccspr_params=ccspr_params,
    )
    metrics_main["phase"] = "main"

    ablation = run_luad_ablations(
        ds["X"],
        ds["y"],
        base_params=ccspr_params,
        var_genes=min(int(cfg["ablation"]["var_genes"]), ds["X"].shape[1]),
        n_splits=int(cfg["ablation"]["n_splits"]),
        test_size=float(cfg["ablation"]["test_size"]),
        seed=seed + 99,
        grids=cfg["ablation"].get("grids"),
        cache_dir=cfg["cache_dir"],
    )
    if not ablation.empty:
        ablation["phase"] = "ablation"

    metrics = pd.concat([metrics_main, ablation], ignore_index=True)
    metrics.to_csv(out / "metrics.csv", index=False)
    metrics.to_csv(global_tables_dir / "cll_rs_scrna_metrics.csv", index=False)

    main_summary = (
        metrics_main.groupby("model")["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    main_summary["ci95"] = 1.96 * main_summary["std"].fillna(0.0) / np.sqrt(main_summary["n"].clip(lower=1))
    _save_table_both(
        main_summary,
        out / "main_results.csv",
        global_tables_dir / "cll_rs_scrna_main_results.csv",
    )
    _save_table_both(
        main_summary,
        tables_dir / "main_results.csv",
        global_tables_dir / "cll_rs_scrna_main_results_table.csv",
    )

    if not ablation.empty:
        ablation_summary = (
            ablation.groupby(["ablation", "ablation_value"])["f1_weighted"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={"count": "n"})
        )
    else:
        ablation_summary = pd.DataFrame(columns=["ablation", "ablation_value", "mean", "std", "n"])

    _save_table_both(
        ablation_summary,
        out / "ablation_summary.csv",
        global_tables_dir / "cll_rs_scrna_ablation_summary.csv",
    )
    _save_table_both(
        ablation_summary,
        tables_dir / "ablation_summary.csv",
        global_tables_dir / "cll_rs_scrna_ablation_summary_table.csv",
    )

    fig_tip = figures_dir / "tip_eu_vs_ricci.png"
    fig_life = figures_dir / "h1_lifetime_eu_vs_ricci.png"
    fig_prom = figures_dir / "h1_prominence_eu_vs_ricci.png"
    fig_f1 = figures_dir / "f1_bars_ci.png"
    fig_lambda = figures_dir / "lambda_robustness.png"

    plot_tip_eu_vs_ricci(aux["tip_eu"], aux["tip_ricci"], str(fig_tip))
    plot_lifetime_prominence(
        aux["lifetime_eu"],
        aux["lifetime_ricci"],
        aux["prominence_eu"],
        aux["prominence_ricci"],
        str(fig_life),
        str(fig_prom),
    )
    plot_f1_bars_with_ci(metrics_main, str(fig_f1))
    if not ablation.empty:
        plot_lambda_robustness(ablation, str(fig_lambda))

    plot_tip_eu_vs_ricci(aux["tip_eu"], aux["tip_ricci"], str(global_fig_dir / "cll_rs_scrna_tip_eu_vs_ricci.png"))
    plot_lifetime_prominence(
        aux["lifetime_eu"],
        aux["lifetime_ricci"],
        aux["prominence_eu"],
        aux["prominence_ricci"],
        str(global_fig_dir / "cll_rs_scrna_h1_lifetime_eu_vs_ricci.png"),
        str(global_fig_dir / "cll_rs_scrna_h1_prominence_eu_vs_ricci.png"),
    )
    plot_f1_bars_with_ci(metrics_main, str(global_fig_dir / "cll_rs_scrna_f1_bars_ci.png"))
    if not ablation.empty:
        plot_lambda_robustness(ablation, str(global_fig_dir / "cll_rs_scrna_lambda_robustness.png"))

    print("Saved:", out / "metrics.csv")


if __name__ == "__main__":
    main()
