#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import ttest_rel, wilcoxon

from ccspr.utils.io import ensure_dir, read_yaml


def _ci95(std: float, n: int) -> float:
    if n <= 1 or not np.isfinite(std):
        return 0.0
    return float(1.96 * std / math.sqrt(n))


def _build_seed_config(base_cfg: dict, seed: int, out_root: Path, n_splits: int, max_samples: int, var_genes: int) -> dict:
    cfg = copy.deepcopy(base_cfg)
    cfg["seed"] = int(seed)

    ds = cfg["datasets"]
    for k in list(ds.keys()):
        ds[k]["enabled"] = (k == "tcga_luad")

    ds["tcga_luad"]["n_splits"] = int(n_splits)
    ds["tcga_luad"]["max_samples"] = int(max_samples)
    ds["tcga_luad"]["var_genes"] = int(var_genes)

    out_dir = out_root / f"seed_{seed}"
    cfg["output"]["results_dir"] = str(out_dir)
    cfg["output"]["global_tables_dir"] = str(out_dir / "tables_global")
    cfg["output"]["global_figures_dir"] = str(out_dir / "figures_global")
    cfg["methods"]["bottleneck_datasets"] = []
    cfg["methods"]["bottleneck_boot"] = 0
    return cfg


def _run_seed(config_path: Path, log_path: Path) -> None:
    env = dict(os.environ)
    py = env.get("PYTHON", sys.executable)
    cmd = [py, "experiments/run_geometry_comparison.py", "--config", str(config_path)]
    ensure_dir(log_path.parent)
    with open(log_path, "w", encoding="utf-8") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"Seed run failed: {config_path} (see {log_path})")


def _paired_stats(df: pd.DataFrame, ref: str, targets: list[str]) -> pd.DataFrame:
    rows = []
    pivot = df.pivot_table(index="replicate", columns="method", values="f1_weighted", aggfunc="mean")
    if ref not in pivot.columns:
        return pd.DataFrame(columns=["reference", "target", "n", "mean_diff", "std_diff", "ci95", "win_rate", "t_stat", "t_pvalue", "wilcoxon_stat", "wilcoxon_pvalue"])

    for tar in targets:
        if tar not in pivot.columns:
            continue
        d = (pivot[tar] - pivot[ref]).dropna()
        n = int(d.shape[0])
        if n == 0:
            continue
        mean_d = float(np.mean(d))
        std_d = float(np.std(d, ddof=1)) if n > 1 else 0.0
        ci = _ci95(std_d, n)

        t_stat, t_p = (np.nan, np.nan)
        if n > 1:
            tr = ttest_rel(pivot.loc[d.index, tar], pivot.loc[d.index, ref], nan_policy="omit")
            t_stat, t_p = float(tr.statistic), float(tr.pvalue)

        w_stat, w_p = (np.nan, np.nan)
        if n > 1 and not np.allclose(d.values, 0.0):
            try:
                wr = wilcoxon(d.values)
                w_stat, w_p = float(wr.statistic), float(wr.pvalue)
            except Exception:
                pass

        rows.append(
            {
                "reference": ref,
                "target": tar,
                "n": n,
                "mean_diff": mean_d,
                "std_diff": std_d,
                "ci95": ci,
                "win_rate": float(np.mean(d.values > 0)),
                "t_stat": t_stat,
                "t_pvalue": t_p,
                "wilcoxon_stat": w_stat,
                "wilcoxon_pvalue": w_p,
            }
        )

    return pd.DataFrame(rows)


def _plot_benchmark(summary: pd.DataFrame, out_path: Path) -> None:
    if summary.empty:
        return
    s = summary.sort_values("mean", ascending=False)
    fig, ax = plt.subplots(figsize=(11.0, 4.2))
    ax.bar(np.arange(len(s)), s["mean"], yerr=s["ci95"], capsize=3)
    ax.set_xticks(np.arange(len(s)))
    ax.set_xticklabels(s["method"], rotation=30, ha="right")
    ax.set_ylabel("Weighted F1")
    ax.set_ylim(0, 1)
    ax.set_title("TCGA-LUAD Multi-seed Geometry Benchmark")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-config", default="configs/geometry_paper_ultrafast.yaml")
    ap.add_argument("--seeds", nargs="+", type=int, default=[11, 23, 42, 77, 101])
    ap.add_argument("--n-splits", type=int, default=3)
    ap.add_argument("--max-samples", type=int, default=36)
    ap.add_argument("--var-genes", type=int, default=200)
    ap.add_argument("--out-root", default="results/luad_multiseed")
    args = ap.parse_args()

    base_cfg = read_yaml(args.base_config)
    out_root = ensure_dir(args.out_root)
    logs_dir = ensure_dir("results/logs")
    ensure_dir("results/tables")
    ensure_dir("results/figures")

    all_metrics = []

    for seed in args.seeds:
        cfg = _build_seed_config(
            base_cfg=base_cfg,
            seed=int(seed),
            out_root=out_root,
            n_splits=int(args.n_splits),
            max_samples=int(args.max_samples),
            var_genes=int(args.var_genes),
        )
        cfg_path = out_root / f"config_seed_{seed}.yaml"
        with open(cfg_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(cfg, f, sort_keys=False)

        log_path = logs_dir / f"luad_multiseed_seed_{seed}.log"
        _run_seed(cfg_path, log_path)

        m_path = Path(cfg["output"]["results_dir"]) / "metrics.csv"
        if not m_path.exists():
            raise FileNotFoundError(f"Missing metrics: {m_path}")
        df = pd.read_csv(m_path)
        df = df[df["dataset"] == "tcga_luad"].copy()
        df["seed"] = int(seed)
        df["replicate"] = [f"{seed}_{int(s)}" for s in df["split"].tolist()]
        all_metrics.append(df)

    metrics = pd.concat(all_metrics, ignore_index=True)
    metrics.to_csv("results/tables/luad_multiseed_geometry_metrics.csv", index=False)

    summary = (
        metrics.groupby("method")["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    summary["ci95"] = [_ci95(float(s), int(n)) for s, n in zip(summary["std"].fillna(0.0), summary["n"])]
    summary.to_csv("results/tables/luad_multiseed_main_benchmark.csv", index=False)

    targets = ["ricci_elastic", "diffusion_elastic", "phate_elastic", "dtm_elastic", "euclidean_l2"]
    paired = _paired_stats(metrics[metrics["method"] != "standard"], ref="euclidean_elastic", targets=targets)
    paired.to_csv("results/tables/luad_multiseed_paired_stats.csv", index=False)

    _plot_benchmark(summary, Path("results/figures/luad_multiseed_benchmark.png"))

    meta = {
        "base_config": args.base_config,
        "seeds": [int(s) for s in args.seeds],
        "n_splits": int(args.n_splits),
        "max_samples": int(args.max_samples),
        "var_genes": int(args.var_genes),
    }
    with open("results/luad_multiseed/run_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print("Saved: results/tables/luad_multiseed_geometry_metrics.csv")


if __name__ == "__main__":
    main()
