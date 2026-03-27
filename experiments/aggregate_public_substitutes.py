#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load_main_table(path: Path, dataset: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["dataset"] = dataset
    return df


def _load_ablation_table(path: Path, dataset: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["dataset"] = dataset
    return df


def _plot_main_benchmark(df: pd.DataFrame, out_path: Path) -> None:
    order = ["standard", "harmonic", "eu_sparse", "ccspr"]
    datasets = sorted(df["dataset"].unique())
    x = np.arange(len(order))
    w = 0.35 if len(datasets) <= 2 else 0.25

    fig, ax = plt.subplots(figsize=(9, 4))
    for i, ds in enumerate(datasets):
        sub = df[df["dataset"] == ds].set_index("model")
        y = [sub["mean"].get(m, np.nan) for m in order]
        e = [sub["ci95"].get(m, 0.0) if "ci95" in sub.columns else 0.0 for m in order]
        ax.bar(x + (i - (len(datasets) - 1) / 2.0) * w, y, width=w, yerr=e, capsize=3, label=ds)

    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Weighted F1")
    ax.set_title("CC-SPR Public Substitute Benchmark")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def _plot_ablation_panels(df: pd.DataFrame, out_path: Path) -> None:
    ablations = ["distance_mode", "iters", "k", "lambda", "top_k", "normalize"]
    datasets = sorted(df["dataset"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes = axes.flatten()

    for ax, abl in zip(axes, ablations):
        sub = df[df["ablation"] == abl].copy()
        if sub.empty:
            ax.set_axis_off()
            continue
        if abl in {"iters", "k", "top_k", "lambda"}:
            try:
                sub["ablation_value_num"] = sub["ablation_value"].astype(float)
            except Exception:
                sub["ablation_value_num"] = np.arange(len(sub))
            for ds in datasets:
                sds = sub[sub["dataset"] == ds].sort_values("ablation_value_num")
                ax.errorbar(
                    sds["ablation_value_num"],
                    sds["mean"],
                    yerr=sds["std"].fillna(0.0),
                    marker="o",
                    label=ds,
                )
            if abl == "lambda":
                ax.set_xscale("log")
            ax.set_xlabel(abl)
        else:
            vals = sorted(sub["ablation_value"].astype(str).unique())
            xpos = np.arange(len(vals))
            width = 0.35 if len(datasets) <= 2 else 0.25
            for i, ds in enumerate(datasets):
                sds = sub[sub["dataset"] == ds].copy()
                sds["ablation_value"] = sds["ablation_value"].astype(str)
                sds = sds.set_index("ablation_value")
                y = [sds["mean"].get(v, np.nan) for v in vals]
                e = [sds["std"].get(v, 0.0) for v in vals]
                ax.bar(xpos + (i - (len(datasets) - 1) / 2.0) * width, y, width=width, yerr=e, label=ds)
            ax.set_xticks(xpos)
            ax.set_xticklabels(vals)
            ax.set_xlabel(abl)

        ax.set_ylim(0, 1)
        ax.set_ylabel("Weighted F1")
        ax.set_title(f"Ablation: {abl}")

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=max(1, len(datasets)))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    tables_dir = Path("results/tables")
    figures_dir = Path("results/figures")
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    main_tables = [
        ("cll_venetoclax", Path("results/cll_venetoclax/main_results.csv")),
        ("cll_rs_scrna", Path("results/cll_rs_scrna/main_results.csv")),
    ]
    ablation_tables = [
        ("cll_venetoclax", Path("results/cll_venetoclax/ablation_summary.csv")),
        ("cll_rs_scrna", Path("results/cll_rs_scrna/ablation_summary.csv")),
    ]

    main_df = pd.concat([_load_main_table(path, ds) for ds, path in main_tables], ignore_index=True)
    abl_df = pd.concat([_load_ablation_table(path, ds) for ds, path in ablation_tables], ignore_index=True)

    main_df.to_csv(tables_dir / "main_results_public_substitutes.csv", index=False)
    abl_df.to_csv(tables_dir / "ablation_summary_public_substitutes.csv", index=False)

    _plot_main_benchmark(main_df, figures_dir / "main_benchmark_public_substitutes.png")
    _plot_ablation_panels(abl_df, figures_dir / "ablation_panels_public_substitutes.png")

    print("Saved:", tables_dir / "main_results_public_substitutes.csv")
    print("Saved:", tables_dir / "ablation_summary_public_substitutes.csv")
    print("Saved:", figures_dir / "main_benchmark_public_substitutes.png")
    print("Saved:", figures_dir / "ablation_panels_public_substitutes.png")


if __name__ == "__main__":
    main()
