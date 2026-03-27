from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ccspr.utils.io import ensure_dir


def plot_tip_eu_vs_ricci(
    tip_eu: np.ndarray,
    tip_ricci: np.ndarray,
    out_path: str,
    top_n: int = 100,
) -> None:
    ensure_dir("/".join(out_path.split("/")[:-1]))

    idx = np.argsort(-(tip_eu + tip_ricci))[: int(top_n)]
    x = np.arange(len(idx))
    w = 0.45

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(x - w / 2, tip_eu[idx], width=w, label="TIP euclidean")
    ax.bar(x + w / 2, tip_ricci[idx], width=w, label="TIP ricci")
    ax.set_title("TIP Comparison: Euclidean vs Ricci")
    ax.set_xlabel("Top features (by combined TIP)")
    ax.set_ylabel("TIP score")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close(fig)


def _box(a: np.ndarray, b: np.ndarray, labels, title: str, ylabel: str, out_path: str):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.boxplot([a, b], labels=labels, showfliers=False)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_lifetime_prominence(
    life_eu: np.ndarray,
    life_ri: np.ndarray,
    prom_eu: np.ndarray,
    prom_ri: np.ndarray,
    out_life: str,
    out_prom: str,
) -> None:
    _box(
        life_eu,
        life_ri,
        labels=["euclidean", "ricci"],
        title="Dominant H1 lifetime: euclidean vs ricci",
        ylabel="lifetime",
        out_path=out_life,
    )
    _box(
        prom_eu,
        prom_ri,
        labels=["euclidean", "ricci"],
        title="Dominant H1 prominence: euclidean vs ricci",
        ylabel="prominence",
        out_path=out_prom,
    )


def plot_f1_bars_with_ci(metrics: pd.DataFrame, out_path: str) -> None:
    if metrics.empty:
        return
    model_order = ["standard", "harmonic", "eu_sparse", "ccspr"]
    sub = metrics[metrics["model"].isin(model_order)].copy()
    g = sub.groupby("model")["f1_weighted"]
    mean = g.mean()
    ci = 1.96 * g.std(ddof=1).fillna(0.0) / np.sqrt(g.count().clip(lower=1))

    order = [m for m in model_order if m in set(sub["model"].unique())]
    m = [mean.get(o, np.nan) for o in order]
    c = [ci.get(o, 0.0) for o in order]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(order, m, yerr=c, capsize=4)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Weighted F1")
    ax.set_title("Weighted F1 with 95% CI")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close(fig)


def plot_lambda_robustness(ablation: pd.DataFrame, out_path: str) -> None:
    sub = ablation[ablation["ablation"] == "lambda"].copy()
    if sub.empty:
        return
    sub["lambda"] = sub["ablation_value"].astype(float)
    g = sub.groupby("lambda")["f1_weighted"].agg(["mean", "std", "count"]).reset_index()
    g["se"] = g["std"].fillna(0.0) / np.sqrt(g["count"].clip(lower=1))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(g["lambda"], g["mean"], yerr=1.96 * g["se"], marker="o")
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlabel("lambda")
    ax.set_ylabel("Weighted F1")
    ax.set_title("Lambda robustness")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close(fig)
