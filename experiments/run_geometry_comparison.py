#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import gudhi
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedShuffleSplit

from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
from ccspr.datasets.cll_venetoclax import load_cll_venetoclax
from ccspr.datasets.tcga_brca import load_tcga_brca_multiomics
from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.geometry.distance import SUPPORTED_MODES, build_distance
from ccspr.preprocess.basic import select_top_variable_genes, standardize
from ccspr.solver.cycle import compute_cycle_representative
from ccspr.stability.tip import tip_bootstrap_topk
from ccspr.topology.persistence import build_cycle_problem, compute_persistence, pick_dominant_h1
from ccspr.utils.io import ensure_dir, read_yaml, save_json


@dataclass
class MethodSpec:
    name: str
    mode: str
    solver: str


@dataclass
class MethodParams:
    k: int
    iters: int
    alpha: float
    lambda_: float
    top_k: int
    normalize: str
    n_boot: int
    subsample_frac: float
    rescale_median: float
    bottleneck_boot: int
    bottleneck_top_bars: int
    bottleneck_datasets: tuple[str, ...]
    bottleneck_solvers: tuple[str, ...]


def _stratified_subsample(x: np.ndarray, y: np.ndarray, n_max: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if n_max <= 0 or x.shape[0] <= n_max:
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


def _topk_indices(scores: np.ndarray, k: int) -> np.ndarray:
    k = min(int(k), scores.shape[0])
    if k <= 0:
        return np.array([], dtype=int)
    return np.argsort(-scores)[:k]


def _fit_predict_logreg(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, seed: int) -> np.ndarray:
    clf = LogisticRegression(max_iter=5000, C=1.0, multi_class="auto", random_state=int(seed), n_jobs=1)
    clf.fit(x_train, y_train)
    return clf.predict(x_test)


def _safe_entropy(w: np.ndarray) -> float:
    a = np.abs(w)
    s = float(np.sum(a))
    if s <= 1e-12:
        return 0.0
    p = a / s
    p = p[p > 0]
    if p.size <= 1:
        return 0.0
    return float(-np.sum(p * np.log(p)) / np.log(p.size))


def _safe_gini(w: np.ndarray) -> float:
    a = np.abs(w)
    if np.sum(a) <= 1e-12:
        return 0.0
    x = np.sort(a)
    n = x.size
    cum = np.cumsum(x)
    return float((n + 1 - 2 * np.sum(cum) / cum[-1]) / n)


def _tip_concentration(tip: np.ndarray, m: int = 10) -> float:
    if tip.size == 0:
        return 0.0
    m = min(int(m), tip.size)
    den = float(np.sum(tip))
    if den <= 1e-12:
        return 0.0
    return float(np.sum(np.sort(tip)[::-1][:m]) / den)


def _bars_to_intervals(bars: list[dict[str, Any]], top_n: int) -> np.ndarray:
    ints = []
    for b in bars[: int(top_n)]:
        bt = float(b.get("bt", np.nan))
        dt = float(b.get("dt", np.nan))
        if np.isfinite(bt) and np.isfinite(dt) and dt > bt:
            ints.append([bt, dt])
    if not ints:
        return np.zeros((0, 2), dtype=float)
    return np.asarray(ints, dtype=float)


def _bootstrap_bottleneck(
    x: np.ndarray,
    method: MethodSpec,
    mp: MethodParams,
    seed: int,
    cache_dir: str,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = x.shape[0]
    n_take = max(6, int(mp.subsample_frac * n))
    diagrams: list[np.ndarray] = []

    for _ in range(int(mp.bottleneck_boot)):
        idx = rng.choice(n, n_take, replace=False)
        x_sub = x[idx]
        dist = build_distance(
            x_sub,
            mode=method.mode,
            k=mp.k,
            iters=mp.iters,
            alpha=mp.alpha,
            rescale_median=mp.rescale_median,
            cache_dir=cache_dir,
        )
        ph = compute_persistence(dist, max_dimension=2, cache_dir=cache_dir)
        diagrams.append(_bars_to_intervals(ph.get("bars", []), top_n=mp.bottleneck_top_bars))

    vals = []
    for i in range(len(diagrams)):
        for j in range(i + 1, len(diagrams)):
            try:
                vals.append(float(gudhi.bottleneck_distance(diagrams[i], diagrams[j])))
            except Exception:
                continue

    if not vals:
        return float("nan"), float("nan")
    return float(np.mean(vals)), float(np.std(vals))


def _single_cycle_stats(
    x: np.ndarray,
    method: MethodSpec,
    mp: MethodParams,
    cache_dir: str,
) -> dict[str, float]:
    dist, gdiag = build_distance(
        x,
        mode=method.mode,
        k=mp.k,
        iters=mp.iters,
        alpha=mp.alpha,
        rescale_median=mp.rescale_median,
        cache_dir=cache_dir,
        return_diagnostics=True,
    )
    ph = compute_persistence(dist, max_dimension=2, cache_dir=cache_dir)
    bars = ph.get("bars", [])

    top1 = float(bars[0]["life"]) if len(bars) >= 1 else 0.0
    top2 = float(bars[1]["life"]) if len(bars) >= 2 else 0.0
    prom = max(0.0, top1 - top2)

    support = 0
    entropy = 0.0
    gini = 0.0
    nz_mass = 0.0

    bar = pick_dominant_h1(ph)
    if bar is not None:
        edges, d_mat, z0 = build_cycle_problem(dist, bar)
        w = compute_cycle_representative(d_mat, z0, solver=method.solver, lambda_=mp.lambda_)
        a = np.abs(np.asarray(w, dtype=float))
        support = int(np.sum(a > 1e-8))
        entropy = _safe_entropy(a)
        gini = _safe_gini(a)
        nz_mass = float(np.sum(a))

    out = {
        "h1_lifetime_single": float(top1),
        "h1_prominence_single": float(prom),
        "h1_bar_count": int(len(bars)),
        "rep_support": int(support),
        "rep_entropy": float(entropy),
        "rep_gini": float(gini),
        "rep_mass": float(nz_mass),
        "curvature_mean": float(gdiag.get("curvature_mean", np.nan)),
        "curvature_std": float(gdiag.get("curvature_std", np.nan)),
        "neg_curvature_frac": float(gdiag.get("neg_curvature_frac", np.nan)),
    }
    return out


def _jaccard_sets(sets: list[set[int]]) -> float:
    if len(sets) < 2:
        return 1.0
    vals = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            vals.append(len(sets[i] & sets[j]) / max(1, len(u)))
    return float(np.mean(vals)) if vals else 1.0


def _build_methods(cfg: dict[str, Any]) -> list[MethodSpec]:
    modes = [str(m).lower() for m in cfg["geometry_modes"]]
    if bool(cfg.get("include_forman_ricci", False)) and "forman_ricci" in SUPPORTED_MODES:
        if "forman_ricci" not in modes:
            modes.append("forman_ricci")

    methods = []
    for mode in modes:
        for solver in cfg["solvers"]:
            methods.append(MethodSpec(name=f"{mode}_{solver}", mode=mode, solver=solver))
    return methods


def _load_dataset(name: str, cfg: dict[str, Any], seed: int) -> tuple[np.ndarray, np.ndarray]:
    t = cfg["type"]
    if t == "luad":
        ds = load_tcga_luad(root=cfg.get("root", "data/tcga_luad"), min_class_size=int(cfg.get("min_class_size", 10)))
        x, y = ds["X"], ds["y"]
        x, y = _stratified_subsample(x, y, n_max=int(cfg.get("max_samples", 0)), seed=seed)
        return x, y
    if t == "brca":
        ds = load_tcga_brca_multiomics(
            root=cfg.get("root", "data/tcga_brca"),
            pca_expr=int(cfg.get("pca_expr", 100)),
            pca_meth=int(cfg.get("pca_meth", 100)),
            max_samples=cfg.get("max_samples"),
            seed=seed,
        )
        x, y = ds["X"], ds["y"]
        x, y = _stratified_subsample(x, y, n_max=int(cfg.get("max_samples", 0)), seed=seed)
        return x, y
    if t == "cll_bulk":
        ds = load_cll_venetoclax(
            processed_matrix_path=cfg.get("processed_matrix_path", "data/cll_venetoclax/processed_matrix.tsv"),
            label_col=cfg.get("label_col", "label"),
        )
        x, y = ds["X"], ds["y"]
        x, y = _stratified_subsample(x, y, n_max=int(cfg.get("max_samples", 0)), seed=seed)
        return x, y
    if t == "cll_scrna":
        ds = load_cll_rs_scrna(
            data_path=cfg["data_path"],
            label_key=cfg.get("label_key", "time_group"),
            min_genes=int(cfg.get("min_genes", 200)),
            min_cells=int(cfg.get("min_cells", 3)),
            n_hvg=int(cfg.get("n_hvg", 2000)),
            n_pcs=int(cfg.get("n_pcs", 50)),
            max_cells=cfg.get("max_cells", 5000),
            use_precomputed_pca=bool(cfg.get("use_precomputed_pca", True)),
            seed=seed,
        )
        x, y = ds["X"], ds["y"]
        x, y = _stratified_subsample(x, y, n_max=int(cfg.get("max_samples", 0)), seed=seed)
        return x, y
    raise ValueError(f"Unsupported dataset type for {name}: {t}")


def _pipeline_figure(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.axis("off")
    boxes = [
        (0.02, 0.28, 0.13, 0.44, "Omics\nMatrix"),
        (0.20, 0.28, 0.18, 0.44, "Geometry Backends\nEu/Ricci/Diffusion\nPHATE/DTM"),
        (0.43, 0.28, 0.13, 0.44, "Persistent\nH1 Bars"),
        (0.61, 0.28, 0.14, 0.44, "Rep Solver\nL2/Elastic"),
        (0.80, 0.28, 0.14, 0.44, "TIP +\nClassifier"),
    ]
    for x, y, w, h, txt in boxes:
        rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", linewidth=1.4, facecolor="#edf4fb")
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=10)
    for i in range(len(boxes) - 1):
        x1 = boxes[i][0] + boxes[i][2]
        x2 = boxes[i + 1][0]
        ax.annotate("", xy=(x2 - 0.01, 0.50), xytext=(x1 + 0.01, 0.50), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.text(0.52, 0.10, "Geometry-sensitive persistent representatives: harmonic baseline + sparse extensions", ha="center", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _ci95(std: float, n: int) -> float:
    if n <= 1 or not np.isfinite(std):
        return 0.0
    return float(1.96 * std / math.sqrt(n))


def _plot_luad_benchmark(summary: pd.DataFrame, out_path: Path) -> None:
    sub = summary[summary["dataset"] == "tcga_luad"].copy()
    if sub.empty:
        return
    sub = sub.sort_values("mean", ascending=False)
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.bar(np.arange(len(sub)), sub["mean"], yerr=sub["ci95"], capsize=3)
    ax.set_xticks(np.arange(len(sub)))
    ax.set_xticklabels(sub["method"], rotation=30, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Weighted F1")
    ax.set_title("LUAD Main Benchmark (Geometry x Solver)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_geometry_ablation(summary: pd.DataFrame, out_path: Path) -> None:
    sub = summary[summary["method"].str.endswith("_elastic")].copy()
    if sub.empty:
        return
    sub["geometry"] = sub["method"].str.replace("_elastic", "", regex=False)
    datasets = sorted(sub["dataset"].unique())
    geoms = ["euclidean", "ricci", "diffusion", "phate", "dtm", "forman_ricci"]
    geoms = [g for g in geoms if g in set(sub["geometry"]) ]

    fig, axes = plt.subplots(1, len(datasets), figsize=(4.8 * len(datasets), 4), sharey=True)
    if len(datasets) == 1:
        axes = [axes]

    for ax, ds in zip(axes, datasets):
        d = sub[sub["dataset"] == ds].set_index("geometry")
        xs = np.arange(len(geoms))
        ys = [float(d["mean"].get(g, np.nan)) for g in geoms]
        es = [float(d["ci95"].get(g, 0.0)) for g in geoms]
        ax.bar(xs, ys, yerr=es, capsize=3)
        ax.set_xticks(xs)
        ax.set_xticklabels(geoms, rotation=30, ha="right")
        ax.set_title(ds)
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("Weighted F1")
    fig.suptitle("Geometry Ablation (Sparse/Elastic)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_h1(summary_top: pd.DataFrame, out_path: Path) -> None:
    sub = summary_top[summary_top["method"].str.endswith("_elastic")].copy()
    if sub.empty:
        return
    sub["geometry"] = sub["method"].str.replace("_elastic", "", regex=False)

    geoms = ["euclidean", "ricci", "diffusion", "phate", "dtm", "forman_ricci"]
    geoms = [g for g in geoms if g in set(sub["geometry"]) ]

    g = sub.groupby("geometry")[["h1_lifetime_single_mean", "h1_prominence_single_mean"]].mean().reindex(geoms)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(np.arange(len(g)), g["h1_lifetime_single_mean"])
    axes[0].set_xticks(np.arange(len(g)))
    axes[0].set_xticklabels(g.index, rotation=30, ha="right")
    axes[0].set_title("Dominant H1 Lifetime")

    axes[1].bar(np.arange(len(g)), g["h1_prominence_single_mean"])
    axes[1].set_xticks(np.arange(len(g)))
    axes[1].set_xticklabels(g.index, rotation=30, ha="right")
    axes[1].set_title("Dominant H1 Prominence")

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_rep(summary_rep: pd.DataFrame, out_path: Path) -> None:
    sub = summary_rep[summary_rep["method"].str.endswith("_elastic")].copy()
    if sub.empty:
        return
    sub["geometry"] = sub["method"].str.replace("_elastic", "", regex=False)

    geoms = ["euclidean", "ricci", "diffusion", "phate", "dtm", "forman_ricci"]
    geoms = [g for g in geoms if g in set(sub["geometry"]) ]
    g = sub.groupby("geometry")[["rep_support_mean", "rep_entropy_mean"]].mean().reindex(geoms)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(np.arange(len(g)), g["rep_support_mean"])
    axes[0].set_xticks(np.arange(len(g)))
    axes[0].set_xticklabels(g.index, rotation=30, ha="right")
    axes[0].set_title("Representative Support Size")

    axes[1].bar(np.arange(len(g)), g["rep_entropy_mean"])
    axes[1].set_xticks(np.arange(len(g)))
    axes[1].set_xticklabels(g.index, rotation=30, ha="right")
    axes[1].set_title("Representative Entropy")

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_tip(summary_rep: pd.DataFrame, out_path: Path) -> None:
    sub = summary_rep[summary_rep["method"].str.endswith("_elastic")].copy()
    if sub.empty:
        return
    sub["geometry"] = sub["method"].str.replace("_elastic", "", regex=False)

    datasets = sorted(sub["dataset"].unique())
    geoms = ["euclidean", "ricci", "diffusion", "phate", "dtm", "forman_ricci"]
    geoms = [g for g in geoms if g in set(sub["geometry"]) ]

    fig, axes = plt.subplots(1, len(datasets), figsize=(4.8 * len(datasets), 4), sharey=True)
    if len(datasets) == 1:
        axes = [axes]

    for ax, ds in zip(axes, datasets):
        d = sub[sub["dataset"] == ds].set_index("geometry")
        ys = [float(d["tip_split_jaccard_mean_mean"].get(g, np.nan)) for g in geoms]
        ax.bar(np.arange(len(geoms)), ys)
        ax.set_xticks(np.arange(len(geoms)))
        ax.set_xticklabels(geoms, rotation=30, ha="right")
        ax.set_title(ds)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Top-K Jaccard")
    fig.suptitle("TIP Stability Across Splits")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _plot_curvature_gain(df_rel: pd.DataFrame, out_path: Path) -> None:
    if df_rel.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for ds in sorted(df_rel["dataset"].unique()):
        d = df_rel[df_rel["dataset"] == ds]
        ax.scatter(d["neg_curvature_frac"], d["prominence_gain_vs_eu_sparse"], label=ds)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Negative curvature edge fraction")
    ax.set_ylabel("Prominence gain vs euclidean_elastic")
    ax.set_title("Curvature Statistic vs Topological Gain")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _concept_table() -> pd.DataFrame:
    rows = [
        ("harmonic_ph", "euclidean", "l2", "parent baseline", "dense representative tendency"),
        ("euclidean_sparse", "euclidean", "elastic", "sparsity control", "no geometry deformation"),
        ("ricci_sparse", "ricci", "elastic", "bridge/bottleneck sharpening", "curvature-flow sensitivity"),
        ("diffusion_sparse", "diffusion", "elastic", "noise-robust manifold smoothing", "time-scale dependence"),
        ("phate_sparse", "phate", "elastic", "potential-distance trajectories", "kernel/time tuning"),
        ("dtm_sparse", "dtm", "elastic", "outlier-robust topology", "weight-scale dependence"),
        ("forman_sparse", "forman_ricci", "elastic", "cheap curvature surrogate", "approximate curvature only"),
    ]
    return pd.DataFrame(rows, columns=["method_family", "geometry", "solver", "strength", "risk"])


def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", message="Solution may be inaccurate")

    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/geometry_paper.yaml")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)

    out_dir = ensure_dir(cfg["output"]["results_dir"])
    fig_dir = ensure_dir(out_dir / "figures")
    tab_dir = ensure_dir(out_dir / "tables")
    global_fig_dir = ensure_dir(cfg["output"].get("global_figures_dir", "results/figures"))
    global_tab_dir = ensure_dir(cfg["output"].get("global_tables_dir", "results/tables"))
    cache_dir = str(cfg.get("cache_dir", "data/cache"))

    save_json({"timestamp": datetime.utcnow().isoformat() + "Z", "config": cfg}, out_dir / "run.json")

    methods = _build_methods(cfg["methods"])
    mp = MethodParams(
        k=int(cfg["methods"]["k"]),
        iters=int(cfg["methods"]["iters"]),
        alpha=float(cfg["methods"]["alpha"]),
        lambda_=float(cfg["methods"]["lambda_"]),
        top_k=int(cfg["methods"]["top_k"]),
        normalize=str(cfg["methods"]["normalize"]),
        n_boot=int(cfg["methods"]["n_boot"]),
        subsample_frac=float(cfg["methods"]["subsample_frac"]),
        rescale_median=float(cfg["methods"]["rescale_median"]),
        bottleneck_boot=int(cfg["methods"].get("bottleneck_boot", 2)),
        bottleneck_top_bars=int(cfg["methods"].get("bottleneck_top_bars", 40)),
        bottleneck_datasets=tuple(str(x) for x in cfg["methods"].get("bottleneck_datasets", ["tcga_luad"])),
        bottleneck_solvers=tuple(str(x) for x in cfg["methods"].get("bottleneck_solvers", ["elastic"])),
    )

    rows: list[dict[str, Any]] = []

    for ds_id, dcfg in cfg["datasets"].items():
        if not bool(dcfg.get("enabled", True)):
            continue

        x, y = _load_dataset(ds_id, dcfg, seed=seed)
        splitter = StratifiedShuffleSplit(
            n_splits=int(dcfg.get("n_splits", 2)),
            test_size=float(dcfg.get("test_size", 0.25)),
            random_state=seed,
        )

        for split_id, (itr, ite) in enumerate(splitter.split(x, y)):
            x_train = x[itr]
            y_train = y[itr]
            x_test = x[ite]
            y_test = y[ite]

            x_train_var, idx_var = select_top_variable_genes(x_train, top_n=min(int(dcfg.get("var_genes", 500)), x.shape[1]))
            x_test_var = x_test[:, idx_var]

            xtr_std, xts_std = standardize(x_train_var, x_test_var)
            y_pred_std = _fit_predict_logreg(xtr_std, y_train, xts_std, seed + split_id)
            f1_std = float(f1_score(y_test, y_pred_std, average="weighted"))
            rows.append(
                {
                    "dataset": ds_id,
                    "split": int(split_id),
                    "method": "standard",
                    "geometry": "none",
                    "solver": "none",
                    "f1_weighted": f1_std,
                    "h1_lifetime_boot": np.nan,
                    "h1_prominence_boot": np.nan,
                    "h1_lifetime_single": np.nan,
                    "h1_prominence_single": np.nan,
                    "bottleneck_mean": np.nan,
                    "bottleneck_std": np.nan,
                    "rep_support": np.nan,
                    "rep_entropy": np.nan,
                    "rep_gini": np.nan,
                    "tip_concentration_top10": np.nan,
                    "tip_split_jaccard_mean": np.nan,
                    "curvature_mean": np.nan,
                    "curvature_std": np.nan,
                    "neg_curvature_frac": np.nan,
                    "topk_signature": "",
                }
            )

            for midx, m in enumerate(methods):
                tip_out = tip_bootstrap_topk(
                    x_train_var,
                    n_boot=mp.n_boot,
                    subsample_frac=mp.subsample_frac,
                    top_k=mp.top_k,
                    mode=m.mode,
                    k=mp.k,
                    iters=mp.iters,
                    alpha=mp.alpha,
                    solver=m.solver,
                    lambda_=mp.lambda_,
                    normalize=mp.normalize,
                    rescale_median=mp.rescale_median,
                    seed=seed + 1000 * (split_id + 1) + midx,
                    cache_dir=cache_dir,
                )

                tip = np.asarray(tip_out["tip"], dtype=float)
                idx = _topk_indices(tip, mp.top_k)

                if idx.size == 0:
                    f1 = float("nan")
                else:
                    xtr = x_train_var[:, idx]
                    xts = x_test_var[:, idx]
                    xtr, xts = standardize(xtr, xts)
                    pred = _fit_predict_logreg(xtr, y_train, xts, seed + 5000 + split_id + midx)
                    f1 = float(f1_score(y_test, pred, average="weighted"))

                cyc = _single_cycle_stats(x_train_var, m, mp, cache_dir=cache_dir)
                do_bottleneck = (ds_id in mp.bottleneck_datasets) and (m.solver in mp.bottleneck_solvers)
                if do_bottleneck:
                    b_mean, b_std = _bootstrap_bottleneck(
                        x_train_var,
                        m,
                        mp,
                        seed=seed + 2000 + split_id + midx,
                        cache_dir=cache_dir,
                    )
                else:
                    b_mean, b_std = float("nan"), float("nan")

                life_boot = tip_out["lifetime"]
                prom_boot = tip_out["prominence"]
                life_boot_mean = float(np.mean(life_boot)) if len(life_boot) else 0.0
                prom_boot_mean = float(np.mean(prom_boot)) if len(prom_boot) else 0.0

                rows.append(
                    {
                        "dataset": ds_id,
                        "split": int(split_id),
                        "method": m.name,
                        "geometry": m.mode,
                        "solver": m.solver,
                        "f1_weighted": f1,
                        "h1_lifetime_boot": life_boot_mean,
                        "h1_prominence_boot": prom_boot_mean,
                        "h1_lifetime_single": cyc["h1_lifetime_single"],
                        "h1_prominence_single": cyc["h1_prominence_single"],
                        "bottleneck_mean": b_mean,
                        "bottleneck_std": b_std,
                        "rep_support": cyc["rep_support"],
                        "rep_entropy": cyc["rep_entropy"],
                        "rep_gini": cyc["rep_gini"],
                        "tip_concentration_top10": _tip_concentration(tip, m=10),
                        "tip_split_jaccard_mean": np.nan,
                        "curvature_mean": cyc["curvature_mean"],
                        "curvature_std": cyc["curvature_std"],
                        "neg_curvature_frac": cyc["neg_curvature_frac"],
                        "topk_signature": ";".join(str(int(i)) for i in idx.tolist()),
                    }
                )

        partial = pd.DataFrame(rows)
        partial.to_csv(out_dir / "metrics_partial.csv", index=False)
        print(f"Completed dataset: {ds_id} (rows={len(partial)})", flush=True)

    metrics = pd.DataFrame(rows)

    # Split-level TIP stability by dataset/method using top-k overlap across splits.
    for (ds, method), block in metrics.groupby(["dataset", "method"]):
        if method == "standard":
            continue
        sets = []
        for sig in block["topk_signature"].tolist():
            if not sig:
                continue
            sets.append(set(int(x) for x in sig.split(";") if x != ""))
        stab = _jaccard_sets(sets) if sets else np.nan
        metrics.loc[(metrics["dataset"] == ds) & (metrics["method"] == method), "tip_split_jaccard_mean"] = stab

    metrics.to_csv(out_dir / "metrics.csv", index=False)
    metrics.to_csv(tab_dir / "geometry_methods_metrics.csv", index=False)
    metrics.to_csv(global_tab_dir / "geometry_methods_metrics.csv", index=False)

    summary = (
        metrics.groupby(["dataset", "method"]) ["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    summary["ci95"] = [ _ci95(float(s), int(n)) for s, n in zip(summary["std"].fillna(0.0), summary["n"]) ]

    top_cols = [
        "h1_lifetime_boot",
        "h1_prominence_boot",
        "h1_lifetime_single",
        "h1_prominence_single",
        "bottleneck_mean",
    ]
    top_summary = metrics.groupby(["dataset", "method"])[top_cols].mean().reset_index()
    top_summary = top_summary.rename(columns={c: f"{c}_mean" for c in top_cols})

    rep_cols = [
        "rep_support",
        "rep_entropy",
        "rep_gini",
        "tip_concentration_top10",
        "tip_split_jaccard_mean",
        "curvature_mean",
        "curvature_std",
        "neg_curvature_frac",
    ]
    rep_summary = metrics.groupby(["dataset", "method"])[rep_cols].mean().reset_index()
    rep_summary = rep_summary.rename(columns={c: f"{c}_mean" for c in rep_cols})

    summary.to_csv(tab_dir / "geometry_main_benchmark.csv", index=False)
    summary.to_csv(global_tab_dir / "geometry_main_benchmark.csv", index=False)
    top_summary.to_csv(tab_dir / "geometry_topology_summary.csv", index=False)
    top_summary.to_csv(global_tab_dir / "geometry_topology_summary.csv", index=False)
    rep_summary.to_csv(tab_dir / "geometry_representative_summary.csv", index=False)
    rep_summary.to_csv(global_tab_dir / "geometry_representative_summary.csv", index=False)

    # Paired split statistics.
    pairs = []
    ref = str(cfg["stats"].get("paired_reference", "euclidean_elastic"))
    targets = [str(x) for x in cfg["stats"].get("paired_targets", [])]

    nonstd = metrics[metrics["method"] != "standard"].copy()
    for ds, block in nonstd.groupby("dataset"):
        piv = block.pivot_table(index="split", columns="method", values="f1_weighted", aggfunc="mean")
        if ref not in piv.columns:
            continue
        for tar in targets:
            if tar not in piv.columns:
                continue
            d = (piv[tar] - piv[ref]).dropna()
            n = int(d.shape[0])
            if n == 0:
                continue
            mean_d = float(np.mean(d))
            std_d = float(np.std(d, ddof=1)) if n > 1 else 0.0
            ci = _ci95(std_d, n)

            t_stat, t_p = (np.nan, np.nan)
            if n > 1:
                t_res = ttest_rel(piv.loc[d.index, tar], piv.loc[d.index, ref], nan_policy="omit")
                t_stat, t_p = float(t_res.statistic), float(t_res.pvalue)

            w_stat, w_p = (np.nan, np.nan)
            if n > 1 and not np.allclose(d.values, 0.0):
                try:
                    wr = wilcoxon(d.values)
                    w_stat, w_p = float(wr.statistic), float(wr.pvalue)
                except Exception:
                    pass

            pairs.append(
                {
                    "dataset": ds,
                    "reference": ref,
                    "target": tar,
                    "n_splits": n,
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

    paired = pd.DataFrame(pairs)
    paired.to_csv(tab_dir / "geometry_paired_stats.csv", index=False)
    paired.to_csv(global_tab_dir / "geometry_paired_stats.csv", index=False)

    # Curvature-vs-prominence relation for sparse methods.
    rel_rows = []
    sparse = metrics[metrics["method"].str.endswith("_elastic")].copy()
    for (ds, split), block in sparse.groupby(["dataset", "split"]):
        row_eu = block[block["method"] == "euclidean_elastic"]
        row_ri = block[block["method"] == "ricci_elastic"]
        if row_eu.empty or row_ri.empty:
            continue
        rel_rows.append(
            {
                "dataset": ds,
                "split": int(split),
                "neg_curvature_frac": float(row_ri.iloc[0]["neg_curvature_frac"]),
                "prominence_gain_vs_eu_sparse": float(row_ri.iloc[0]["h1_prominence_single"] - row_eu.iloc[0]["h1_prominence_single"]),
                "f1_gain_vs_eu_sparse": float(row_ri.iloc[0]["f1_weighted"] - row_eu.iloc[0]["f1_weighted"]),
            }
        )
    rel = pd.DataFrame(rel_rows)
    rel.to_csv(tab_dir / "geometry_curvature_prominence_relation.csv", index=False)
    rel.to_csv(global_tab_dir / "geometry_curvature_prominence_relation.csv", index=False)

    concept = _concept_table()
    concept.to_csv(tab_dir / "geometry_method_concept_table.csv", index=False)
    concept.to_csv(global_tab_dir / "geometry_method_concept_table.csv", index=False)

    _pipeline_figure(fig_dir / "geometry_comparison_pipeline.png")
    _pipeline_figure(global_fig_dir / "geometry_comparison_pipeline.png")
    _plot_luad_benchmark(summary, fig_dir / "luad_geometry_benchmark.png")
    _plot_luad_benchmark(summary, global_fig_dir / "luad_geometry_benchmark.png")
    _plot_geometry_ablation(summary, fig_dir / "geometry_ablation_all_datasets.png")
    _plot_geometry_ablation(summary, global_fig_dir / "geometry_ablation_all_datasets.png")
    _plot_h1(top_summary, fig_dir / "geometry_h1_lifetime_prominence.png")
    _plot_h1(top_summary, global_fig_dir / "geometry_h1_lifetime_prominence.png")
    _plot_rep(rep_summary, fig_dir / "representative_sparsity_entropy.png")
    _plot_rep(rep_summary, global_fig_dir / "representative_sparsity_entropy.png")
    _plot_tip(rep_summary, fig_dir / "tip_stability_geometry.png")
    _plot_tip(rep_summary, global_fig_dir / "tip_stability_geometry.png")
    _plot_curvature_gain(rel, fig_dir / "curvature_prominence_gain.png")
    _plot_curvature_gain(rel, global_fig_dir / "curvature_prominence_gain.png")

    print("Saved:", out_dir / "metrics.csv")
    print("Saved summaries:", tab_dir)
    print("Saved figures:", fig_dir)


if __name__ == "__main__":
    main()
