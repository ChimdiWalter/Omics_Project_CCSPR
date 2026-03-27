#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
from ccspr.datasets.cll_venetoclax import load_cll_venetoclax
from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.eval.classification import FeatureParams
from ccspr.geometry.distance import build_distance
from ccspr.preprocess.basic import select_top_variable_genes
from ccspr.solver.cycle import compute_cycle_representative
from ccspr.solver.projection import project_to_features
from ccspr.stability.tip import tip_bootstrap_topk
from ccspr.topology.persistence import build_cycle_problem, compute_persistence, pick_dominant_h1
from ccspr.utils.io import ensure_dir, read_yaml


def _stratified_small(x: np.ndarray, y: np.ndarray, n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if len(y) <= n:
        return x, y
    rng = np.random.default_rng(seed)
    keep = []
    cls, cnt = np.unique(y, return_counts=True)
    for c, k in zip(cls, cnt):
        idx = np.where(y == c)[0]
        take = max(2, int(round((k / len(y)) * n)))
        take = min(take, len(idx))
        keep.extend(rng.choice(idx, size=take, replace=False).tolist())
    keep = np.array(sorted(set(keep)), dtype=int)
    if len(keep) > n:
        keep = rng.choice(keep, size=n, replace=False)
    return x[keep], y[keep]


def _single_scores(x: np.ndarray, p: FeatureParams, cache_dir: str) -> np.ndarray:
    dist = build_distance(
        x,
        mode=p.mode,
        k=p.k,
        iters=p.iters,
        alpha=p.alpha,
        rescale_median=p.rescale_median,
        cache_dir=cache_dir,
    )
    ph = compute_persistence(dist, cache_dir=cache_dir)
    bar = pick_dominant_h1(ph)
    if bar is None:
        return np.zeros(x.shape[1])
    edges, d_mat, z0 = build_cycle_problem(dist, bar)
    w = compute_cycle_representative(d_mat, z0, solver=p.solver, lambda_=p.lambda_)
    return project_to_features(x, edges, w, normalize=p.normalize)


def _jaccard(a: set[int], b: set[int]) -> float:
    return len(a & b) / max(1, len(a | b))


def _pipeline_schematic(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 3))
    ax.axis("off")
    labels = ["X", "Distance\n(Eu/Ricci)", "PH + H1", "Cycle\n(L2/Elastic)", "Projection", "TIP + Classifier"]
    xs = np.linspace(0.05, 0.92, len(labels))
    for i, (x, label) in enumerate(zip(xs, labels)):
        r = patches.FancyBboxPatch((x - 0.06, 0.35), 0.12, 0.3, boxstyle="round,pad=0.02", facecolor="#eef3fb")
        ax.add_patch(r)
        ax.text(x, 0.5, label, ha="center", va="center", fontsize=10)
        if i < len(labels) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.07, 0.5), xytext=(x + 0.07, 0.5), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.text(0.5, 0.18, "CC-SPR extends harmonic PH with Ricci geometry and sparse cycle representatives", ha="center")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _synthetic(fig_dir: Path, tab_dir: Path, cache_dir: str) -> None:
    rng = np.random.default_rng(7)
    n, p = 70, 100
    t = np.linspace(0, 2 * np.pi, n, endpoint=False) + rng.normal(0, 0.05, n)
    x = rng.normal(0, 0.3, size=(n, p))
    informative = np.arange(10)
    x[:, 0] = np.cos(t) + rng.normal(0, 0.05, n)
    x[:, 1] = np.sin(t) + rng.normal(0, 0.05, n)
    for j in informative[2:]:
        x[:, j] = x[:, j % 2] + rng.normal(0, 0.08, n)

    methods = {
        "harmonic": FeatureParams("euclidean", "l2", 10, 0, 0.5, 1e-6, 20, "std", 1, 0.8, 1.0),
        "eu_sparse": FeatureParams("euclidean", "elastic", 10, 0, 0.5, 1e-6, 20, "std", 1, 0.8, 1.0),
        "ccspr": FeatureParams("ricci", "elastic", 10, 5, 0.5, 1e-6, 20, "std", 1, 0.8, 1.0),
    }
    rows = []
    curves = {}
    for name, pm in methods.items():
        s = _single_scores(x, pm, cache_dir)
        curves[name] = np.sort(s)[::-1][:40]
        top = set(np.argsort(-s)[:20].tolist())
        rec = len(top & set(informative.tolist())) / len(informative)
        rows.append({"method": name, "recall_at20": rec, "support_gt_1e_2": int(np.sum(s > 1e-2))})
    df = pd.DataFrame(rows)
    df.to_csv(tab_dir / "synthetic_metrics.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(np.cos(t), np.sin(t), c=t, cmap="viridis", s=20)
    axes[0].set_title("Synthetic Ring")
    axes[0].axis("equal")
    for name in ["harmonic", "eu_sparse", "ccspr"]:
        axes[1].plot(curves[name], label=name)
    axes[1].set_title("Top-40 Feature Score Decay")
    axes[1].set_xlabel("rank")
    axes[1].set_ylabel("score")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "synthetic_proof_of_concept.png", dpi=220)
    plt.close(fig)


def _load_specs(cache_dir: str):
    luad_cfg = read_yaml("configs/tcga_luad.yaml")
    luad = load_tcga_luad(root=luad_cfg["dataset"].get("root", "data/tcga_luad"))
    x_luad, y_luad = _stratified_small(luad["X"], luad["y"], 36, 42)
    luad_h = FeatureParams(**luad_cfg["protocol"]["harmonic_params"])
    luad_c = FeatureParams(**luad_cfg["protocol"]["ccspr_params"])

    ven_cfg = read_yaml("configs/cll_venetoclax.yaml")
    ven = load_cll_venetoclax(ven_cfg["dataset"]["processed_matrix_path"], ven_cfg["dataset"]["label_col"])
    x_ven, y_ven = _stratified_small(ven["X"], ven["y"], 36, 43)
    ven_h = FeatureParams(**ven_cfg["protocol"]["harmonic_params"])
    ven_c = FeatureParams(**ven_cfg["protocol"]["ccspr_params"])

    scr_cfg = read_yaml("configs/cll_rs_scrna.yaml")
    scr = load_cll_rs_scrna(
        data_path=scr_cfg["dataset"]["data_path"],
        label_key=scr_cfg["dataset"]["label_key"],
        min_genes=int(scr_cfg["dataset"]["min_genes"]),
        min_cells=int(scr_cfg["dataset"]["min_cells"]),
        n_hvg=int(scr_cfg["dataset"]["n_hvg"]),
        n_pcs=int(scr_cfg["dataset"]["n_pcs"]),
        max_cells=scr_cfg["dataset"]["max_cells"],
        use_precomputed_pca=bool(scr_cfg["dataset"].get("use_precomputed_pca", True)),
        seed=int(scr_cfg["seed"]),
    )
    x_scr, y_scr = _stratified_small(scr["X"], scr["y"], 36, 44)
    scr_h = FeatureParams(**scr_cfg["protocol"]["harmonic_params"])
    scr_c = FeatureParams(**scr_cfg["protocol"]["ccspr_params"])

    return [
        ("tcga_luad", x_luad, y_luad, 500, luad_h, luad_c),
        ("gse161711_bulk", x_ven, y_ven, 500, ven_h, ven_c),
        ("gse165087_scrna", x_scr, y_scr, 50, scr_h, scr_c),
    ]


def _stability_sparsity(fig_dir: Path, tab_dir: Path, cache_dir: str) -> None:
    stab_rows, sparse_rows = [], []
    luad_h, luad_c = None, None
    for name, x, y, vg, h, c in _load_specs(cache_dir):
        xv, _ = select_top_variable_genes(x, top_n=min(vg, x.shape[1]))
        for method, pm in [("harmonic", h), ("ccspr", c)]:
            outs = []
            sets = []
            for sd in [11, 23]:
                out = tip_bootstrap_topk(
                    xv,
                    n_boot=1,
                    subsample_frac=0.75,
                    top_k=20,
                    mode=pm.mode,
                    k=pm.k,
                    iters=pm.iters,
                    alpha=pm.alpha,
                    solver=pm.solver,
                    lambda_=pm.lambda_,
                    normalize=pm.normalize,
                    rescale_median=pm.rescale_median,
                    seed=sd,
                    cache_dir=cache_dir,
                )
                outs.append(out)
                sets.append(set(np.argsort(-out["tip"])[:20].tolist()))
            score = _single_scores(xv, pm, cache_dir)
            stab_rows.append(
                {
                    "dataset": name,
                    "method": method,
                    "tip_topk_jaccard_mean": _jaccard(sets[0], sets[1]),
                    "lifetime_mean": float(np.mean([np.mean(o["lifetime"]) if o["lifetime"].size else 0 for o in outs])),
                    "prominence_mean": float(np.mean([np.mean(o["prominence"]) if o["prominence"].size else 0 for o in outs])),
                }
            )
            sparse_rows.append(
                {
                    "dataset": name,
                    "method": method,
                    "support_gt_1e_2": int(np.sum(score > 1e-2)),
                    "support_gt_5e_2": int(np.sum(score > 5e-2)),
                    "top20_mass": float(np.sum(np.sort(score)[::-1][:20]) / max(np.sum(score), 1e-12)),
                }
            )
            if name == "tcga_luad":
                if method == "harmonic":
                    luad_h = np.sort(score)[::-1][:50]
                else:
                    luad_c = np.sort(score)[::-1][:50]

    stab = pd.DataFrame(stab_rows)
    sparse = pd.DataFrame(sparse_rows)
    stab.to_csv(tab_dir / "stability_metrics.csv", index=False)
    sparse.to_csv(tab_dir / "sparsity_metrics.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    stab.pivot(index="dataset", columns="method", values="tip_topk_jaccard_mean").plot(kind="bar", ax=axes[0], ylim=(0, 1))
    axes[0].set_title("TIP Stability (Jaccard)")
    axes[0].set_ylabel("Jaccard")
    axes[0].set_xlabel("")
    sparse.pivot(index="dataset", columns="method", values="support_gt_1e_2").plot(kind="bar", ax=axes[1])
    axes[1].set_title("Feature Support >1e-2")
    axes[1].set_ylabel("count")
    axes[1].set_xlabel("")
    fig.tight_layout()
    fig.savefig(fig_dir / "harmonic_vs_ccspr_sparsity_stability.png", dpi=220)
    plt.close(fig)

    if luad_h is not None and luad_c is not None:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(luad_h, label="harmonic")
        ax.plot(luad_c, label="ccspr")
        ax.set_title("LUAD Feature Support: Harmonic vs CC-SPR")
        ax.set_xlabel("feature rank")
        ax.set_ylabel("score")
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "luad_feature_support_harmonic_vs_ccspr.png", dpi=220)
        plt.close(fig)


def _robustness(fig_dir: Path) -> None:
    src = {
        "tcga_luad": "results/ablation_summary.csv",
        "gse161711_bulk": "results/cll_venetoclax/ablation_summary.csv",
        "gse165087_scrna": "results/cll_rs_scrna/ablation_summary.csv",
    }
    rows = []
    for ds, p in src.items():
        d = pd.read_csv(p)
        d["dataset"] = ds
        rows.append(d)
    df = pd.concat(rows, ignore_index=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, abl in zip(axes, ["lambda", "k", "iters"]):
        sub = df[df["ablation"] == abl].copy()
        sub["v"] = pd.to_numeric(sub["ablation_value"], errors="coerce")
        for ds in sorted(sub["dataset"].unique()):
            s = sub[sub["dataset"] == ds].sort_values("v")
            ax.plot(s["v"], s["mean"], marker="o", label=ds)
        if abl == "lambda":
            ax.set_xscale("log")
        ax.set_title(abl)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Weighted F1")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "robustness_lambda_k_iters_all.png", dpi=220)
    plt.close(fig)


def _concept_table(tab_dir: Path) -> None:
    rows = [
        ("metric assumption", "flat Euclidean", "Euclidean or Ricci-corrected"),
        ("representative criterion", "L2 harmonic representative", "elastic sparse representative"),
        ("sparsity", "implicit", "explicit L1 component"),
        ("stability", "not TIP-centric", "TIP bootstrap top-K"),
        ("feature attribution", "can smear", "sharper support"),
        ("geometry sensitivity", "limited", "curvature-aware"),
        ("implementation stack", "maTilDA/R-centric", "Python stack here; optional R parity hook"),
        ("datasets in this repo", "LUAD framing", "LUAD, BRCA, GSE161711, GSE165087"),
    ]
    pd.DataFrame(rows, columns=["dimension", "harmonic_ph", "ccspr"]).to_csv(
        tab_dir / "harmonic_vs_ccspr_comparison.csv", index=False
    )


def main() -> None:
    fig_dir = ensure_dir("results/figures")
    tab_dir = ensure_dir("results/tables")
    cache_dir = "data/cache"
    _pipeline_schematic(fig_dir / "pipeline_schematic_ccspr.png")
    _synthetic(fig_dir, tab_dir, cache_dir)
    _stability_sparsity(fig_dir, tab_dir, cache_dir)
    _robustness(fig_dir)
    _concept_table(tab_dir)
    print("saved fast diagnostics")


if __name__ == "__main__":
    main()
