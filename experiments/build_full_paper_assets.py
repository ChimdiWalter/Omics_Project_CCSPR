#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
from ccspr.datasets.cll_venetoclax import load_cll_venetoclax
from ccspr.datasets.tcga_brca import load_tcga_brca_multiomics
from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.eval.classification import FeatureParams, evaluate_luad_protocol
from ccspr.geometry.distance import build_distance
from ccspr.preprocess.basic import select_top_variable_genes
from ccspr.solver.cycle import compute_cycle_representative
from ccspr.solver.projection import project_to_features
from ccspr.stability.tip import tip_bootstrap_topk
from ccspr.topology.persistence import build_cycle_problem, compute_persistence, pick_dominant_h1
from ccspr.utils.io import ensure_dir, read_yaml


@dataclass
class DatasetSpec:
    name: str
    x: np.ndarray
    y: np.ndarray
    var_genes: int
    harmonic: FeatureParams
    ccspr: FeatureParams


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


def _summary_from_metrics(metrics: pd.DataFrame, dataset: str) -> pd.DataFrame:
    g = (
        metrics.groupby("model")["f1_weighted"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    g["ci95"] = 1.96 * g["std"].fillna(0.0) / np.sqrt(g["n"].clip(lower=1))
    g["dataset"] = dataset
    return g


def _run_extended_main_luad(cache_dir: str) -> pd.DataFrame:
    cfg = read_yaml("configs/tcga_luad.yaml")
    ds = load_tcga_luad(
        root=cfg["dataset"].get("root", "data/tcga_luad"),
        min_class_size=int(cfg["dataset"].get("min_class_size", 10)),
    )
    x, y = ds["X"], ds["y"]
    n_max = int(cfg["dataset"].get("max_samples", 0))
    if n_max > 0:
        x, y = _stratified_subsample(x, y, n_max=n_max, seed=int(cfg["seed"]))

    harmonic = FeatureParams(**cfg["protocol"]["harmonic_params"])
    ccspr = FeatureParams(**cfg["protocol"]["ccspr_params"])
    metrics, _ = evaluate_luad_protocol(
        x,
        y,
        var_genes=int(cfg["protocol"]["var_genes"]),
        n_splits=int(cfg["protocol"]["n_splits"]),
        test_size=float(cfg["protocol"]["test_size"]),
        seed=int(cfg["seed"]),
        cache_dir=cache_dir,
        harmonic_params=harmonic,
        ccspr_params=ccspr,
    )
    return metrics


def _run_extended_main_brca(cache_dir: str) -> pd.DataFrame:
    cfg = read_yaml("configs/tcga_brca.yaml")
    ds = load_tcga_brca_multiomics(
        root=cfg["dataset"].get("root", "data/tcga_brca"),
        pca_expr=int(cfg["dataset"].get("pca_expr", 100)),
        pca_meth=int(cfg["dataset"].get("pca_meth", 100)),
        max_samples=cfg["dataset"].get("max_samples"),
        seed=int(cfg["seed"]),
    )
    harmonic = FeatureParams(**cfg["protocol"]["harmonic_params"])
    ccspr = FeatureParams(**cfg["protocol"]["ccspr_params"])
    metrics, _ = evaluate_luad_protocol(
        ds["X"],
        ds["y"],
        var_genes=min(int(cfg["protocol"]["var_genes"]), ds["X"].shape[1]),
        n_splits=int(cfg["protocol"]["n_splits"]),
        test_size=float(cfg["protocol"]["test_size"]),
        seed=int(cfg["seed"]),
        cache_dir=cache_dir,
        harmonic_params=harmonic,
        ccspr_params=ccspr,
    )
    return metrics


def _plot_main_benchmark_all(df: pd.DataFrame, out_path: Path) -> None:
    order = ["standard", "harmonic", "eu_sparse", "ccspr"]
    datasets = list(df["dataset"].unique())
    n_d = len(datasets)
    fig, axes = plt.subplots(1, n_d, figsize=(4.4 * n_d, 4), sharey=True)
    if n_d == 1:
        axes = [axes]
    for ax, ds in zip(axes, datasets):
        sub = df[df["dataset"] == ds].set_index("model")
        x = np.arange(len(order))
        y = [sub["mean"].get(m, np.nan) for m in order]
        e = [sub["ci95"].get(m, 0.0) for m in order]
        ax.bar(x, y, yerr=e, capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=25, ha="right")
        ax.set_ylim(0, 1.0)
        ax.set_title(ds)
        ax.set_ylabel("Weighted F1")
    fig.suptitle("Main Benchmark Across Datasets")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _pipeline_schematic(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.axis("off")
    boxes = [
        (0.02, 0.3, 0.14, 0.4, "Omics Matrix\nX"),
        (0.20, 0.3, 0.14, 0.4, "Distance\nEuclidean/Ricci"),
        (0.38, 0.3, 0.14, 0.4, "PH + H1 Bar\nSelection"),
        (0.56, 0.3, 0.14, 0.4, "Cycle Solve\nL2/Elastic"),
        (0.74, 0.3, 0.14, 0.4, "Feature\nProjection"),
        (0.90, 0.3, 0.08, 0.4, "TIP\n+ CLS"),
    ]
    for x, y, w, h, txt in boxes:
        rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", linewidth=1.4, facecolor="#f3f6fb")
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=10)
    for i in range(len(boxes) - 1):
        x1 = boxes[i][0] + boxes[i][2]
        x2 = boxes[i + 1][0]
        ax.annotate("", xy=(x2 - 0.01, 0.5), xytext=(x1 + 0.01, 0.5), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.text(0.56, 0.1, "CC-SPR differs from harmonic PH baseline by Ricci geometry + elastic sparsity", ha="center", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _single_score_vector(x: np.ndarray, params: FeatureParams, cache_dir: str):
    dist = build_distance(
        x,
        mode=params.mode,
        k=params.k,
        iters=params.iters,
        alpha=params.alpha,
        rescale_median=params.rescale_median,
        cache_dir=cache_dir,
    )
    ph = compute_persistence(dist, max_dimension=2, cache_dir=cache_dir)
    bar = pick_dominant_h1(ph)
    if bar is None:
        return np.zeros(x.shape[1]), np.array([]), np.array([])
    edges, d_mat, z0 = build_cycle_problem(dist, bar)
    w = compute_cycle_representative(d_mat, z0, solver=params.solver, lambda_=params.lambda_)
    scores = project_to_features(x, edges, w, normalize=params.normalize)
    return scores, np.array(w), np.array(z0)


def _synthetic_data(seed: int = 7, n: int = 140, p: int = 180):
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    t = t + rng.normal(0, 0.07, size=n)
    z = np.stack([np.cos(t), np.sin(t)], axis=1)
    x = rng.normal(0, 0.35, size=(n, p))
    informative = np.arange(12)
    x[:, 0] = np.cos(t) + rng.normal(0, 0.08, size=n)
    x[:, 1] = np.sin(t) + rng.normal(0, 0.08, size=n)
    x[:, 2] = np.cos(2 * t) + rng.normal(0, 0.08, size=n)
    x[:, 3] = np.sin(2 * t) + rng.normal(0, 0.08, size=n)
    for j in informative[4:]:
        x[:, j] = x[:, j % 4] + rng.normal(0, 0.05, size=n)
    y = (np.sin(t) > 0).astype(int).astype(str)
    return x, y, z, informative


def _recall_at_k(scores: np.ndarray, true_idx: np.ndarray, k: int = 20):
    idx = np.argsort(-scores)[:k]
    hit = len(set(idx).intersection(set(true_idx.tolist())))
    return hit / max(1, min(k, len(true_idx)))


def _pairwise_jaccard(sets: list[set[int]]) -> float:
    if len(sets) < 2:
        return 1.0
    vals = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            inter = sets[i] & sets[j]
            vals.append(len(inter) / max(1, len(u)))
    return float(np.mean(vals))


def _synthetic_bundle(fig_dir: Path, tab_dir: Path, cache_dir: str):
    x, y, z, informative = _synthetic_data()
    harmonic = FeatureParams(
        mode="euclidean",
        solver="l2",
        k=10,
        iters=0,
        alpha=0.5,
        lambda_=1e-6,
        top_k=20,
        normalize="std",
        n_boot=6,
        subsample_frac=0.8,
        rescale_median=1.0,
    )
    eu_sparse = FeatureParams(
        mode="euclidean",
        solver="elastic",
        k=10,
        iters=0,
        alpha=0.5,
        lambda_=1e-6,
        top_k=20,
        normalize="std",
        n_boot=6,
        subsample_frac=0.8,
        rescale_median=1.0,
    )
    ccspr = FeatureParams(
        mode="ricci",
        solver="elastic",
        k=10,
        iters=5,
        alpha=0.5,
        lambda_=1e-6,
        top_k=20,
        normalize="std",
        n_boot=6,
        subsample_frac=0.8,
        rescale_median=1.0,
    )
    methods = [("harmonic", harmonic), ("eu_sparse", eu_sparse), ("ccspr", ccspr)]

    rows = []
    score_map = {}
    for name, p in methods:
        s, w, _ = _single_score_vector(x, p, cache_dir=cache_dir)
        out = tip_bootstrap_topk(
            x,
            n_boot=1,
            subsample_frac=p.subsample_frac,
            top_k=p.top_k,
            mode=p.mode,
            k=p.k,
            iters=p.iters,
            alpha=p.alpha,
            solver=p.solver,
            lambda_=p.lambda_,
            normalize=p.normalize,
            rescale_median=p.rescale_median,
            seed=77,
            cache_dir=cache_dir,
        )
        score_map[name] = s
        rows.append(
            {
                "method": name,
                "recall_at20": _recall_at_k(s, informative, 20),
                "support_gt_1e_2": int(np.sum(s > 1e-2)),
                "support_gt_5e_2": int(np.sum(s > 5e-2)),
                "tip_n_ok": int(out["n_ok"]),
                "life_mean": float(np.mean(out["lifetime"])) if out["lifetime"].size else 0.0,
                "prom_mean": float(np.mean(out["prominence"])) if out["prominence"].size else 0.0,
                "w_nnz": int(np.sum(np.abs(w) > 1e-8)) if w.size else 0,
            }
        )

    syn_df = pd.DataFrame(rows)
    syn_df.to_csv(tab_dir / "synthetic_metrics.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax = axes[0, 0]
    ax.scatter(z[:, 0], z[:, 1], c=(y.astype(int)), cmap="coolwarm", s=20, alpha=0.8)
    ax.set_title("Synthetic Ring Geometry")
    ax.set_xlabel("z1")
    ax.set_ylabel("z2")
    ax.axis("equal")

    ax = axes[0, 1]
    for name, color in [("harmonic", "tab:blue"), ("eu_sparse", "tab:orange"), ("ccspr", "tab:green")]:
        s = np.sort(score_map[name])[::-1][:40]
        ax.plot(np.arange(len(s)), s, marker="o", ms=2, label=name, color=color)
    ax.set_title("Feature Score Decay (Top 40)")
    ax.set_xlabel("rank")
    ax.set_ylabel("score")
    ax.legend()

    ax = axes[1, 0]
    syn_df.plot(x="method", y=["recall_at20"], kind="bar", ax=ax, legend=False, color=["#4c72b0", "#dd8452", "#55a868"])
    ax.set_ylim(0, 1)
    ax.set_title("Synthetic Informative Recall@20")
    ax.set_ylabel("recall")
    ax.set_xlabel("")

    ax = axes[1, 1]
    syn_df.plot(x="method", y=["support_gt_5e_2", "support_gt_1e_2"], kind="bar", ax=ax)
    ax.set_title("Synthetic Feature Support Size")
    ax.set_ylabel("feature count")
    ax.set_xlabel("")

    fig.suptitle("Synthetic Proof-of-Concept: Harmonic vs Sparse vs Ricci+Sparse")
    fig.tight_layout()
    fig.savefig(fig_dir / "synthetic_proof_of_concept.png", dpi=220)
    plt.close(fig)


def _load_dataset_specs(cache_dir: str) -> list[DatasetSpec]:
    luad_cfg = read_yaml("configs/tcga_luad.yaml")
    luad = load_tcga_luad(
        root=luad_cfg["dataset"].get("root", "data/tcga_luad"),
        min_class_size=int(luad_cfg["dataset"].get("min_class_size", 10)),
    )
    x_luad, y_luad = luad["X"], luad["y"]
    nmax = int(luad_cfg["dataset"].get("max_samples", 0))
    if nmax > 0:
        x_luad, y_luad = _stratified_subsample(x_luad, y_luad, nmax, int(luad_cfg["seed"]))
    luad_spec = DatasetSpec(
        name="tcga_luad",
        x=x_luad,
        y=y_luad,
        var_genes=int(luad_cfg["protocol"]["var_genes"]),
        harmonic=FeatureParams(**luad_cfg["protocol"]["harmonic_params"]),
        ccspr=FeatureParams(**luad_cfg["protocol"]["ccspr_params"]),
    )

    ven_cfg = read_yaml("configs/cll_venetoclax.yaml")
    ven = load_cll_venetoclax(
        processed_matrix_path=ven_cfg["dataset"]["processed_matrix_path"],
        label_col=ven_cfg["dataset"].get("label_col", "label"),
    )
    ven_spec = DatasetSpec(
        name="gse161711_bulk",
        x=ven["X"],
        y=ven["y"],
        var_genes=int(ven_cfg["protocol"]["var_genes"]),
        harmonic=FeatureParams(**ven_cfg["protocol"]["harmonic_params"]),
        ccspr=FeatureParams(**ven_cfg["protocol"]["ccspr_params"]),
    )

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
    scr_spec = DatasetSpec(
        name="gse165087_scrna",
        x=scr["X"],
        y=scr["y"],
        var_genes=int(scr_cfg["protocol"]["var_genes"]),
        harmonic=FeatureParams(**scr_cfg["protocol"]["harmonic_params"]),
        ccspr=FeatureParams(**scr_cfg["protocol"]["ccspr_params"]),
    )

    return [luad_spec, ven_spec, scr_spec]


def _stability_and_sparsity(specs: list[DatasetSpec], tab_dir: Path, fig_dir: Path, cache_dir: str) -> None:
    stab_rows = []
    sparse_rows = []
    luad_score_cache = {}

    for spec in specs:
        x_var, _ = select_top_variable_genes(spec.x, top_n=min(spec.var_genes, spec.x.shape[1], 500))
        for method_name, params in [("harmonic", spec.harmonic), ("ccspr", spec.ccspr)]:
            tip_sets = []
            life_vals = []
            prom_vals = []
            n_ok_vals = []
            for sd in [11, 23]:
                out = tip_bootstrap_topk(
                    x_var,
                    n_boot=1,
                    subsample_frac=float(params.subsample_frac),
                    top_k=int(params.top_k),
                    mode=params.mode,
                    k=int(params.k),
                    iters=int(params.iters),
                    alpha=float(params.alpha),
                    solver=params.solver,
                    lambda_=float(params.lambda_),
                    normalize=params.normalize,
                    rescale_median=float(params.rescale_median),
                    seed=int(sd),
                    cache_dir=cache_dir,
                )
                idx = np.argsort(-out["tip"])[: int(params.top_k)]
                tip_sets.append(set(idx.tolist()))
                life_vals.append(float(np.mean(out["lifetime"])) if out["lifetime"].size else 0.0)
                prom_vals.append(float(np.mean(out["prominence"])) if out["prominence"].size else 0.0)
                n_ok_vals.append(int(out["n_ok"]))

            stab_rows.append(
                {
                    "dataset": spec.name,
                    "method": method_name,
                    "tip_topk_jaccard_mean": _pairwise_jaccard(tip_sets),
                    "lifetime_mean": float(np.mean(life_vals)),
                    "prominence_mean": float(np.mean(prom_vals)),
                    "n_ok_mean": float(np.mean(n_ok_vals)),
                }
            )

            scores, _, _ = _single_score_vector(x_var, params, cache_dir=cache_dir)
            sparse_rows.append(
                {
                    "dataset": spec.name,
                    "method": method_name,
                    "support_gt_1e_3": int(np.sum(scores > 1e-3)),
                    "support_gt_1e_2": int(np.sum(scores > 1e-2)),
                    "support_gt_5e_2": int(np.sum(scores > 5e-2)),
                    "score_l1": float(np.sum(np.abs(scores))),
                    "score_top20_mass": float(np.sum(np.sort(scores)[::-1][:20]) / max(np.sum(scores), 1e-12)),
                }
            )
            if spec.name == "tcga_luad":
                luad_score_cache[method_name] = scores

    stab_df = pd.DataFrame(stab_rows)
    sparse_df = pd.DataFrame(sparse_rows)
    stab_df.to_csv(tab_dir / "stability_metrics.csv", index=False)
    sparse_df.to_csv(tab_dir / "sparsity_metrics.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    p = sparse_df.pivot(index="dataset", columns="method", values="support_gt_1e_2")
    p.plot(kind="bar", ax=axes[0])
    axes[0].set_title("Feature Support >1e-2")
    axes[0].set_ylabel("count")
    axes[0].set_xlabel("")
    q = stab_df.pivot(index="dataset", columns="method", values="tip_topk_jaccard_mean")
    q.plot(kind="bar", ax=axes[1], ylim=(0, 1))
    axes[1].set_title("TIP Top-K Jaccard Stability")
    axes[1].set_ylabel("Jaccard")
    axes[1].set_xlabel("")
    fig.tight_layout()
    fig.savefig(fig_dir / "harmonic_vs_ccspr_sparsity_stability.png", dpi=220)
    plt.close(fig)

    if "harmonic" in luad_score_cache and "ccspr" in luad_score_cache:
        s1 = np.sort(luad_score_cache["harmonic"])[::-1][:50]
        s2 = np.sort(luad_score_cache["ccspr"])[::-1][:50]
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(np.arange(50), s1, label="harmonic", marker="o", ms=2)
        ax.plot(np.arange(50), s2, label="ccspr", marker="o", ms=2)
        ax.set_title("LUAD Feature Support Profile: Harmonic vs CC-SPR")
        ax.set_xlabel("feature rank")
        ax.set_ylabel("score")
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "luad_feature_support_harmonic_vs_ccspr.png", dpi=220)
        plt.close(fig)


def _robustness_panel(out_path: Path) -> None:
    paths = {
        "tcga_luad": Path("results/ablation_summary.csv"),
        "gse161711_bulk": Path("results/cll_venetoclax/ablation_summary.csv"),
        "gse165087_scrna": Path("results/cll_rs_scrna/ablation_summary.csv"),
    }
    rows = []
    for name, p in paths.items():
        if p.exists():
            d = pd.read_csv(p)
            d["dataset"] = name
            rows.append(d)
    if not rows:
        return
    df = pd.concat(rows, ignore_index=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, abl in zip(axes, ["lambda", "k", "iters"]):
        sub = df[df["ablation"] == abl].copy()
        sub["ablation_value_num"] = pd.to_numeric(sub["ablation_value"], errors="coerce")
        for ds in sorted(sub["dataset"].unique()):
            s = sub[sub["dataset"] == ds].sort_values("ablation_value_num")
            ax.errorbar(s["ablation_value_num"], s["mean"], yerr=s["std"].fillna(0), marker="o", label=ds)
        if abl == "lambda":
            ax.set_xscale("log")
        ax.set_title(f"Robustness: {abl}")
        ax.set_xlabel(abl)
        ax.set_ylabel("Weighted F1")
        ax.set_ylim(0, 1)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _concept_table(tab_dir: Path) -> None:
    rows = [
        ("Metric assumption", "Euclidean distance default", "Euclidean or Ricci-corrected graph distance"),
        ("Representative criterion", "L2-minimum norm harmonic representative", "Elastic representative (L1 + lambda L2)"),
        ("Sparsity", "No explicit sparsity", "Explicit sparsity encouragement"),
        ("Stability", "Not TIP-centric by default", "TIP bootstrap top-K stability"),
        ("Feature attribution", "Can smear over many features", "Sharper support via sparse representative"),
        ("Geometry sensitivity", "Flat metric assumption", "Curvature-aware via Ollivier-Ricci flow"),
        ("Implementation dependency", "maTilDA/R-centric ecosystem", "Python stack in this repo (optional R parity hook)"),
        ("Datasets evaluated here", "LUAD baseline framing", "LUAD + BRCA + GSE161711 + GSE165087 in this repo"),
    ]
    df = pd.DataFrame(rows, columns=["dimension", "harmonic_ph", "ccspr"])
    df.to_csv(tab_dir / "harmonic_vs_ccspr_comparison.csv", index=False)


def main() -> None:
    fig_dir = ensure_dir("results/figures")
    tab_dir = ensure_dir("results/tables")
    cache_dir = "data/cache"

    _pipeline_schematic(fig_dir / "pipeline_schematic_ccspr.png")

    luad_main_path = tab_dir / "luad_main_extended_metrics.csv"
    brca_main_path = tab_dir / "brca_main_extended_metrics.csv"
    if luad_main_path.exists():
        luad_metrics = pd.read_csv(luad_main_path)
    else:
        luad_metrics = _run_extended_main_luad(cache_dir=cache_dir)
        luad_metrics.to_csv(luad_main_path, index=False)
    if brca_main_path.exists():
        brca_metrics = pd.read_csv(brca_main_path)
    else:
        brca_metrics = _run_extended_main_brca(cache_dir=cache_dir)
        brca_metrics.to_csv(brca_main_path, index=False)

    bench_rows = [
        _summary_from_metrics(luad_metrics, "tcga_luad"),
        _summary_from_metrics(brca_metrics, "tcga_brca"),
        pd.read_csv("results/cll_venetoclax/main_results.csv").assign(dataset="gse161711_bulk"),
        pd.read_csv("results/cll_rs_scrna/main_results.csv").assign(dataset="gse165087_scrna"),
    ]
    main_all = pd.concat(bench_rows, ignore_index=True)
    main_all.to_csv(tab_dir / "main_benchmark_all_datasets.csv", index=False)
    _plot_main_benchmark_all(main_all, fig_dir / "main_benchmark_all_datasets.png")

    _synthetic_bundle(fig_dir=fig_dir, tab_dir=tab_dir, cache_dir=cache_dir)

    specs = _load_dataset_specs(cache_dir=cache_dir)
    _stability_and_sparsity(specs, tab_dir=tab_dir, fig_dir=fig_dir, cache_dir=cache_dir)

    _robustness_panel(fig_dir / "robustness_lambda_k_iters_all.png")
    _concept_table(tab_dir=tab_dir)

    print("Saved:", fig_dir / "pipeline_schematic_ccspr.png")
    print("Saved:", fig_dir / "synthetic_proof_of_concept.png")
    print("Saved:", fig_dir / "main_benchmark_all_datasets.png")
    print("Saved:", fig_dir / "harmonic_vs_ccspr_sparsity_stability.png")
    print("Saved:", fig_dir / "luad_feature_support_harmonic_vs_ccspr.png")
    print("Saved:", fig_dir / "robustness_lambda_k_iters_all.png")
    print("Saved:", tab_dir / "synthetic_metrics.csv")
    print("Saved:", tab_dir / "stability_metrics.csv")
    print("Saved:", tab_dir / "sparsity_metrics.csv")
    print("Saved:", tab_dir / "harmonic_vs_ccspr_comparison.csv")
    print("Saved:", tab_dir / "main_benchmark_all_datasets.csv")


if __name__ == "__main__":
    main()
