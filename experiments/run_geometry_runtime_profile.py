#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
from ccspr.datasets.cll_venetoclax import load_cll_venetoclax
from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.preprocess.basic import select_top_variable_genes
from ccspr.stability.tip import tip_bootstrap_topk
from ccspr.topology.persistence import build_cycle_problem, compute_persistence, pick_dominant_h1
from ccspr.solver.cycle import compute_cycle_representative
from ccspr.geometry.distance import build_distance
from ccspr.utils.io import ensure_dir, read_yaml


def _ci95(std: float, n: int) -> float:
    if n <= 1 or not np.isfinite(std):
        return 0.0
    return float(1.96 * std / math.sqrt(n))


def _subsample(x: np.ndarray, y: np.ndarray, n_max: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if n_max <= 0 or x.shape[0] <= n_max:
        return x, y
    rng = np.random.default_rng(seed)
    idx = rng.choice(x.shape[0], n_max, replace=False)
    return x[idx], y[idx]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/geometry_paper_ultrafast.yaml")
    ap.add_argument("--repeats", type=int, default=2)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)

    ensure_dir("results/tables")
    ensure_dir("results/figures")

    methods = [
        ("euclidean_l2", "euclidean", "l2"),
        ("euclidean_elastic", "euclidean", "elastic"),
        ("ricci_elastic", "ricci", "elastic"),
        ("diffusion_elastic", "diffusion", "elastic"),
        ("phate_elastic", "phate", "elastic"),
        ("dtm_elastic", "dtm", "elastic"),
    ]

    datasets = {}

    luad_cfg = cfg["datasets"]["tcga_luad"]
    ds_l = load_tcga_luad(root=luad_cfg.get("root", "data/tcga_luad"), min_class_size=int(luad_cfg.get("min_class_size", 10)))
    x_l, y_l = _subsample(ds_l["X"], ds_l["y"], n_max=int(luad_cfg.get("max_samples", 24)), seed=seed)
    x_l, _ = select_top_variable_genes(x_l, top_n=min(int(luad_cfg.get("var_genes", 120)), x_l.shape[1], 120))
    datasets["tcga_luad"] = (x_l, y_l)

    bcfg = cfg["datasets"].get("gse161711_bulk")
    if bcfg and bool(bcfg.get("enabled", True)):
        ds_b = load_cll_venetoclax(processed_matrix_path=bcfg.get("processed_matrix_path", "data/cll_venetoclax/processed_matrix.tsv"), label_col=bcfg.get("label_col", "label"))
        x_b, y_b = _subsample(ds_b["X"], ds_b["y"], n_max=int(bcfg.get("max_samples", 24)), seed=seed + 1)
        x_b, _ = select_top_variable_genes(x_b, top_n=min(int(bcfg.get("var_genes", 120)), x_b.shape[1], 120))
        datasets["gse161711_bulk"] = (x_b, y_b)

    scfg = cfg["datasets"].get("gse165087_scrna")
    if scfg and bool(scfg.get("enabled", True)):
        ds_s = load_cll_rs_scrna(
            data_path=scfg.get("data_path", "data/cll_rs_scrna/data.h5ad"),
            label_key=scfg.get("label_key", "time_group"),
            min_genes=int(scfg.get("min_genes", 200)),
            min_cells=int(scfg.get("min_cells", 3)),
            n_hvg=int(scfg.get("n_hvg", 2000)),
            n_pcs=int(scfg.get("n_pcs", 50)),
            max_cells=scfg.get("max_cells", 120),
            use_precomputed_pca=bool(scfg.get("use_precomputed_pca", True)),
            seed=seed,
        )
        x_s, y_s = _subsample(ds_s["X"], ds_s["y"], n_max=int(scfg.get("max_samples", 24)), seed=seed + 2)
        x_s, _ = select_top_variable_genes(x_s, top_n=min(int(scfg.get("var_genes", 50)), x_s.shape[1], 50))
        datasets["gse165087_scrna"] = (x_s, y_s)

    mcfg = cfg["methods"]
    rows = []

    for ds_name, (x, y) in datasets.items():
        for rep in range(int(args.repeats)):
            for mname, mode, solver in methods:
                t0 = time.perf_counter()
                tip_out = tip_bootstrap_topk(
                    x,
                    n_boot=1,
                    subsample_frac=float(mcfg.get("subsample_frac", 0.75)),
                    top_k=int(mcfg.get("top_k", 15)),
                    mode=mode,
                    k=int(mcfg.get("k", 8)),
                    iters=int(mcfg.get("iters", 2)),
                    alpha=float(mcfg.get("alpha", 0.5)),
                    solver=solver,
                    lambda_=float(mcfg.get("lambda_", 1e-6)),
                    normalize=str(mcfg.get("normalize", "std")),
                    rescale_median=float(mcfg.get("rescale_median", 1.0)),
                    seed=seed + 10000 * rep,
                    cache_dir=str(cfg.get("cache_dir", "data/cache")),
                )
                t1 = time.perf_counter()

                d = build_distance(
                    x,
                    mode=mode,
                    k=int(mcfg.get("k", 8)),
                    iters=int(mcfg.get("iters", 2)),
                    alpha=float(mcfg.get("alpha", 0.5)),
                    rescale_median=float(mcfg.get("rescale_median", 1.0)),
                    cache_dir=str(cfg.get("cache_dir", "data/cache")),
                )
                ph = compute_persistence(d, max_dimension=2, cache_dir=str(cfg.get("cache_dir", "data/cache")))
                bar = pick_dominant_h1(ph)
                support = 0
                if bar is not None:
                    edges, dm, z0 = build_cycle_problem(d, bar)
                    w = compute_cycle_representative(dm, z0, solver=solver, lambda_=float(mcfg.get("lambda_", 1e-6)))
                    support = int(np.sum(np.abs(np.asarray(w)) > 1e-8))
                t2 = time.perf_counter()

                rows.append(
                    {
                        "dataset": ds_name,
                        "method": mname,
                        "geometry": mode,
                        "solver": solver,
                        "repeat": int(rep),
                        "n_samples": int(x.shape[0]),
                        "n_features": int(x.shape[1]),
                        "tip_sec": float(t1 - t0),
                        "topology_cycle_sec": float(t2 - t1),
                        "total_sec": float(t2 - t0),
                        "tip_lifetime": float(np.mean(tip_out.get("lifetime", [0.0]))),
                        "support": int(support),
                    }
                )

    raw = pd.DataFrame(rows)
    raw.to_csv("results/tables/runtime_cost_geometry_raw.csv", index=False)

    summary = (
        raw.groupby(["dataset", "method", "geometry", "solver"]) [["tip_sec", "topology_cycle_sec", "total_sec", "support"]]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    summary.columns = ["_".join([c for c in col if c]) for col in summary.columns]
    summary = summary.rename(columns={"dataset_": "dataset", "method_": "method", "geometry_": "geometry", "solver_": "solver"})
    summary["total_sec_ci95"] = [_ci95(float(s), int(n)) for s, n in zip(summary["total_sec_std"].fillna(0.0), summary["total_sec_count"])]
    summary.to_csv("results/tables/runtime_cost_geometry.csv", index=False)

    overall = (
        raw.groupby(["method", "geometry", "solver"])[["tip_sec", "topology_cycle_sec", "total_sec", "support"]]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    overall.columns = ["_".join([c for c in col if c]) for col in overall.columns]
    overall = overall.rename(columns={"method_": "method", "geometry_": "geometry", "solver_": "solver"})
    overall["total_sec_ci95"] = [_ci95(float(s), int(n)) for s, n in zip(overall["total_sec_std"].fillna(0.0), overall["total_sec_count"])]
    overall.to_csv("results/tables/runtime_cost_geometry_overall.csv", index=False)

    # Figure 1: overall runtime per method
    o = overall.sort_values("total_sec_mean", ascending=False)
    fig, ax = plt.subplots(figsize=(10.8, 4.2))
    ax.bar(np.arange(len(o)), o["total_sec_mean"], yerr=o["total_sec_ci95"], capsize=3)
    ax.set_xticks(np.arange(len(o)))
    ax.set_xticklabels(o["method"], rotation=30, ha="right")
    ax.set_ylabel("Seconds (mean)")
    ax.set_title("Runtime/Cost Comparison Across Geometry Methods")
    fig.tight_layout()
    fig.savefig("results/figures/runtime_cost_geometry.png", dpi=220)
    plt.close(fig)

    # Figure 2: per-dataset runtime heatmap-like bars
    piv = summary.pivot_table(index="dataset", columns="method", values="total_sec_mean", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(12, 4.0))
    im = ax.imshow(piv.values, aspect="auto")
    ax.set_xticks(np.arange(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels(piv.index)
    ax.set_title("Runtime (sec) by Dataset and Method")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("seconds")
    fig.tight_layout()
    fig.savefig("results/figures/runtime_cost_geometry_heatmap.png", dpi=220)
    plt.close(fig)

    print("Saved runtime/cost tables and figures")


if __name__ == "__main__":
    main()
