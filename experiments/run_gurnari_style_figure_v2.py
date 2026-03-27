#!/usr/bin/env python3
"""Gurnari-style hierarchical clustering figure using CC-SPR cycle weights (v2).

Improved version with:
  - Ward linkage for balanced clusters
  - Side-by-side CC-SPR vs Harmonic comparison
  - Sparsity statistics overlay
  - LUAD molecular subtypes (Bronchioid/Magnoid/Squamoid)
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

_src = str(Path(__file__).resolve().parent.parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)


def log(msg: str):
    print(f"[{datetime.utcnow().isoformat()}Z] {msg}", flush=True)


def build_sample_bar_matrix(x, dist, bars, solver="elastic", lambda_=1e-6, max_bars=20):
    """Build samples x bars matrix of cycle weight participation."""
    from ccspr.solver.cycle import compute_cycle_representative
    from ccspr.topology.persistence import build_cycle_problem
    import gudhi

    n = x.shape[0]
    use_bars = bars[:max_bars]
    W = np.zeros((n, len(use_bars)), dtype=np.float64)
    bar_info = []

    for j, bar in enumerate(use_bars):
        try:
            edges, d_mat, z0 = build_cycle_problem(dist, bar)
            w = compute_cycle_representative(d_mat, z0, solver=solver, lambda_=lambda_)
            active = np.where(np.abs(w) > 0)[0]
            for idx in active:
                u, v = edges[idx]
                aw = np.abs(w[idx])
                if u < n:
                    W[u, j] += aw
                if v < n:
                    W[v, j] += aw
            bar_info.append({
                "bar_idx": j, "birth": bar["bt"], "death": bar["dt"],
                "life": bar["life"], "n_active_edges": len(active),
                "n_total_edges": len(edges), "solver": solver,
                "sparsity": 1.0 - len(active) / max(len(edges), 1),
            })
            log(f"  Bar {j+1}/{len(use_bars)}: life={bar['life']:.4f}, "
                f"active={len(active)}/{len(edges)} edges "
                f"(sparsity={bar_info[-1]['sparsity']:.1%})")
        except Exception as e:
            log(f"  Bar {j+1} failed: {e}")
            bar_info.append({
                "bar_idx": j, "birth": bar["bt"], "death": bar["dt"],
                "life": bar["life"], "n_active_edges": 0,
                "n_total_edges": 0, "solver": solver, "sparsity": 1.0,
            })
    return W, bar_info


def plot_comparison_figure(
    W_sparse, W_harm, labels, bar_info_sparse, bar_info_harm,
    out_path, dataset_name="LUAD"
):
    """Side-by-side Gurnari-style figure: CC-SPR vs Harmonic."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import pdist

    unique_labels = sorted(np.unique(labels))
    label_to_int = {l: i for i, l in enumerate(unique_labels)}
    n_classes = len(unique_labels)
    label_colors_map = ["#2196F3", "#FF5722", "#4CAF50", "#FFC107"][:n_classes]
    label_cmap = ListedColormap(label_colors_map)

    fig, axes = plt.subplots(2, 3, figsize=(18, 14),
                             gridspec_kw={"width_ratios": [0.5, 3, 1.5]})

    for row, (W, name, binfo) in enumerate([
        (W_sparse, "CC-SPR (sparse)", bar_info_sparse),
        (W_harm, "Harmonic (L2)", bar_info_harm),
    ]):
        # Filter zero-variance bars
        bar_var = np.var(W, axis=0)
        keep = bar_var > 1e-12
        W_f = W[:, keep]
        n_bars = W_f.shape[1]

        if n_bars == 0:
            for ax in axes[row]:
                ax.text(0.5, 0.5, "No data", ha="center", va="center")
            continue

        # Row-normalize
        row_max = np.max(np.abs(W_f), axis=1, keepdims=True)
        row_max[row_max < 1e-12] = 1.0
        W_n = W_f / row_max

        # Ward linkage
        Z = linkage(pdist(W_n, "euclidean"), method="ward")
        dn = dendrogram(Z, no_plot=True)
        order = np.array(dn["leaves"])
        clusters = fcluster(Z, t=n_classes, criterion="maxclust")

        W_ordered = W_n[order]
        labels_ordered = labels[order]
        clusters_ordered = clusters[order]

        # Panel 1: Dendrogram + label bar
        ax_dendro = axes[row, 0]
        dendrogram(Z, orientation="left", ax=ax_dendro,
                   no_labels=True, color_threshold=0,
                   above_threshold_color="gray")
        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])

        # Add label color strip on the right edge of dendrogram
        label_ints = np.array([label_to_int[l] for l in labels_ordered])
        # Create narrow color strip
        divider_width = 0.02
        ax_strip = ax_dendro.inset_axes([1.0, 0, divider_width * 8, 1.0])
        ax_strip.imshow(label_ints.reshape(-1, 1), aspect="auto",
                        cmap=label_cmap, interpolation="nearest")
        ax_strip.set_xticks([])
        ax_strip.set_yticks([])

        ax_dendro.set_ylabel(f"{name}\n({W.shape[0]} samples)", fontsize=11)

        # Panel 2: Heatmap
        ax_heat = axes[row, 1]
        im = ax_heat.imshow(W_ordered, aspect="auto", cmap="inferno",
                            interpolation="nearest", vmin=0)
        ax_heat.set_yticks([])
        if n_bars <= 30:
            ax_heat.set_xticks(range(n_bars))
            ax_heat.set_xticklabels([f"B{i+1}" for i in range(n_bars)],
                                    fontsize=7, rotation=45, ha="right")
        ax_heat.set_xlabel(f"H$_1$ bars (n={n_bars})", fontsize=10)

        # Sparsity annotation
        nonzero_frac = (W_f > 1e-12).mean()
        mean_sparsity = np.mean([b["sparsity"] for b in binfo if b["n_total_edges"] > 0])
        ax_heat.set_title(
            f"{name}: {nonzero_frac:.1%} nonzero entries, "
            f"mean edge sparsity {mean_sparsity:.1%}",
            fontsize=10
        )
        plt.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02, shrink=0.8)

        # Panel 3: Cluster composition bars
        ax_comp = axes[row, 2]
        for c in sorted(np.unique(clusters)):
            mask = clusters == c
            n_c = np.sum(mask)
            bottom = 0
            for i, ul in enumerate(unique_labels):
                frac = np.sum(labels[mask] == ul) / n_c
                ax_comp.barh(f"C{c}\n(n={n_c})", frac, left=bottom,
                             color=label_colors_map[i], height=0.6)
                bottom += frac

        ax_comp.set_xlim(0, 1)
        ax_comp.set_xlabel("Fraction", fontsize=10)
        ax_comp.set_title("Cluster composition", fontsize=10)

        # Purity calculation
        purity = sum(
            pd.Series(labels[clusters == c]).value_counts().max()
            for c in np.unique(clusters)
        ) / len(labels)
        ax_comp.text(0.95, 0.02, f"Purity: {purity:.2f}",
                     transform=ax_comp.transAxes, ha="right", fontsize=9,
                     bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # Legend
    legend_elements = [Patch(facecolor=label_colors_map[i], label=ul)
                       for i, ul in enumerate(unique_labels)]
    fig.legend(handles=legend_elements, loc="upper center",
               ncol=n_classes, fontsize=10, title="LUAD Subtypes",
               bbox_to_anchor=(0.5, 0.99))

    fig.suptitle(
        f"{dataset_name}: Sample-level cycle-weight heatmaps with Ward clustering\n"
        f"(cf. Gurnari et al. Fig. 4 — sparse CC-SPR vs dense harmonic representatives)",
        fontsize=13, y=1.02
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved: {out_path}")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-samples", type=int, default=150)
    ap.add_argument("--max-bars", type=int, default=15)
    ap.add_argument("--pca-dim", type=int, default=50)
    ap.add_argument("--lambda_", type=float, default=1e-6)
    ap.add_argument("--mode", default="euclidean")
    ap.add_argument("--out-dir", default="results/figures")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log("Loading TCGA-LUAD data...")
    from ccspr.datasets.tcga_luad import load_tcga_luad
    data = load_tcga_luad()
    X, y = data["X"], data["y"]
    sample_ids = data["sample_ids"]
    log(f"  {X.shape[0]} samples, labels: {dict(zip(*np.unique(y, return_counts=True)))}")

    # Balanced subsample
    np.random.seed(42)
    if X.shape[0] > args.max_samples:
        classes = np.unique(y)
        per_class = args.max_samples // len(classes)
        keep = []
        for c in classes:
            idx = np.where(y == c)[0]
            keep.extend(np.random.choice(idx, min(per_class, len(idx)), replace=False))
        keep = np.sort(keep)
        X, y, sample_ids = X[keep], y[keep], sample_ids[keep]
        log(f"  Subsampled to {len(keep)}")

    # PCA
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    X = np.nan_to_num(X, nan=0.0)
    X_pca = PCA(n_components=min(args.pca_dim, *X.shape)).fit_transform(
        StandardScaler().fit_transform(X)
    )
    log(f"  PCA → {X_pca.shape[1]} dims")

    # Distance
    from ccspr.geometry.distance import build_distance
    dist = build_distance(X_pca, mode=args.mode)

    # Persistence
    from ccspr.topology.persistence import compute_persistence, top_h1_bars
    ph = compute_persistence(dist, cache_dir="data/cache")
    bars = top_h1_bars(ph, top_k=args.max_bars)
    log(f"  {len(bars)} H1 bars")

    if not bars:
        log("No H1 bars — aborting.")
        sys.exit(1)

    # Build matrices
    log("CC-SPR (elastic) weights...")
    W_sp, bi_sp = build_sample_bar_matrix(X_pca, dist, bars, solver="elastic",
                                          lambda_=args.lambda_, max_bars=args.max_bars)
    log("Harmonic (L2) weights...")
    W_hm, bi_hm = build_sample_bar_matrix(X_pca, dist, bars, solver="l2",
                                          lambda_=1e-6, max_bars=args.max_bars)

    # Save data
    np.savez(out_dir / "gurnari_v2_matrices.npz",
             W_sparse=W_sp, W_harmonic=W_hm, labels=y, sample_ids=sample_ids)
    pd.DataFrame(bi_sp).to_csv(out_dir / "gurnari_v2_sparse_bars.csv", index=False)
    pd.DataFrame(bi_hm).to_csv(out_dir / "gurnari_v2_harmonic_bars.csv", index=False)

    # Plot
    plot_comparison_figure(
        W_sp, W_hm, y, bi_sp, bi_hm,
        out_path=out_dir / "gurnari_style_comparison.png",
        dataset_name="TCGA-LUAD",
    )

    # Print sparsity summary
    log("\n=== Sparsity comparison ===")
    sp_nz = (np.abs(W_sp) > 1e-12).mean()
    hm_nz = (np.abs(W_hm) > 1e-12).mean()
    log(f"  CC-SPR:    {sp_nz:.1%} nonzero entries")
    log(f"  Harmonic:  {hm_nz:.1%} nonzero entries")
    log(f"  Sparsity ratio: {hm_nz / max(sp_nz, 1e-12):.1f}x denser for harmonic")

    log("=== DONE ===")


if __name__ == "__main__":
    main()
