#!/usr/bin/env python3
"""Gurnari-style hierarchical clustering figure using CC-SPR cycle weights.

Analogous to Gurnari et al. Fig. 4:
  (a) Heatmap of samples × H1-bars cycle-weight matrix with single-linkage
      hierarchical clustering on samples (rows).
  (b) Distribution of descriptors (tumor/normal, TIP scores) within the
      detected cluster.

Uses LUAD TCGA data with CC-SPR sparse cycle representatives.
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# -- project path setup --
_src = str(Path(__file__).resolve().parent.parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)


def log(msg: str):
    print(f"[{datetime.utcnow().isoformat()}Z] {msg}", flush=True)


def build_sample_bar_matrix(
    x: np.ndarray,
    dist: np.ndarray,
    bars: list,
    solver: str = "elastic",
    lambda_: float = 1e-6,
    max_bars: int = 20,
) -> tuple[np.ndarray, list]:
    """Build samples × bars matrix of cycle weight participation.

    For each H1 bar, compute cycle weights and sum absolute edge weights
    incident to each sample. Returns (n_samples, n_bars) matrix.
    """
    from ccspr.solver.cycle import compute_cycle_representative
    from ccspr.topology.persistence import build_cycle_problem

    import gudhi

    n = x.shape[0]
    use_bars = bars[:max_bars]
    W = np.zeros((n, len(use_bars)), dtype=np.float64)
    bar_info = []

    rips = gudhi.RipsComplex(distance_matrix=dist)
    st = rips.create_simplex_tree(max_dimension=2)

    for j, bar in enumerate(use_bars):
        log(f"  Bar {j+1}/{len(use_bars)}: life={bar['life']:.4f}")
        try:
            edges, d_mat, z0 = build_cycle_problem(dist, bar)
            w = compute_cycle_representative(d_mat, z0, solver=solver, lambda_=lambda_)

            # Sum |w| for edges incident to each sample
            active = np.where(np.abs(w) > 0)[0]
            for idx in active:
                u, v = edges[idx]
                aw = np.abs(w[idx])
                if u < n:
                    W[u, j] += aw
                if v < n:
                    W[v, j] += aw

            bar_info.append({
                "bar_idx": j,
                "birth": bar["bt"],
                "death": bar["dt"],
                "life": bar["life"],
                "n_active_edges": len(active),
                "solver": solver,
            })
        except Exception as e:
            log(f"    Failed: {e}")
            bar_info.append({
                "bar_idx": j,
                "birth": bar["bt"],
                "death": bar["dt"],
                "life": bar["life"],
                "n_active_edges": 0,
                "solver": solver,
            })

    return W, bar_info


def plot_gurnari_figure(
    W: np.ndarray,
    labels: np.ndarray,
    sample_ids: np.ndarray,
    bar_info: list,
    out_path: Path,
    title: str = "CC-SPR",
):
    """Plot Gurnari Fig. 4-style panel: heatmap + descriptor distributions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import pdist

    # Remove zero-variance bars
    bar_var = np.var(W, axis=0)
    keep_bars = bar_var > 1e-12
    W_filt = W[:, keep_bars]
    n_bars_used = W_filt.shape[1]
    log(f"  Using {n_bars_used} non-trivial bars out of {W.shape[1]}")

    if n_bars_used == 0:
        log("  No non-trivial bars — skipping figure.")
        return

    # Normalize rows for visualization
    row_max = np.max(np.abs(W_filt), axis=1, keepdims=True)
    row_max[row_max < 1e-12] = 1.0
    W_norm = W_filt / row_max

    # Hierarchical clustering
    dist_vec = pdist(W_norm, metric="euclidean")
    Z = linkage(dist_vec, method="single")

    # Determine clusters (2 clusters to match binary labels)
    n_clusters = min(len(np.unique(labels)), 4)
    clusters = fcluster(Z, t=n_clusters, criterion="maxclust")

    # Reorder samples by dendrogram
    dn = dendrogram(Z, no_plot=True)
    order = np.array(dn["leaves"])

    W_ordered = W_norm[order]
    labels_ordered = labels[order]
    clusters_ordered = clusters[order]

    # -- Create figure --
    fig = plt.figure(figsize=(14, 8))

    # Panel (a): Heatmap with dendrogram
    ax_dendro = fig.add_axes([0.02, 0.1, 0.08, 0.8])
    dn_plot = dendrogram(Z, orientation="left", ax=ax_dendro,
                         no_labels=True, color_threshold=0,
                         above_threshold_color="gray")
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])
    ax_dendro.spines[:].set_visible(False)

    # Color bar for labels
    ax_labels = fig.add_axes([0.11, 0.1, 0.02, 0.8])
    unique_labels = sorted(np.unique(labels))
    label_to_int = {l: i for i, l in enumerate(unique_labels)}
    label_colors = np.array([label_to_int[l] for l in labels_ordered]).reshape(-1, 1)

    label_cmap = ListedColormap(["#2196F3", "#FF5722", "#4CAF50", "#FFC107"][:len(unique_labels)])
    ax_labels.imshow(label_colors, aspect="auto", cmap=label_cmap,
                     interpolation="nearest")
    ax_labels.set_xticks([])
    ax_labels.set_yticks([])
    ax_labels.set_xlabel("Label", fontsize=8, rotation=0)

    # Cluster color bar
    ax_clust = fig.add_axes([0.135, 0.1, 0.02, 0.8])
    clust_colors = clusters_ordered.reshape(-1, 1)
    clust_cmap = ListedColormap(["#FF9800", "#9C27B0", "#009688", "#E91E63"][:n_clusters])
    ax_clust.imshow(clust_colors, aspect="auto", cmap=clust_cmap,
                    interpolation="nearest")
    ax_clust.set_xticks([])
    ax_clust.set_yticks([])
    ax_clust.set_xlabel("Clust", fontsize=8, rotation=0)

    # Heatmap
    ax_heat = fig.add_axes([0.17, 0.1, 0.45, 0.8])
    im = ax_heat.imshow(W_ordered, aspect="auto", cmap="viridis",
                        interpolation="nearest")
    ax_heat.set_xlabel(f"H$_1$ bars (n={n_bars_used})", fontsize=11)
    ax_heat.set_ylabel(f"Samples (n={W.shape[0]})", fontsize=11)
    ax_heat.set_yticks([])
    bar_labels = [f"B{i+1}" for i in range(n_bars_used)]
    if n_bars_used <= 30:
        ax_heat.set_xticks(range(n_bars_used))
        ax_heat.set_xticklabels(bar_labels, fontsize=7, rotation=45, ha="right")
    ax_heat.set_title(f"(a) {title} cycle-weight heatmap\nwith single-linkage clustering",
                      fontsize=12)
    plt.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02, label="Normalized weight")

    # Panel (b): Descriptor distributions per cluster
    ax_dist = fig.add_axes([0.72, 0.1, 0.26, 0.8])

    # For each cluster, show label composition
    cluster_data = []
    for c in sorted(np.unique(clusters)):
        mask = clusters == c
        for ul in unique_labels:
            count = np.sum(labels[mask] == ul)
            frac = count / np.sum(mask) if np.sum(mask) > 0 else 0
            cluster_data.append({
                "Cluster": f"C{c} (n={np.sum(mask)})",
                "Label": ul,
                "Count": count,
                "Fraction": frac,
            })

    cdf = pd.DataFrame(cluster_data)
    clusters_unique = cdf["Cluster"].unique()
    x_pos = np.arange(len(clusters_unique))
    bar_width = 0.8 / len(unique_labels)
    colors = ["#2196F3", "#FF5722", "#4CAF50", "#FFC107"][:len(unique_labels)]

    for i, ul in enumerate(unique_labels):
        sub = cdf[cdf["Label"] == ul]
        fracs = [sub[sub["Cluster"] == c]["Fraction"].values[0]
                 if len(sub[sub["Cluster"] == c]) > 0 else 0
                 for c in clusters_unique]
        ax_dist.bar(x_pos + i * bar_width, fracs, bar_width,
                    label=ul, color=colors[i], alpha=0.85)

    ax_dist.set_xticks(x_pos + bar_width * (len(unique_labels) - 1) / 2)
    ax_dist.set_xticklabels(clusters_unique, fontsize=9)
    ax_dist.set_ylabel("Fraction", fontsize=11)
    ax_dist.set_title("(b) Label composition\nper cluster", fontsize=12)
    ax_dist.legend(fontsize=8, loc="upper right")
    ax_dist.set_ylim(0, 1.05)

    # Legend for label colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors[i], label=ul)
                       for i, ul in enumerate(unique_labels)]

    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    log(f"  Saved: {out_path}")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-samples", type=int, default=120,
                    help="Max samples to use (for tractable computation)")
    ap.add_argument("--max-bars", type=int, default=15,
                    help="Max H1 bars to compute cycles for")
    ap.add_argument("--pca-dim", type=int, default=50,
                    help="PCA dimension for distance computation")
    ap.add_argument("--solver", default="elastic",
                    choices=["elastic", "l2"])
    ap.add_argument("--lambda_", type=float, default=1e-6)
    ap.add_argument("--mode", default="euclidean",
                    choices=["euclidean", "ricci", "dtm", "diffusion", "phate_like"])
    ap.add_argument("--out-dir", default="results/figures")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load LUAD data
    log("Loading TCGA-LUAD data...")
    from ccspr.datasets.tcga_luad import load_tcga_luad
    data = load_tcga_luad()
    X, y = data["X"], data["y"]
    sample_ids = data["sample_ids"]
    feature_names = data["feature_names"]
    log(f"  Loaded: {X.shape[0]} samples × {X.shape[1]} features")
    log(f"  Labels: {dict(zip(*np.unique(y, return_counts=True)))}")

    # 2. Subsample if needed (balanced)
    np.random.seed(42)
    if X.shape[0] > args.max_samples:
        classes = np.unique(y)
        per_class = args.max_samples // len(classes)
        keep_idx = []
        for c in classes:
            c_idx = np.where(y == c)[0]
            sel = np.random.choice(c_idx, size=min(per_class, len(c_idx)), replace=False)
            keep_idx.extend(sel)
        keep_idx = np.sort(keep_idx)
        X = X[keep_idx]
        y = y[keep_idx]
        sample_ids = sample_ids[keep_idx]
        log(f"  Subsampled to {len(keep_idx)} samples (balanced)")

    # 3. PCA reduction
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Handle NaN
    X = np.nan_to_num(X, nan=0.0)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=min(args.pca_dim, X.shape[1], X.shape[0]))
    X_pca = pca.fit_transform(X_scaled)
    log(f"  PCA: {X_pca.shape[1]} components, {pca.explained_variance_ratio_.sum():.2%} variance")

    # 4. Build distance matrix
    log(f"Building distance matrix (mode={args.mode})...")
    from ccspr.geometry.distance import build_distance
    dist = build_distance(X_pca, mode=args.mode)
    log(f"  Distance matrix: {dist.shape}")

    # 5. Compute persistence
    log("Computing persistence...")
    from ccspr.topology.persistence import compute_persistence, top_h1_bars
    ph = compute_persistence(dist, cache_dir="data/cache")
    bars = top_h1_bars(ph, top_k=args.max_bars)
    top_lives = [round(b["life"], 4) for b in bars[:5]]
    log(f"  Found {len(bars)} H1 bars (top lifetimes: {top_lives})")

    if len(bars) == 0:
        log("ERROR: No H1 bars found. Cannot create figure.")
        sys.exit(1)

    # 6. Build samples × bars cycle weight matrix
    log(f"Building sample-bar weight matrix ({args.solver} solver)...")
    W_sparse, bar_info_sparse = build_sample_bar_matrix(
        X_pca, dist, bars,
        solver=args.solver, lambda_=args.lambda_, max_bars=args.max_bars,
    )
    log(f"  Sparse matrix: {W_sparse.shape}, nonzero frac: {(W_sparse > 0).mean():.3f}")

    # Also compute harmonic (l2) for comparison
    log("Building harmonic (L2) weight matrix for comparison...")
    W_harm, bar_info_harm = build_sample_bar_matrix(
        X_pca, dist, bars,
        solver="l2", lambda_=1e-6, max_bars=args.max_bars,
    )
    log(f"  Harmonic matrix: {W_harm.shape}, nonzero frac: {(W_harm > 0).mean():.3f}")

    # 7. Save matrices
    np.savez(
        out_dir / "gurnari_style_matrices.npz",
        W_sparse=W_sparse, W_harmonic=W_harm,
        labels=y, sample_ids=sample_ids,
    )
    pd.DataFrame(bar_info_sparse).to_csv(
        out_dir / "gurnari_style_bar_info.csv", index=False
    )
    log(f"  Saved matrices to {out_dir}")

    # 8. Plot CC-SPR (sparse) figure
    log("Plotting CC-SPR (sparse) figure...")
    plot_gurnari_figure(
        W_sparse, y, sample_ids, bar_info_sparse,
        out_path=out_dir / "gurnari_style_ccspr.png",
        title="CC-SPR (sparse)",
    )

    # 9. Plot harmonic figure for comparison
    log("Plotting harmonic figure...")
    plot_gurnari_figure(
        W_harm, y, sample_ids, bar_info_harm,
        out_path=out_dir / "gurnari_style_harmonic.png",
        title="Harmonic (L2)",
    )

    # 10. Summary statistics
    log("\n=== Cluster purity summary ===")
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import pdist

    for name, W in [("CC-SPR", W_sparse), ("Harmonic", W_harm)]:
        bar_var = np.var(W, axis=0)
        W_f = W[:, bar_var > 1e-12]
        if W_f.shape[1] == 0:
            log(f"  {name}: no non-trivial bars")
            continue
        row_max = np.max(np.abs(W_f), axis=1, keepdims=True)
        row_max[row_max < 1e-12] = 1.0
        W_n = W_f / row_max
        Z = linkage(pdist(W_n, "euclidean"), "single")
        for nc in [2, 3]:
            clusters = fcluster(Z, t=nc, criterion="maxclust")
            # Purity
            purity = 0
            for c in np.unique(clusters):
                mask = clusters == c
                label_counts = pd.Series(y[mask]).value_counts()
                purity += label_counts.max()
            purity /= len(y)
            log(f"  {name} ({nc} clusters): purity = {purity:.3f}")

    log("=== DONE ===")


if __name__ == "__main__":
    main()
