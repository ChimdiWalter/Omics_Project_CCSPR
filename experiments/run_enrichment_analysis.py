#!/usr/bin/env python3
"""Pathway enrichment analysis for CC-SPR TIP-selected gene sets.

Uses gseapy ORA (over-representation analysis) against MSigDB Hallmark,
KEGG, and Reactome gene sets for human datasets, and GO Biological Process
for Arabidopsis.

Reads existing TIP gene lists from results/tables/luad_top_features_by_method.csv
or recomputes if directed. Outputs enrichment tables and heatmap figures.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def load_tip_genes(csv_path: str | Path) -> dict[str, list[str]]:
    """Load per-method gene lists from top_features CSV."""
    df = pd.read_csv(csv_path)
    result = {}
    for method, grp in df.groupby("method"):
        genes = grp.sort_values("rank")["gene"].astype(str).tolist()
        # Filter out non-standard gene names
        genes = [g for g in genes if not g.startswith("?|") and g != "nan"]
        result[str(method)] = genes
    return result


def load_gene_universe(dataset: str = "luad") -> list[str]:
    """Load full gene universe for background in ORA.

    For LUAD: loads the full TCGA-LUAD expression matrix feature names.
    """
    if dataset == "luad":
        from ccspr.datasets.tcga_luad import load_tcga_luad
        ds = load_tcga_luad(root="data/tcga_luad", min_class_size=2)
        return ds["feature_names"].astype(str).tolist()
    elif dataset == "arabidopsis":
        # Will be populated when Arabidopsis loader is ready
        try:
            from ccspr.datasets.arabidopsis_root import load_arabidopsis_root
            ds = load_arabidopsis_root(data_path="data/arabidopsis_root")
            return ds["feature_names"].astype(str).tolist()
        except ImportError:
            return []
    return []


def run_gseapy_ora(
    gene_list: list[str],
    gene_sets: str | list[str],
    background: list[str] | int | None = None,
    organism: str = "human",
    cutoff: float = 0.25,
) -> pd.DataFrame:
    """Run gseapy ORA enrichment.

    Parameters
    ----------
    gene_list : list of gene symbols
    gene_sets : MSigDB name(s) or path to GMT, e.g.
        'MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022'
    background : gene universe (list of symbols or int count)
    organism : 'human' or 'arabidopsis' (for Enrichr libraries)
    cutoff : adjusted p-value cutoff for reporting
    """
    import gseapy as gp

    if isinstance(gene_sets, str):
        gene_sets = [gene_sets]

    all_results = []
    for gs in gene_sets:
        try:
            enr = gp.enrich(
                gene_list=gene_list,
                gene_sets=gs,
                background=background,
                outdir=None,
                no_plot=True,
                cutoff=cutoff,
                verbose=False,
            )
            if enr.results is not None and len(enr.results) > 0:
                df = enr.results.copy()
                df["gene_set_library"] = gs
                all_results.append(df)
        except Exception as e:
            print(f"  Warning: {gs} failed: {e}", flush=True)

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()


def run_luad_enrichment(
    tip_genes: dict[str, list[str]],
    background: list[str] | None = None,
    out_dir: Path = Path("results/enrichment"),
) -> pd.DataFrame:
    """Run enrichment for all LUAD methods."""
    gene_set_libs = [
        "MSigDB_Hallmark_2020",
        "KEGG_2021_Human",
        "Reactome_2022",
    ]

    all_rows = []
    for method, genes in tip_genes.items():
        print(f"  Enrichment: {method} ({len(genes)} genes)...", flush=True)
        df = run_gseapy_ora(
            gene_list=genes,
            gene_sets=gene_set_libs,
            background=background,
            organism="human",
            cutoff=0.5,  # Keep more results, filter later
        )
        if not df.empty:
            df["method"] = method
            all_rows.append(df)
        else:
            print(f"    No enrichment results for {method}", flush=True)

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        combined.to_csv(out_dir / "luad_enrichment_full.csv", index=False)
        return combined
    return pd.DataFrame()


def make_enrichment_summary(
    df: pd.DataFrame,
    top_n: int = 5,
) -> pd.DataFrame:
    """Extract top N enriched terms per method."""
    if df.empty:
        return pd.DataFrame()

    rows = []
    for method, grp in df.groupby("method"):
        grp_sorted = grp.sort_values("Adjusted P-value", ascending=True)
        for _, row in grp_sorted.head(top_n).iterrows():
            rows.append({
                "method": method,
                "term": row.get("Term", ""),
                "gene_set_library": row.get("gene_set_library", ""),
                "overlap": row.get("Overlap", ""),
                "p_value": row.get("P-value", 1.0),
                "adjusted_p": row.get("Adjusted P-value", 1.0),
                "neglog10_padj": float(-np.log10(max(1e-300, row.get("Adjusted P-value", 1.0)))),
                "genes": row.get("Genes", ""),
            })
    return pd.DataFrame(rows)


def plot_enrichment_heatmap(
    summary_df: pd.DataFrame,
    out_path: Path,
    title: str = "Enrichment: -log10(adj p)",
    max_terms: int = 15,
):
    """Plot enrichment heatmap: methods × terms."""
    if summary_df.empty:
        print("  No enrichment data to plot.", flush=True)
        return

    # Pivot: get top terms by minimum p across methods
    term_min_p = summary_df.groupby("term")["adjusted_p"].min().sort_values()
    top_terms = term_min_p.head(max_terms).index.tolist()

    sub = summary_df[summary_df["term"].isin(top_terms)].copy()
    pivot = sub.pivot_table(
        index="term", columns="method", values="neglog10_padj", fill_value=0.0,
    )
    # Reorder terms by significance
    pivot = pivot.reindex(top_terms)

    fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 1.2), max(5, len(pivot) * 0.4)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(pivot.index)))
    # Truncate long term names
    ylabels = [t[:50] + "..." if len(t) > 50 else t for t in pivot.index]
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, label="-log10(adj p)", shrink=0.7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}", flush=True)


def plot_backend_gene_overlap(
    tip_genes: dict[str, list[str]],
    out_path: Path,
):
    """Plot gene overlap across backends as a matrix heatmap."""
    methods = sorted(tip_genes.keys())
    n = len(methods)
    overlap_matrix = np.zeros((n, n), dtype=int)
    for i, m1 in enumerate(methods):
        for j, m2 in enumerate(methods):
            overlap_matrix[i, j] = len(set(tip_genes[m1]) & set(tip_genes[m2]))

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(overlap_matrix, cmap="Blues", interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_xticklabels(methods, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(methods, fontsize=8)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, str(overlap_matrix[i, j]),
                    ha="center", va="center", fontsize=8,
                    color="white" if overlap_matrix[i, j] > n * 0.6 else "black")
    ax.set_title("Gene Overlap Across Backends (Top-20 TIP Genes)")
    plt.colorbar(im, ax=ax, label="Shared genes", shrink=0.7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}", flush=True)


def compute_shared_vs_specific(
    tip_genes: dict[str, list[str]],
) -> tuple[list[str], dict[str, list[str]]]:
    """Identify genes shared across all backends vs backend-specific."""
    all_sets = {m: set(g) for m, g in tip_genes.items()}
    if not all_sets:
        return [], {}

    shared = set.intersection(*all_sets.values())
    specific = {}
    for m, gs in all_sets.items():
        others = set.union(*(s for k, s in all_sets.items() if k != m))
        specific[m] = sorted(gs - others)

    return sorted(shared), specific


def main():
    ap = argparse.ArgumentParser(description="Enrichment analysis for CC-SPR TIP genes")
    ap.add_argument("--tip-csv", default="results/tables/luad_top_features_by_method.csv",
                     help="Path to TIP gene list CSV")
    ap.add_argument("--dataset", default="luad", choices=["luad", "arabidopsis"])
    ap.add_argument("--out-dir", default="results/enrichment")
    ap.add_argument("--no-background", action="store_true",
                     help="Skip loading full gene universe (use default background)")
    ap.add_argument("--top-summary", type=int, default=5,
                     help="Top N enriched terms per method in summary")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)

    # 1. Load TIP genes
    print(f"Loading TIP genes from {args.tip_csv}...", flush=True)
    tip_genes = load_tip_genes(args.tip_csv)
    print(f"  Methods: {list(tip_genes.keys())}", flush=True)
    for m, g in tip_genes.items():
        print(f"    {m}: {len(g)} genes — {g[:5]}...", flush=True)

    # 2. Load gene universe for background
    background = None
    if not args.no_background:
        print(f"Loading gene universe for {args.dataset}...", flush=True)
        bg = load_gene_universe(args.dataset)
        if bg:
            background = bg
            print(f"  Universe: {len(bg)} genes", flush=True)
        else:
            print("  Warning: Could not load gene universe, using default", flush=True)

    # 3. Run enrichment
    if args.dataset == "luad":
        print("Running LUAD enrichment...", flush=True)
        full_df = run_luad_enrichment(tip_genes, background=background, out_dir=out_dir)
    else:
        print(f"Enrichment for {args.dataset} not yet implemented", flush=True)
        return

    # 4. Summary table
    summary = make_enrichment_summary(full_df, top_n=args.top_summary)
    if not summary.empty:
        summary.to_csv(out_dir / "luad_enrichment_summary.csv", index=False)
        print(f"\nEnrichment Summary (top {args.top_summary} per method):")
        print(summary[["method", "term", "adjusted_p", "overlap"]].to_string(index=False))
    else:
        print("\nNo significant enrichment found across any method.", flush=True)

    # 5. Enrichment heatmap
    plot_enrichment_heatmap(
        summary, out_dir / "luad_enrichment_heatmap.png",
        title="LUAD: CC-SPR Backend Enrichment (-log10 adj p)",
    )

    # 6. Gene overlap analysis
    print("\nBackend gene overlap analysis...", flush=True)
    plot_backend_gene_overlap(tip_genes, out_dir / "luad_backend_gene_overlap.png")

    shared, specific = compute_shared_vs_specific(tip_genes)
    print(f"  Shared across ALL backends: {shared}")
    for m, sp in specific.items():
        print(f"  Unique to {m}: {sp}")

    overlap_rows = []
    overlap_rows.append({"category": "shared_all", "method": "all", "genes": ";".join(shared), "count": len(shared)})
    for m, sp in specific.items():
        overlap_rows.append({"category": "unique", "method": m, "genes": ";".join(sp), "count": len(sp)})
    pd.DataFrame(overlap_rows).to_csv(out_dir / "luad_backend_gene_overlap.csv", index=False)

    # 7. Run enrichment on SHARED genes if enough
    if len(shared) >= 3:
        print(f"\nRunning enrichment on {len(shared)} shared genes...", flush=True)
        shared_enr = run_gseapy_ora(
            shared,
            ["MSigDB_Hallmark_2020", "KEGG_2021_Human"],
            background=background,
            cutoff=0.5,
        )
        if not shared_enr.empty:
            shared_enr.to_csv(out_dir / "luad_shared_genes_enrichment.csv", index=False)
            print("  Shared gene enrichment:")
            print(shared_enr[["Term", "Adjusted P-value", "Overlap"]].head(10).to_string(index=False))

    print("\n=== Enrichment analysis complete ===", flush=True)


if __name__ == "__main__":
    main()
