#!/usr/bin/env python3
"""GO/KEGG enrichment analysis for Arabidopsis TIP genes.

Reads the extracted TIP genes from arabidopsis_top_features_by_backend.csv
and runs gseapy ORA against plant-relevant databases.
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

def log(msg: str):
    print(f"[{datetime.utcnow().isoformat()}Z] {msg}", flush=True)


def run_enrichment(gene_list, gene_sets, organism="plant", description=""):
    """Run gseapy ORA enrichment."""
    import gseapy as gp
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism if organism != "plant" else None,
            outdir=None,
            no_plot=True,
        )
        df = enr.results
        if df is not None and len(df) > 0:
            df = df.sort_values("Adjusted P-value")
        return df
    except Exception as e:
        log(f"  enrichment failed for {gene_sets}: {e}")
        return pd.DataFrame()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-csv",
                     default="results/arabidopsis/tables/arabidopsis_top_features_by_backend.csv")
    ap.add_argument("--out-dir", default="results/arabidopsis/enrichment")
    ap.add_argument("--top-k", type=int, default=20)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_csv = Path(args.gene_csv)
    if not gene_csv.exists():
        log(f"ERROR: {gene_csv} not found. Run gene extraction first.")
        sys.exit(1)

    gene_df = pd.read_csv(gene_csv)
    backends = gene_df["backend"].unique()
    log(f"Backends: {list(backends)}")

    # Arabidopsis gene names may be AGI codes (AT1G...) or symbols
    # Use multiple databases that accept Arabidopsis genes
    # Enrichr has limited plant support, so we try multiple approaches
    gene_set_libs = [
        "GO_Biological_Process_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023",
        "KEGG_2021_Human",  # Some gene symbols overlap
    ]

    all_results = []
    for backend in backends:
        genes = gene_df[gene_df["backend"] == backend]["gene"].tolist()[:args.top_k]
        log(f"\n{backend}: {len(genes)} genes — {genes[:5]}...")

        # Many Arabidopsis genes are AGI codes (AT1G12345) which Enrichr won't recognize.
        # Filter to gene symbols (non-AT*G* format) for Enrichr queries
        symbol_genes = [g for g in genes if not g.startswith("AT") or not "G" in g[2:6]]
        agi_genes = [g for g in genes if g not in symbol_genes]
        log(f"  Symbols: {len(symbol_genes)}, AGI codes: {len(agi_genes)}")

        if len(symbol_genes) < 3:
            log(f"  Too few symbol genes for enrichment, trying all genes...")
            query_genes = genes
        else:
            query_genes = symbol_genes

        for lib in gene_set_libs:
            log(f"  {lib}...")
            df = run_enrichment(query_genes, lib)
            if df is not None and len(df) > 0:
                df["backend"] = backend
                df["library"] = lib
                sig = df[df["Adjusted P-value"] < 0.1]
                log(f"    {len(sig)} terms with adj p < 0.1")
                if len(sig) > 0:
                    for _, row in sig.head(3).iterrows():
                        log(f"    {row['Term']}: adj_p={row['Adjusted P-value']:.4f}")
                all_results.append(df)

    if all_results:
        full_df = pd.concat(all_results, ignore_index=True)
        full_df.to_csv(out_dir / "arabidopsis_enrichment_full.csv", index=False)

        # Summary: top 5 per backend
        summary_rows = []
        for backend in backends:
            sub = full_df[full_df["backend"] == backend].sort_values("Adjusted P-value")
            summary_rows.append(sub.head(5))
        summary_df = pd.concat(summary_rows, ignore_index=True)
        summary_df.to_csv(out_dir / "arabidopsis_enrichment_summary.csv", index=False)
        log(f"\nSaved: {out_dir / 'arabidopsis_enrichment_full.csv'}")
        log(f"Saved: {out_dir / 'arabidopsis_enrichment_summary.csv'}")

        # Heatmap
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np

            sig_df = full_df[full_df["Adjusted P-value"] < 0.2].copy()
            if len(sig_df) > 0:
                top_terms = sig_df.groupby("Term")["Adjusted P-value"].min().nsmallest(15).index
                pivot = sig_df[sig_df["Term"].isin(top_terms)].pivot_table(
                    index="Term", columns="backend",
                    values="Adjusted P-value", aggfunc="min"
                ).fillna(1.0)
                fig, ax = plt.subplots(figsize=(8, max(4, len(pivot) * 0.4)))
                im = ax.imshow(-np.log10(pivot.values + 1e-10), aspect="auto",
                               cmap="YlOrRd")
                ax.set_xticks(range(len(pivot.columns)))
                ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
                ax.set_yticks(range(len(pivot.index)))
                ax.set_yticklabels(pivot.index, fontsize=8)
                plt.colorbar(im, ax=ax, label="-log10(adj p)")
                ax.set_title("Arabidopsis Root: Enrichment by Backend")
                fig.tight_layout()
                fig.savefig(out_dir / "arabidopsis_enrichment_heatmap.png", dpi=220)
                plt.close(fig)
                log(f"Saved: {out_dir / 'arabidopsis_enrichment_heatmap.png'}")
        except Exception as e:
            log(f"Heatmap generation failed: {e}")
    else:
        log("No enrichment results found.")

    log("=== DONE ===")


if __name__ == "__main__":
    main()
