#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, norm

from ccspr.datasets.common import normalize_sample_id, read_tsv_auto
from ccspr.datasets.tcga_luad import load_tcga_luad
from ccspr.preprocess.basic import select_top_variable_genes, standardize
from ccspr.stability.tip import tip_bootstrap_topk
from ccspr.utils.io import ensure_dir, read_yaml


def _stratified_subsample(x: np.ndarray, y: np.ndarray, sample_ids: np.ndarray, n_max: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if n_max <= 0 or x.shape[0] <= n_max:
        return x, y, sample_ids
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
    keep = np.asarray(sorted(keep), dtype=int)
    return x[keep], y[keep], sample_ids[keep]


def _pathway_sets() -> dict[str, set[str]]:
    return {
        "LUNG_ALVEOLAR_AT2": {
            "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "NAPSA", "CLDN18", "EPCAM", "KRT8", "KRT18",
        },
        "CLUB_SECRETORY": {
            "SCGB1A1", "SCGB3A1", "KRT19", "EPCAM", "KRT8", "KRT18", "MUC1",
        },
        "MUCINOUS_LINEAGE": {
            "MUC1", "MUC5B", "MUC5AC", "SPINK1", "AGR2", "ERN2", "TFF1", "TFF3", "KRT7",
        },
        "FIBRINOGEN_COAG": {
            "FGA", "FGB", "FGG", "SERPINA1", "APOH", "PLG", "F2", "PROC",
        },
        "SEX_CHR_MARKERS": {
            "RPS4Y1", "DDX3Y", "KDM5D", "USP9Y", "EIF1AY", "UTY", "XIST",
        },
        "EGFR_ERBB": {
            "EGFR", "ERBB2", "ERBB3", "ERBB4", "GRB2", "SHC1", "SOS1", "PIK3CA", "AKT1", "PTEN", "KRAS", "BRAF", "MAP2K1", "MAPK1",
        },
        "KRAS_MAPK": {
            "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "DUSP6", "ELK1",
        },
        "PI3K_AKT_MTOR": {
            "PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "MTOR", "PTEN", "TSC1", "TSC2", "RHEB",
        },
        "CELL_CYCLE": {
            "CDK1", "CDK2", "CDK4", "CDK6", "CCNA2", "CCNB1", "CCNE1", "E2F1", "RB1", "MKI67", "PCNA",
        },
        "EMT": {
            "VIM", "CDH1", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "ITGA5", "COL1A1",
        },
        "DNA_REPAIR": {
            "BRCA1", "BRCA2", "RAD51", "ATM", "ATR", "CHEK1", "CHEK2", "TP53", "PARP1", "XRCC1", "MRE11",
        },
    }


def _top_idx(scores: np.ndarray, k: int) -> np.ndarray:
    k = min(int(k), scores.shape[0])
    if k <= 0:
        return np.array([], dtype=int)
    return np.argsort(-scores)[:k]


def _coherence_for_method(top_genes: list[str], universe_genes: list[str], pathways: dict[str, set[str]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = []
    U = {g.upper() for g in universe_genes}
    T = {g.upper() for g in top_genes}
    M = len(U)
    N = len(T)
    if M == 0 or N == 0:
        return rows, {"best_pathway": "NA", "best_neglog10_p": 0.0, "mean_jaccard": 0.0, "max_overlap": 0}

    jac_vals = []
    best = ("NA", 1.0, 0)

    for name, gset in pathways.items():
        G = U & {g.upper() for g in gset}
        n = len(G)
        k = len(T & G)
        p = float(hypergeom.sf(k - 1, M, n, N)) if n > 0 and k > 0 else 1.0
        jac = float(k / max(1, len(T | G)))
        jac_vals.append(jac)
        if p < best[1]:
            best = (name, p, k)
        rows.append(
            {
                "pathway": name,
                "universe_hits": n,
                "overlap": k,
                "jaccard": jac,
                "hypergeom_p": p,
                "neglog10_p": float(-np.log10(max(1e-300, p))),
            }
        )

    summary = {
        "best_pathway": best[0],
        "best_neglog10_p": float(-np.log10(max(1e-300, best[1]))),
        "mean_jaccard": float(np.mean(jac_vals)) if jac_vals else 0.0,
        "max_overlap": int(best[2]),
    }
    return rows, summary


def _km_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uniq = np.unique(time[event == 1])
    if uniq.size == 0:
        return np.array([0.0]), np.array([1.0])
    surv = [1.0]
    tvals = [0.0]
    s = 1.0
    for t in np.sort(uniq):
        at_risk = np.sum(time >= t)
        d = np.sum((time == t) & (event == 1))
        if at_risk > 0:
            s *= (1.0 - d / at_risk)
        tvals.append(float(t))
        surv.append(float(s))
    return np.asarray(tvals), np.asarray(surv)


def _median_survival(time: np.ndarray, event: np.ndarray) -> float:
    t, s = _km_curve(time, event)
    idx = np.where(s <= 0.5)[0]
    if idx.size == 0:
        return float(np.nan)
    return float(t[idx[0]])


def _logrank(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> tuple[float, float]:
    times = np.unique(time[event == 1])
    o_minus_e = 0.0
    var = 0.0
    for t in times:
        r1 = np.sum((group == 1) & (time >= t))
        r0 = np.sum((group == 0) & (time >= t))
        r = r1 + r0
        if r <= 1:
            continue
        d1 = np.sum((group == 1) & (time == t) & (event == 1))
        d0 = np.sum((group == 0) & (time == t) & (event == 1))
        d = d1 + d0
        if d == 0:
            continue
        e1 = d * (r1 / r)
        v1 = (r1 * r0 * d * (r - d)) / (r * r * (r - 1))
        o_minus_e += (d1 - e1)
        var += v1
    if var <= 1e-12:
        return 0.0, 1.0
    z = float(o_minus_e / math.sqrt(var))
    p = float(2.0 * (1.0 - norm.cdf(abs(z))))
    return z, p


def _c_index(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    conc = 0.0
    ties = 0.0
    comp = 0.0
    n = len(time)
    for i in range(n):
        for j in range(i + 1, n):
            ti, tj = time[i], time[j]
            ei, ej = event[i], event[j]
            ri, rj = risk[i], risk[j]

            if ei == 1 and ti < tj:
                comp += 1
                if ri > rj:
                    conc += 1
                elif ri == rj:
                    ties += 1
            elif ej == 1 and tj < ti:
                comp += 1
                if rj > ri:
                    conc += 1
                elif ri == rj:
                    ties += 1
    if comp == 0:
        return float("nan")
    return float((conc + 0.5 * ties) / comp)


def _parse_survival(labels_path: str | Path) -> pd.DataFrame:
    df = read_tsv_auto(labels_path)
    id_col = "sampleID" if "sampleID" in df.columns else df.columns[0]

    out = pd.DataFrame({"sample_id": df[id_col].astype(str).map(normalize_sample_id)})

    if "days_to_death" in df.columns:
        d_death = pd.to_numeric(df["days_to_death"], errors="coerce")
    else:
        d_death = pd.Series(np.nan, index=df.index)

    if "days_to_last_followup" in df.columns:
        d_last = pd.to_numeric(df["days_to_last_followup"], errors="coerce")
    else:
        d_last = pd.Series(np.nan, index=df.index)

    vital = df["vital_status"].astype(str).str.lower() if "vital_status" in df.columns else pd.Series("", index=df.index)

    event = (vital.str.contains("dead|deceased", regex=True)).astype(int)
    time = d_death.where(event == 1, d_last)
    time = time.where(np.isfinite(time), d_death.fillna(d_last))

    out["time_days"] = pd.to_numeric(time, errors="coerce")
    out["event"] = event.astype(int)
    out = out.dropna(subset=["time_days"]).drop_duplicates(subset=["sample_id"])
    out = out[out["time_days"] > 0]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/geometry_paper_ultrafast.yaml")
    ap.add_argument("--top-k", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=5)
    ap.add_argument("--max-samples", type=int, default=80)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    seed = int(cfg.get("seed", 42))
    np.random.seed(seed)

    ensure_dir("results/tables")
    ensure_dir("results/figures")

    dcfg = cfg["datasets"]["tcga_luad"]
    ds = load_tcga_luad(root=dcfg.get("root", "data/tcga_luad"), min_class_size=int(dcfg.get("min_class_size", 10)))

    x = ds["X"]
    y = ds["y"]
    sample_ids = ds["sample_ids"]
    feat = ds["feature_names"].astype(str)
    x, y, sample_ids = _stratified_subsample(x, y, sample_ids, n_max=int(args.max_samples), seed=seed)

    top_n = min(int(dcfg.get("var_genes", 500)), x.shape[1], 500)
    x_var, idx_var = select_top_variable_genes(x, top_n=top_n)
    genes_var = feat[idx_var]

    methods = [
        ("euclidean_l2", "euclidean", "l2"),
        ("euclidean_elastic", "euclidean", "elastic"),
        ("ricci_elastic", "ricci", "elastic"),
        ("diffusion_elastic", "diffusion", "elastic"),
        ("phate_elastic", "phate", "elastic"),
        ("dtm_elastic", "dtm", "elastic"),
    ]

    mcfg = cfg["methods"]
    top_gene_rows = []
    pathway_rows = []
    pathway_summary_rows = []

    tip_cache: dict[str, np.ndarray] = {}
    top_genes_cache: dict[str, list[str]] = {}

    pathways = _pathway_sets()

    for i, (name, mode, solver) in enumerate(methods):
        print(f"Running LUAD bio diagnostics for {name} on n={x_var.shape[0]}...", flush=True)
        tip_out = tip_bootstrap_topk(
            x_var,
            n_boot=int(args.n_boot),
            subsample_frac=float(mcfg.get("subsample_frac", 0.8)),
            top_k=int(args.top_k),
            mode=mode,
            k=int(mcfg.get("k", 10)),
            iters=int(mcfg.get("iters", 5)),
            alpha=float(mcfg.get("alpha", 0.5)),
            solver=solver,
            lambda_=float(mcfg.get("lambda_", 1e-6)),
            normalize=str(mcfg.get("normalize", "std")),
            rescale_median=float(mcfg.get("rescale_median", 1.0)),
            seed=seed + 100 * i,
            cache_dir=str(cfg.get("cache_dir", "data/cache")),
        )

        tip = np.asarray(tip_out["tip"], dtype=float)
        tip_cache[name] = tip
        idx = _top_idx(tip, int(args.top_k))
        genes = genes_var[idx].astype(str).tolist()
        top_genes_cache[name] = genes

        for rank, g in enumerate(genes, start=1):
            top_gene_rows.append({"method": name, "rank": int(rank), "gene": g, "tip": float(tip[idx[rank - 1]])})

        detail, summ = _coherence_for_method(genes, genes_var.astype(str).tolist(), pathways)
        for r in detail:
            r2 = {"method": name, **r}
            pathway_rows.append(r2)
        pathway_summary_rows.append({"method": name, **summ})

    pd.DataFrame(top_gene_rows).to_csv("results/tables/luad_top_features_by_method.csv", index=False)
    pd.DataFrame(pathway_rows).to_csv("results/tables/luad_pathway_overlap_detail.csv", index=False)
    path_summary = pd.DataFrame(pathway_summary_rows)
    path_summary.to_csv("results/tables/luad_pathway_coherence.csv", index=False)

    surv = _parse_survival(dcfg.get("labels_path", "data/tcga_luad/labels.tsv.gz"))
    sample_df = pd.DataFrame(
        {
            "sample_id": [normalize_sample_id(s) for s in sample_ids],
            "sample_row": np.arange(len(sample_ids), dtype=int),
        }
    )
    merged = sample_df.merge(surv, on="sample_id", how="inner")

    sidx = merged["sample_row"].to_numpy(dtype=int)
    x_surv = x_var[sidx]
    surv_time = merged["time_days"].to_numpy(dtype=float)
    surv_event = merged["event"].to_numpy(dtype=int)

    surv_rows = []
    km_lines = []

    for name, genes in top_genes_cache.items():
        if len(genes) == 0:
            continue
        gene_to_idx = {g: i for i, g in enumerate(genes_var.astype(str).tolist())}
        idx = np.array([gene_to_idx[g] for g in genes if g in gene_to_idx], dtype=int)
        if idx.size == 0:
            continue

        z, _ = standardize(x_surv[:, idx], x_surv[:, idx])
        risk = np.mean(z, axis=1)
        thr = float(np.median(risk))
        group = (risk >= thr).astype(int)

        z_stat, p_val = _logrank(surv_time, surv_event, group)
        cidx = _c_index(surv_time, surv_event, risk)

        med_hi = _median_survival(surv_time[group == 1], surv_event[group == 1])
        med_lo = _median_survival(surv_time[group == 0], surv_event[group == 0])

        surv_rows.append(
            {
                "method": name,
                "n": int(len(group)),
                "n_high": int(np.sum(group == 1)),
                "n_low": int(np.sum(group == 0)),
                "logrank_z": float(z_stat),
                "logrank_p": float(p_val),
                "neglog10_p": float(-np.log10(max(1e-300, p_val))),
                "c_index": float(cidx),
                "median_survival_high": float(med_hi) if np.isfinite(med_hi) else np.nan,
                "median_survival_low": float(med_lo) if np.isfinite(med_lo) else np.nan,
            }
        )

        t0, s0 = _km_curve(surv_time[group == 0], surv_event[group == 0])
        t1, s1 = _km_curve(surv_time[group == 1], surv_event[group == 1])
        km_lines.append((name, t0, s0, t1, s1, p_val))

    surv_df = pd.DataFrame(surv_rows)
    if not surv_df.empty:
        surv_df = surv_df.sort_values("logrank_p", ascending=True)
    surv_df.to_csv("results/tables/luad_survival_stratification.csv", index=False)

    # Figure: pathway coherence + survival significance
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    p = path_summary.sort_values("best_neglog10_p", ascending=False)
    axes[0].bar(np.arange(len(p)), p["best_neglog10_p"])
    axes[0].set_xticks(np.arange(len(p)))
    axes[0].set_xticklabels(p["method"], rotation=30, ha="right")
    axes[0].set_title("LUAD Pathway Coherence")
    axes[0].set_ylabel("Best -log10 hypergeom p")

    if not surv_df.empty:
        s = surv_df.sort_values("neglog10_p", ascending=False)
        axes[1].bar(np.arange(len(s)), s["neglog10_p"])
        axes[1].set_xticks(np.arange(len(s)))
        axes[1].set_xticklabels(s["method"], rotation=30, ha="right")
        axes[1].set_title("LUAD Survival Stratification")
        axes[1].set_ylabel("-log10 logrank p")
    else:
        axes[1].axis("off")
    fig.tight_layout()
    fig.savefig("results/figures/luad_biological_diagnostics.png", dpi=220)
    plt.close(fig)

    # KM for best two methods by logrank p (if available)
    if km_lines:
        km_sorted = sorted(km_lines, key=lambda z: z[-1])[:2]
        fig, axes = plt.subplots(1, len(km_sorted), figsize=(5.4 * len(km_sorted), 4.0), sharey=True)
        if len(km_sorted) == 1:
            axes = [axes]
        for ax, (name, t0, s0, t1, s1, pval) in zip(axes, km_sorted):
            ax.step(t0, s0, where="post", label="Low-risk")
            ax.step(t1, s1, where="post", label="High-risk")
            ax.set_title(f"{name} (p={pval:.3g})")
            ax.set_xlabel("Days")
            ax.set_ylabel("KM survival")
            ax.set_ylim(0, 1.02)
            ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig("results/figures/luad_survival_km_by_method.png", dpi=220)
        plt.close(fig)

    print("Saved LUAD biological diagnostics tables/figures")


if __name__ == "__main__":
    main()
