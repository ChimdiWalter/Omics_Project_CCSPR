#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

from ccspr.datasets.tcga_brca import download_tcga_brca, load_tcga_brca_multiomics
from ccspr.eval.classification import FeatureParams, evaluate_luad_protocol
from ccspr.utils.io import ensure_dir, read_yaml, save_json


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/tcga_brca.yaml")
    args = ap.parse_args()
    cfg = read_yaml(args.config)

    out = Path(cfg["output"]["results_dir"])
    ensure_dir(out)
    save_json({"timestamp": datetime.utcnow().isoformat() + "Z", "config": cfg}, out / "run.json")

    if cfg["dataset"].get("download", True):
        download_tcga_brca(cfg["dataset"].get("root", "data/tcga_brca"))

    ds = load_tcga_brca_multiomics(
        root=cfg["dataset"].get("root", "data/tcga_brca"),
        pca_expr=int(cfg["dataset"].get("pca_expr", 100)),
        pca_meth=int(cfg["dataset"].get("pca_meth", 100)),
        max_samples=cfg["dataset"].get("max_samples"),
        seed=int(cfg.get("seed", 42)),
    )

    harmonic_params = None
    ccspr_params = None
    if "harmonic_params" in cfg.get("protocol", {}):
        harmonic_params = FeatureParams(**cfg["protocol"]["harmonic_params"])
    if "ccspr_params" in cfg.get("protocol", {}):
        ccspr_params = FeatureParams(**cfg["protocol"]["ccspr_params"])

    metrics, _ = evaluate_luad_protocol(
        ds["X"],
        ds["y"],
        var_genes=min(int(cfg["protocol"].get("var_genes", 1000)), ds["X"].shape[1]),
        n_splits=int(cfg["protocol"].get("n_splits", 3)),
        test_size=float(cfg["protocol"].get("test_size", 0.25)),
        seed=int(cfg.get("seed", 42)),
        cache_dir=cfg["cache_dir"],
        harmonic_params=harmonic_params,
        ccspr_params=ccspr_params,
    )

    metrics.to_csv(out / "metrics.csv", index=False)
    pd.DataFrame([ds["meta"]]).to_csv(out / "dataset_meta.csv", index=False)
    print("Saved:", out / "metrics.csv")


if __name__ == "__main__":
    main()
