#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _load(path: Path, dataset_fallback: str | None) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "dataset" not in df.columns and dataset_fallback is not None:
        df["dataset"] = dataset_fallback
    df["source_file"] = path.name
    return df


def _merge_family(results_dir: Path, family: str, candidates: list[tuple[str, str | None]], out_name: str) -> None:
    frames: list[pd.DataFrame] = []
    for rel_name, dataset_fallback in candidates:
        p = results_dir / rel_name
        if p.exists() and p.stat().st_size > 0:
            frames.append(_load(p, dataset_fallback))

    if not frames:
        print(f"[skip] no files found for {family}")
        return

    merged = pd.concat(frames, ignore_index=True)
    out_path = results_dir / out_name
    merged.to_csv(out_path, index=False)
    print(f"[ok] wrote {out_path} ({len(merged)} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description="Merge Hellbender CC-SPR per-dataset outputs")
    ap.add_argument("--results-dir", default="results")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        raise FileNotFoundError(f"results directory not found: {results_dir}")

    _merge_family(
        results_dir,
        family="metrics",
        candidates=[
            ("metrics.csv", "tcga_luad_existing"),
            ("metrics_cll_venetoclax_gse161711.csv", None),
            ("metrics_cll_rs_scrna_gse165087.csv", None),
        ],
        out_name="metrics_all_datasets.csv",
    )

    _merge_family(
        results_dir,
        family="main_results",
        candidates=[
            ("main_results.csv", "tcga_luad_existing"),
            ("main_results_cll_venetoclax_gse161711.csv", None),
            ("main_results_cll_rs_scrna_gse165087.csv", None),
        ],
        out_name="main_results_all_datasets.csv",
    )

    _merge_family(
        results_dir,
        family="ablation_summary",
        candidates=[
            ("ablation_summary.csv", "tcga_luad_existing"),
            ("ablation_summary_cll_venetoclax_gse161711.csv", None),
            ("ablation_summary_cll_rs_scrna_gse165087.csv", None),
        ],
        out_name="ablation_summary_all_datasets.csv",
    )


if __name__ == "__main__":
    main()
