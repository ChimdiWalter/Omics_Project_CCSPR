# CC-SPR on Public CLL Substitutes

## Abstract

This report presents a reproducible implementation of CC-SPR, a curvature-corrected sparse extension of harmonic persistent homology for omics data analysis. The originally requested CLL venetoclax bulk matrix and CLL/RS single-cell cohort were unavailable due access constraints, so we used public GEO substitutes chosen to preserve the CLL and longitudinal-treatment context: GSE161711 (bulk RNA-seq) and GSE165087 (scRNA-seq). We implemented end-to-end data ingestion, preprocessing, baselines, CC-SPR, full ablation grids, and figure/table generation in a local-first pipeline.

## 1. Problem Framing

We keep the original scientific framing: harmonic PH-style omics probing extended with:

1. Curvature-aware geometry (Ricci-corrected distances),
2. Sparse cycle representatives (elastic cycle solver),
3. TIP stability scoring for feature ranking.

The goal is not to replace standard classifiers, but to evaluate whether topological-sparse features are stable and informative in CLL-like datasets.

## 2. Datasets

### 2.1 Public Substitutes

- Bulk substitute: **GSE161711** (`Interdependencies of Intratumoral Genetic Heterogeneity and the Immune Microenvironment in CLL`), converted to:
  - `data/cll_venetoclax/processed_matrix.tsv`
- scRNA substitute: **GSE165087** (`Longitudinal CLL single-cell RNA-seq`), converted to:
  - `data/cll_rs_scrna/data.h5ad`

These substitutes were selected because they preserve CLL disease context and bulk/single-cell modality alignment while being public.

### 2.2 Data Construction in This Repo

- `scripts/download_cll_venetoclax.py` downloads GEO supplementary and series metadata, then builds a sample-by-feature matrix with `sample_id` and `label`.
- `scripts/download_cll_rs_scrna.py` downloads GEO RAW 10X archives, merges all samples, attaches metadata (including time-based labels), applies Scanpy QC/normalization/log1p/HVG/PCA, and writes `data.h5ad`.

## 3. Methods

### 3.1 Core Pipeline

For each train split:

1. Select top variable genes/features.
2. Build distance matrix (`euclidean` or `ricci`).
3. Compute H1 persistence bars.
4. Build cycle problem and solve representative (`l2` or `elastic`).
5. Project representative to feature scores.
6. Repeat by bootstrap and compute TIP top-K stability vectors.
7. Train/evaluate logistic regression on selected features.

### 3.2 Models Compared

- `standard`: direct feature baseline.
- `harmonic`: euclidean + L2 cycle representative.
- `eu_sparse`: euclidean + elastic cycle representative.
- `ccspr`: ricci + elastic cycle representative.

### 3.3 Ablation Grid (full grid retained)

- `distance_mode`: {euclidean, ricci}
- `iters`: {0, 2, 5, 10}
- `k`: {5, 10, 15}
- `lambda`: {1e-8, 1e-7, 1e-6, 1e-5, 1e-4}
- `top_k`: {10, 20, 50}
- `normalize`: {none, std}

## 4. Experimental Setup

- Deterministic seed: 42.
- Caching enabled for distance and persistence (`data/cache`).
- Local-scale scRNA execution used for feasibility (full grid preserved, reduced cell count / bootstrap depth for local runtime).

## 5. Results

### 5.1 Main Benchmark (Weighted F1)

| Dataset | Model | Mean F1 | 95% CI |
|---|---|---:|---:|
| GSE161711 bulk | standard | 0.8961 | 0.0404 |
| GSE161711 bulk | harmonic | 0.4804 | 0.1107 |
| GSE161711 bulk | eu_sparse | 0.5885 | 0.0931 |
| GSE161711 bulk | ccspr | 0.4818 | 0.0028 |
| GSE165087 scRNA | standard | 0.7819 | 0.0994 |
| GSE165087 scRNA | harmonic | 0.5214 | 0.1554 |
| GSE165087 scRNA | eu_sparse | 0.4704 | 0.0466 |
| GSE165087 scRNA | ccspr | 0.4734 | 0.1438 |

### 5.2 Ablation Highlights

- Bulk (GSE161711): best means in this local run favored `k=15`, `top_k=50`, `normalize=std`, and `lambda=1e-4`.
- scRNA (GSE165087): best means favored `iters=10`, `k=5`, `top_k=50`, and `lambda=1e-7`.

### 5.3 Figures Produced

- TIP Euclidean vs Ricci:
  - `results/figures/cll_venetoclax_tip_eu_vs_ricci.png`
  - `results/figures/cll_rs_scrna_tip_eu_vs_ricci.png`
- H1 lifetime/prominence Euclidean vs Ricci:
  - `results/figures/cll_venetoclax_h1_lifetime_eu_vs_ricci.png`
  - `results/figures/cll_venetoclax_h1_prominence_eu_vs_ricci.png`
  - `results/figures/cll_rs_scrna_h1_lifetime_eu_vs_ricci.png`
  - `results/figures/cll_rs_scrna_h1_prominence_eu_vs_ricci.png`
- Main summary / ablation over substitutes:
  - `results/figures/main_benchmark_public_substitutes.png`
  - `results/figures/ablation_panels_public_substitutes.png`

## 6. Limitations

1. Public substitutes are biologically aligned but not identical to the originally requested proprietary datasets.
2. Current reported scRNA run is local-scale (subsampled cells) to keep full ablation feasible on local compute.
3. Split count and bootstrap depth are intentionally conservative for runtime; larger values are recommended for publication-grade confidence intervals.

## 7. Future Work

1. Re-run full-scale scRNA with increased cell budget and bootstrap count.
2. Add patient-level / timepoint-level aggregation and survival-aware endpoints.
3. Add explicit R/maTilDA parity experiments in the same benchmark table.
4. Expand uncertainty estimates with nested cross-validation and repeated seeds.

## 8. Reproducibility Commands

```bash
source ~/.venvs/lesegenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
export NUMBA_DISABLE_COVERAGE=1
export NUMBA_CACHE_DIR=/tmp/numba_cache

python3 scripts/download_cll_venetoclax.py --root data/cll_venetoclax
python3 scripts/download_cll_rs_scrna.py --root data/cll_rs_scrna --max-cells-raw 40000 --max-cells-processed 12000 --seed 42

python3 experiments/run_cll_venetoclax.py --config configs/cll_venetoclax.yaml
python3 experiments/run_cll_rs_scrna.py --config configs/cll_rs_scrna.yaml
python3 experiments/aggregate_public_substitutes.py
```
