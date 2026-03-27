# Reproducing CC-SPR Results

Step-by-step instructions to reproduce all results reported in the CC-SPR manuscript.

## Prerequisites

- **OS:** Linux (tested on Ubuntu 20.04+)
- **Python:** 3.10 or later
- **RAM:** 32 GB minimum; 125 GB recommended for full-scale LUAD
- **CPU:** Multi-core recommended (experiments used Intel i9-9920X, 24 threads)
- **GPU:** Not required (CPU-only pipeline)
- **Disk:** ~10 GB for all datasets and results

## Step 0: Environment Setup

```bash
git clone https://github.com/ChimdiWalter/cc-spr.git
cd cc-spr

python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip setuptools wheel
pip install -e .
```

Set recommended environment variables:

```bash
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
export NUMBA_DISABLE_COVERAGE=1
export NUMBA_CACHE_DIR=/tmp/numba_cache
```

## Step 1: Download Datasets

### TCGA-LUAD (primary benchmark)

```bash
python3 scripts/download_tcga_luad.py
```

Downloads and preprocesses TCGA lung adenocarcinoma expression + clinical data into `data/`.

### Arabidopsis Root Atlas (cross-domain test)

```bash
python3 scripts/download_arabidopsis_root.py
```

Downloads GSE123013 scRNA-seq data for plant root developmental analysis.

### TCGA-BRCA

```bash
python3 scripts/download_tcga_brca.py
```

### CLL Public Substitutes (optional)

```bash
python3 scripts/download_cll_venetoclax.py --root data/cll_venetoclax
python3 scripts/download_cll_rs_scrna.py \
  --root data/cll_rs_scrna \
  --max-cells-raw 40000 \
  --max-cells-processed 12000 \
  --seed 42
```

## Step 2: LUAD Pilot Run (Quick Validation)

Estimated time: ~15 minutes.

```bash
python3 experiments/run_luad_bio_diagnostics.py \
  --config configs/geometry_paper_ultrafast.yaml \
  --max-samples 60 \
  --n-boot 1 \
  --top-k 20
```

This runs a small-scale version of the LUAD benchmark (n=60) to verify the pipeline works before committing to the full run.

## Step 3: Full-Scale LUAD Benchmark (Table 1 in Manuscript)

Estimated time: 4-8 hours depending on hardware.

### 3a. Full-scale evaluation (n=275, all backends)

```bash
python3 experiments/run_luad_full_scale.py
```

Expected results in `results/luad_full/tables/luad_full_benchmark.csv`:

| Backend | Weighted F1 |
|---|---|
| Ricci | 0.6339 |
| DTM | 0.5814 |
| Euclidean | 0.5655 |
| Standard baseline | 0.7501 |

### 3b. Multi-geometry comparison

```bash
python3 experiments/run_geometry_comparison.py \
  --config configs/geometry_paper_local.yaml
```

### 3c. Multi-seed paired statistics

```bash
python3 experiments/run_luad_multiseed_geometry.py \
  --base-config configs/geometry_paper_ultrafast.yaml \
  --seeds 11 23 42 77 101 \
  --n-splits 2 \
  --max-samples 36 \
  --var-genes 120
```

### 3d. Runtime/cost profiling

```bash
python3 experiments/run_geometry_runtime_profile.py \
  --config configs/geometry_paper_ultrafast.yaml \
  --repeats 2
```

## Step 4: Arabidopsis Root Experiment (Table 2 in Manuscript)

Estimated time: 1-3 hours.

```bash
python3 experiments/run_arabidopsis.py
```

Expected results:

| Backend | Weighted F1 |
|---|---|
| DTM | 0.6664 |
| Euclidean | 0.5933 |

Results are written to `results/arabidopsis/`.

## Step 5: CLL Public Substitute Experiments (Optional)

```bash
python3 experiments/run_cll_venetoclax.py --config configs/cll_venetoclax.yaml
python3 experiments/run_cll_rs_scrna.py --config configs/cll_rs_scrna.yaml
python3 experiments/aggregate_public_substitutes.py
```

## Step 6: Generate Paper Assets and Compile Manuscript

```bash
# Generate diagnostic figures and tables
python3 experiments/generate_fast_diagnostics.py

# Compile the manuscript (run twice for cross-references)
cd manuscript_claude
pdflatex -interaction=nonstopmode -halt-on-error cc_spr_revised_full.tex
pdflatex -interaction=nonstopmode -halt-on-error cc_spr_revised_full.tex
cd ..
```

The compiled PDF will be at `manuscript_claude/cc_spr_revised_full.pdf`.

## Output Summary

After full reproduction, key outputs are:

| File | Contents |
|---|---|
| `results/luad_full/tables/luad_full_benchmark.csv` | Main LUAD benchmark (Table 1) |
| `results/luad_full/tables/luad_full_top_features_by_method.csv` | Top selected features per backend |
| `results/arabidopsis/tables/` | Arabidopsis benchmark tables |
| `results/tables/luad_multiseed_paired_stats.csv` | Paired statistical tests |
| `results/tables/runtime_cost_geometry_overall.csv` | Runtime profiling |
| `results/figures/` | All generated figures |
| `manuscript_claude/cc_spr_revised_full.pdf` | Compiled manuscript |

## One-Command Reproduction

To run everything automatically:

```bash
make run-all
make manuscript
```

## Troubleshooting

- **Out of memory on full LUAD run:** Reduce to pilot mode with `--max-samples 60`, or use `configs/geometry_paper_ultrafast.yaml` which limits sample/gene counts.
- **GUDHI installation fails:** Ensure C++ build tools are available (`apt install build-essential cmake`).
- **GraphRicciCurvature slow:** The Ricci backend is the most expensive geometry. On machines with fewer than 8 cores, expect longer runtimes for this backend.
- **Missing LaTeX:** Install with `apt install texlive-latex-base texlive-latex-extra texlive-fonts-recommended` or equivalent.
