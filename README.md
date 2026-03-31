# CC-SPR: Cycle-Consistent Sparse Persistent Representatives

A topological data analysis (TDA) framework for geometry-sensitive feature selection in omics data. CC-SPR uses persistent homology with multiple curvature-aware distance backends to identify biologically coherent gene modules from high-dimensional expression matrices.

**Target journal:** Bioinformatics (Oxford) — manuscript in preparation

## Key Results

| Experiment | Backend | Weighted F1 |
|---|---|---|
| LUAD full-scale (n=275) | Ricci | 0.6339 |
| LUAD full-scale (n=275) | DTM | 0.5814 |
| LUAD full-scale (n=275) | Euclidean | 0.5655 |
| LUAD full-scale (n=275) | Standard baseline | 0.7501 |
| Arabidopsis root (n=500) | DTM | 0.6664 |
| Arabidopsis root (n=500) | Euclidean | 0.5933 |

Five geometry backends are supported: `euclidean`, `ricci`, `diffusion`, `phate_like`, `dtm`.

## Installation

Requires Python 3.10+.

```bash
# Clone the repository
git clone https://github.com/ChimdiWalter/cc-spr.git
cd cc-spr

# Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install in editable mode
pip install -e .
```

Or use the Makefile:

```bash
make setup
```

### Dependencies

Core dependencies are declared in `pyproject.toml` and include:
numpy, pandas, scipy, scikit-learn, matplotlib, networkx, gudhi, cvxpy, GraphRicciCurvature, PyYAML, anndata, scanpy, joblib, requests.

A frozen environment is available in `requirements-frozen.txt`.

## Quick Start

```bash
# 1. Download TCGA-LUAD data
python3 scripts/download_tcga_luad.py

# 2. Run the pilot experiment (n=60, fast)
python3 experiments/run_luad_bio_diagnostics.py \
  --config configs/geometry_paper_ultrafast.yaml \
  --max-samples 60 --n-boot 1 --top-k 20

# 3. Run the full-scale LUAD evaluation (n=275, all backends)
python3 experiments/run_luad_full_scale.py
```

## Full Reproduction

See [REPRODUCE.md](REPRODUCE.md) for detailed step-by-step reproduction instructions, or use the Makefile:

```bash
make run-all        # Run all experiments end-to-end
make manuscript     # Compile the PDF
```

### Individual Makefile Targets

| Target | Description |
|---|---|
| `make setup` | Create virtualenv and install dependencies |
| `make download-data` | Download all datasets (LUAD, BRCA, CLL, Arabidopsis) |
| `make run-pilot` | LUAD pilot run (n=60, ~15 min) |
| `make run-full` | Full-scale LUAD with all 5 backends (~4-8 h) |
| `make run-arabidopsis` | Arabidopsis root atlas experiment |
| `make run-all` | All experiments sequentially |
| `make manuscript` | Compile LaTeX manuscript to PDF |
| `make clean` | Remove generated artifacts |

## Project Structure

```
cc-spr/
├── src/ccspr/                 # Core library
│   ├── datasets/              # Data loaders (TCGA-LUAD, BRCA, CLL, Arabidopsis)
│   ├── preprocess/            # Normalization, variable gene selection
│   ├── geometry/              # Distance backends (Euclidean, Ricci, DTM, ...)
│   ├── topology/              # Persistent homology, H1 bar extraction, cycles
│   ├── solver/                # L2/elastic cycle representative solvers
│   ├── stability/             # TIP bootstrap
│   ├── eval/                  # Weighted-F1 evaluation, ablations
│   ├── plots/                 # Figure generation
│   └── utils/                 # I/O, caching
├── configs/                   # YAML experiment configurations
├── experiments/               # Runnable experiment scripts
│   ├── run_luad_full_scale.py
│   ├── run_arabidopsis.py
│   ├── run_geometry_comparison.py
│   ├── run_luad_multiseed_geometry.py
│   └── ...
├── scripts/                   # Data download and utility scripts
│   ├── download_tcga_luad.py
│   ├── download_arabidopsis_root.py
│   └── ...
├── manuscript_new/            # LaTeX manuscript sources
│   └── cc_spr_revised_full.tex  # Main manuscript (37 pp)
├── results/                   # Generated results (tables, figures)
├── pyproject.toml             # Package metadata and dependencies
├── Makefile                   # Reproducibility automation
├── REPRODUCE.md               # Step-by-step reproduction guide
└── LICENSE                    # MIT License
```

## Hardware

All experiments were run on a single workstation:
- **CPU:** Intel Core i9-9920X (24 threads)
- **RAM:** 125 GB
- **GPU:** Not required (CPU-only pipeline)

## Citation

If you use CC-SPR in your research, please cite:

```bibtex
@article{ndubuisi2026ccspr,
  title={CC-SPR: Cycle-Consistent Sparse Persistent Representatives for
         Geometry-Sensitive Topological Feature Selection in Omics Data},
  author={Ndubuisi, Chimdi Walter},
  year={2026}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
