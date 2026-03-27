# CC-SPR Reproducibility Makefile
# ================================
# Usage: make <target>
# Run `make help` for a list of targets.

SHELL := /bin/bash
PYTHON := python3
VENV := .venv
PIP := $(VENV)/bin/pip
PYTHON_VENV := $(VENV)/bin/python3
PDFLATEX := pdflatex
MANUSCRIPT_DIR := manuscript_claude

.PHONY: help setup download-data run-pilot run-full run-arabidopsis run-all manuscript clean

help:  ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

# ── Environment ──────────────────────────────────────────────────────

setup: $(VENV)/bin/activate  ## Create virtualenv and install all dependencies
$(VENV)/bin/activate:
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip setuptools wheel
	$(PIP) install -e .
	@echo "Activate with: source $(VENV)/bin/activate"

# ── Data Download ────────────────────────────────────────────────────

download-data: setup  ## Download all datasets (LUAD, BRCA, CLL, Arabidopsis)
	$(PYTHON_VENV) scripts/download_tcga_luad.py
	$(PYTHON_VENV) scripts/download_tcga_brca.py
	$(PYTHON_VENV) scripts/download_cll_venetoclax.py --root data/cll_venetoclax
	$(PYTHON_VENV) scripts/download_cll_rs_scrna.py \
		--root data/cll_rs_scrna \
		--max-cells-raw 40000 \
		--max-cells-processed 12000 \
		--seed 42
	$(PYTHON_VENV) scripts/download_arabidopsis_root.py

download-luad: setup  ## Download TCGA-LUAD only
	$(PYTHON_VENV) scripts/download_tcga_luad.py

download-arabidopsis: setup  ## Download Arabidopsis root atlas only
	$(PYTHON_VENV) scripts/download_arabidopsis_root.py

# ── Experiments ──────────────────────────────────────────────────────

run-pilot: setup  ## LUAD pilot run (n=60, ~15 min)
	$(PYTHON_VENV) experiments/run_luad_bio_diagnostics.py \
		--config configs/geometry_paper_ultrafast.yaml \
		--max-samples 60 \
		--n-boot 1 \
		--top-k 20

run-full: setup  ## Full-scale LUAD with all 5 geometry backends (~4-8 h)
	$(PYTHON_VENV) experiments/run_luad_full_scale.py
	$(PYTHON_VENV) experiments/run_geometry_comparison.py \
		--config configs/geometry_paper_local.yaml
	$(PYTHON_VENV) experiments/run_luad_multiseed_geometry.py \
		--base-config configs/geometry_paper_ultrafast.yaml \
		--seeds 11 23 42 77 101 \
		--n-splits 2 \
		--max-samples 36 \
		--var-genes 120
	$(PYTHON_VENV) experiments/run_geometry_runtime_profile.py \
		--config configs/geometry_paper_ultrafast.yaml \
		--repeats 2

run-arabidopsis: setup  ## Arabidopsis root atlas experiment
	$(PYTHON_VENV) experiments/run_arabidopsis.py

run-cll: setup  ## CLL public substitute experiments
	$(PYTHON_VENV) experiments/run_cll_venetoclax.py --config configs/cll_venetoclax.yaml
	$(PYTHON_VENV) experiments/run_cll_rs_scrna.py --config configs/cll_rs_scrna.yaml
	$(PYTHON_VENV) experiments/aggregate_public_substitutes.py

run-all: download-data run-pilot run-full run-arabidopsis run-cll  ## Run everything end-to-end

# ── Manuscript ───────────────────────────────────────────────────────

manuscript: setup  ## Compile LaTeX manuscript to PDF
	$(PYTHON_VENV) experiments/generate_fast_diagnostics.py
	cd $(MANUSCRIPT_DIR) && $(PDFLATEX) -interaction=nonstopmode -halt-on-error cc_spr_revised_full.tex
	cd $(MANUSCRIPT_DIR) && $(PDFLATEX) -interaction=nonstopmode -halt-on-error cc_spr_revised_full.tex

# ── Cleanup ──────────────────────────────────────────────────────────

clean:  ## Remove generated artifacts (keeps data and results CSVs)
	rm -rf $(VENV)
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name '*.pyc' -delete 2>/dev/null || true
	find . -name '*.egg-info' -type d -exec rm -rf {} + 2>/dev/null || true
	cd $(MANUSCRIPT_DIR) && rm -f *.aux *.log *.out *.toc *.bbl *.blg 2>/dev/null || true
	@echo "Clean complete. Data and result CSVs preserved."
