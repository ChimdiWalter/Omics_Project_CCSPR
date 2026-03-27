#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src

mkdir -p results/logs

scRNA_ready=0
if [[ -s results/metrics_cll_rs_scrna_gse165087.csv && -s results/main_results_cll_rs_scrna_gse165087.csv && -s results/ablation_summary_cll_rs_scrna_gse165087.csv ]]; then
  scRNA_ready=1
fi

luad_ms_ready=0
if [[ -s results/metrics_luad_multiseed.csv && -s results/paired_stats_luad_multiseed.csv ]]; then
  luad_ms_ready=1
fi

# Merge consolidated tables when scRNA outputs are available.
if [[ "$scRNA_ready" -eq 1 ]]; then
  python3.11 scripts/merge_hellbender_results.py --results-dir results
fi

# Recompile full paper/supplement.
(
  cd manuscript
  pdflatex -interaction=nonstopmode -halt-on-error cc_spr_full_paper.tex >/tmp/ccspr_full_paper_refresh_1.log 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error cc_spr_full_paper.tex >/tmp/ccspr_full_paper_refresh_2.log 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error cc_spr_full_supplement.tex >/tmp/ccspr_full_supp_refresh_1.log 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error cc_spr_full_supplement.tex >/tmp/ccspr_full_supp_refresh_2.log 2>&1
)

cat > results/manuscript_upgrade_report.md <<'REPORT'
# Manuscript Upgrade Report

Date: __DATE__

## Reused Existing Outputs
- LUAD baseline + ablations (`results/metrics.csv`, `results/main_results.csv`, `results/ablation_summary.csv`)
- Bulk substitute GSE161711 outputs (`results/*cll_venetoclax_gse161711*`)
- Existing figures under `results/figures/`
- Existing substitute manuscript as starting source

## New Work Added
- Full-paper assets:
  - `manuscript/cc_spr_full_paper.tex`
  - `manuscript/cc_spr_full_paper.pdf`
- Full-supplement assets:
  - `manuscript/cc_spr_full_supplement.tex`
  - `manuscript/cc_spr_full_supplement.pdf`
- Gap report:
  - `results/manuscript_gap_report.md`
- LUAD multiseed paired-statistics pipeline:
  - `experiments/run_luad_multiseed_stats.py`
  - `slurm/hellbender_run_luad_multiseed_stats.sbatch`

## Status of Pending Data-Dependent Inclusions
- scRNA outputs ready: __SCRNA_READY__
- LUAD multiseed outputs ready: __LUAD_READY__

## Tables/Figures Added in This Upgrade Pass
- New manuscript PDFs were compiled from full-paper/full-supplement sources.
- If scRNA outputs were ready, consolidated tables were refreshed:
  - `results/metrics_all_datasets.csv`
  - `results/main_results_all_datasets.csv`
  - `results/ablation_summary_all_datasets.csv`

## Caveats Remaining
- Backend-family claims beyond Euclidean/Ricci remain manuscript-scoped unless corresponding code and experiments are added.
- If `scRNA outputs ready=0`, final scRNA quantitative integration remains pending completion of job 12833143 and downstream finalizer.
- If `LUAD multiseed outputs ready=0`, paired LUAD multiseed statistics remain pending completion of job 12839604.
REPORT

sed -i "s/__DATE__/$(date --iso-8601=seconds)/" results/manuscript_upgrade_report.md
sed -i "s/__SCRNA_READY__/${scRNA_ready}/" results/manuscript_upgrade_report.md
sed -i "s/__LUAD_READY__/${luad_ms_ready}/" results/manuscript_upgrade_report.md

echo "refresh_complete $(date --iso-8601=seconds) scRNA_ready=${scRNA_ready} luad_ms_ready=${luad_ms_ready}"
