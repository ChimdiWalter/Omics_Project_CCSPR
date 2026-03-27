#!/usr/bin/env bash
set -euo pipefail
cd /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src

TAG="cll_rs_scrna_gse165087_ablation_safe"
if [[ -f "results/metrics_${TAG}.csv" && -f "results/main_results_${TAG}.csv" && -f "results/ablation_summary_${TAG}.csv" ]]; then
  echo "[$(date --iso-8601=seconds)] Final scRNA ablation-safe outputs already exist; skipping"
  exit 0
fi

python3.11 -u experiments/run_cll_rs_scrna_ablation_safe.py --config configs/cll_rs_scrna_ablation_safe.yaml
