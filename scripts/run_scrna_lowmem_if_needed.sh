#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-8}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-8}

REQ=(
  "results/metrics_cll_rs_scrna_gse165087.csv"
  "results/main_results_cll_rs_scrna_gse165087.csv"
  "results/ablation_summary_cll_rs_scrna_gse165087.csv"
)

ready=1
for f in "${REQ[@]}"; do
  if [[ ! -s "$f" ]]; then
    ready=0
    break
  fi
done

if [[ "$ready" -eq 1 ]]; then
  echo "scRNA outputs already present; low-memory fallback skipped"
  exit 0
fi

echo "scRNA outputs missing; running low-memory fallback config"
python3.11 experiments/run_cll_rs_scrna.py --config configs/cll_rs_scrna_lowmem.yaml

missing=0
for f in "${REQ[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "MISSING_AFTER_LOWMEM: $f"
    missing=1
  fi
done
if [[ "$missing" -ne 0 ]]; then
  echo "ERROR: required scRNA outputs still missing after low-memory run"
  exit 1
fi

echo "Merging consolidated tables and refreshing full manuscript"
python3.11 scripts/merge_hellbender_results.py --results-dir results
bash scripts/refresh_full_manuscript_after_runs.sh

echo "low-memory fallback completed"
