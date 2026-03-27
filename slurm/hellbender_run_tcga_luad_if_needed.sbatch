#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_run_luad
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

set -euo pipefail
mkdir -p results/logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib

if [[ -s results/main_results.csv && -s results/ablation_summary.csv && -s results/metrics.csv ]]; then
  echo "TCGA-LUAD results already present; skipping rerun to preserve completed work."
  exit 0
fi

python3.11 experiments/run_tcga_luad.py --config configs/tcga_luad.yaml
