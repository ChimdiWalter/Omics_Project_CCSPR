#!/usr/bin/env bash
#SBATCH --job-name=ccspr_cll_rs
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

set -euo pipefail
mkdir -p logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
python3.11 experiments/run_cll_rs_scrna.py --config configs/cll_rs_scrna.yaml
