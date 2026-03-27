#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_run_scrna165087
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G

set -euo pipefail
mkdir -p results/logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
python3.11 experiments/run_cll_rs_scrna.py --config configs/cll_rs_scrna.yaml
