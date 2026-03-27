#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_run_luad_multiseed
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G

set -euo pipefail
mkdir -p results/logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
python3.11 experiments/run_luad_multiseed_stats.py --config configs/tcga_luad.yaml --seeds 42,52,62,72,82,92 --n-splits 3 --out-prefix luad_multiseed
