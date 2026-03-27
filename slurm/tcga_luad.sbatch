#!/usr/bin/env bash
#SBATCH --job-name=ccspr_luad
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail
mkdir -p logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
python3.11 experiments/run_tcga_luad.py --config configs/tcga_luad.yaml
