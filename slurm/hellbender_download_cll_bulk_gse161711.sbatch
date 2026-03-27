#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_dl_bulk161711
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

set -euo pipefail
mkdir -p results/logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
python3.11 scripts/download_cll_venetoclax.py --root data/cll_venetoclax
