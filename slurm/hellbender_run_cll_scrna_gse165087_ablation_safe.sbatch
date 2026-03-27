#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_run_scrna165087_abs
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G

set -euo pipefail
mkdir -p results/logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib

bash scripts/run_scrna_ablation_safe_if_needed.sh 2>&1 | tee results/logs/scrna_ablation_safe_if_needed_${SLURM_JOB_ID}.log
