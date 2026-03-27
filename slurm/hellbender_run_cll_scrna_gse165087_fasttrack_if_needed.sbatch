#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_run_scrna165087_ft
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G

set -euo pipefail
mkdir -p results/logs
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src

bash scripts/run_scrna_fasttrack_if_needed.sh | tee results/logs/scrna_fasttrack_if_needed_${SLURM_JOB_ID}.log
