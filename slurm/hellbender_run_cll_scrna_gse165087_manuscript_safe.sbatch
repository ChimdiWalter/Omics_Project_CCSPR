#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_run_scrna165087_ms
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
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-6}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-6}

python3.11 -u experiments/run_cll_rs_scrna_manuscript_safe.py --config configs/cll_rs_scrna_manuscript_safe.yaml 2>&1 | tee results/logs/ccspr_run_scrna165087_ms_${SLURM_JOB_ID}.log
