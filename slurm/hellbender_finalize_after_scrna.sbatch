#!/usr/bin/env bash
#SBATCH -J ccspr_finalize_post
#SBATCH -p requeue
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr/results/logs/ccspr_finalize_%j.out
#SBATCH -e /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr/results/logs/ccspr_finalize_%j.err

set -euo pipefail

cd /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src

bash scripts/finalize_hellbender_after_jobs.sh 12833129,12833143 | tee results/logs/finalize_hellbender_auto_%j.log
