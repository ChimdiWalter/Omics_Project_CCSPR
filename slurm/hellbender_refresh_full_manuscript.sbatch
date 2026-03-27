#!/usr/bin/env bash
#SBATCH -J ccspr_refresh_manuscript
#SBATCH -p requeue
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr/results/logs/ccspr_refresh_manuscript_%j.out
#SBATCH -e /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr/results/logs/ccspr_refresh_manuscript_%j.err

set -euo pipefail
cd /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate

bash scripts/refresh_full_manuscript_after_runs.sh | tee results/logs/refresh_full_manuscript_%j.log
