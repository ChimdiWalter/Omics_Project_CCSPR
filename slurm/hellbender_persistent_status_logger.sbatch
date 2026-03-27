#!/usr/bin/env bash
#SBATCH -J ccspr_status_watch
#SBATCH -p requeue
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 2-00:00:00
#SBATCH -o /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr/results/logs/ccspr_status_watch_%j.out
#SBATCH -e /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr/results/logs/ccspr_status_watch_%j.err

set -euo pipefail

cd /cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate

bash scripts/persistent_status_logger.sh \
  "12833143,12833568" \
  "results/logs/persistent_status_watch_${SLURM_JOB_ID}.log" \
  "60"
