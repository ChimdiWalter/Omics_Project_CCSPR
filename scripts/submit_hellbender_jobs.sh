#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"
mkdir -p results/logs

j_dl_bulk=$(sbatch --parsable slurm/hellbender_download_cll_bulk_gse161711.sbatch)
j_dl_scrna=$(sbatch --parsable slurm/hellbender_download_cll_scrna_gse165087.sbatch)
j_run_bulk=$(sbatch --parsable --dependency=afterok:${j_dl_bulk} slurm/hellbender_run_cll_bulk_gse161711.sbatch)
j_run_scrna=$(sbatch --parsable --dependency=afterok:${j_dl_scrna} slurm/hellbender_run_cll_scrna_gse165087.sbatch)
j_luad=$(sbatch --parsable slurm/hellbender_run_tcga_luad_if_needed.sbatch)

ts=$(date +%Y%m%d_%H%M%S)
out="results/logs/submitted_jobs_${ts}.txt"
{
  echo "timestamp=${ts}"
  echo "download_bulk=${j_dl_bulk}"
  echo "download_scrna=${j_dl_scrna}"
  echo "run_bulk=${j_run_bulk}"
  echo "run_scrna=${j_run_scrna}"
  echo "tcga_luad_if_needed=${j_luad}"
} | tee "$out"

echo "JOB_IDS ${j_dl_bulk} ${j_dl_scrna} ${j_run_bulk} ${j_run_scrna} ${j_luad}"
