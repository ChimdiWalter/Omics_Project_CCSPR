#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_geometry_family
#SBATCH --output=results/logs/geometry_family_%j.out
#SBATCH --error=results/logs/geometry_family_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

###############################################################################
# CC-SPR Geometry-Family Empirical Study
#
# Runs all 5 geometry backends (euclidean, ricci, diffusion, phate_like, dtm)
# across LUAD (60 samples), GSE161711 (96 samples), GSE165087 (320 bounded).
#
# Checkpoints saved after each phase. Safe to re-run after interruption.
#
# Resource budget: 48 GB RAM, 8 CPUs, 12h wall clock
# Expected usage: ~12-20 GB peak, ~2-6h total
###############################################################################

set -euo pipefail

mkdir -p results/logs results/geometry_family/tables

echo "============================================================"
echo "[$(date -Iseconds)] CC-SPR Geometry-Family starting"
echo "[$(date -Iseconds)] Job ID: ${SLURM_JOB_ID:-local}"
echo "[$(date -Iseconds)] Node: $(hostname)"
echo "[$(date -Iseconds)] CPUs: ${SLURM_CPUS_PER_TASK:-8}"
echo "[$(date -Iseconds)] Memory: ${SLURM_MEM_PER_NODE:-48G}"
echo "============================================================"

# Activate environment
source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src
export MPLCONFIGDIR=/tmp/matplotlib
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export PYTHONUNBUFFERED=1

echo "[$(date -Iseconds)] Python: $(which python3.11 2>/dev/null || which python3)"
echo "[$(date -Iseconds)] Working directory: $(pwd)"

# Run geometry family experiments with full logging
python3.11 -u experiments/run_geometry_family.py \
    2>&1 | tee results/logs/geometry_family_experiments_${SLURM_JOB_ID:-local}.log

EXIT_CODE=$?
echo "[$(date -Iseconds)] Exit code: ${EXIT_CODE}"

# Print final outputs
echo ""
echo "Final outputs:"
for f in \
    results/geometry_family/checkpoint.json \
    results/geometry_family/geometry_family_results.json \
    results/geometry_family/tables/geometry_family_benchmark.csv \
    results/geometry_family/tables/geometry_family_diagnostics.csv \
    results/geometry_family/tables/geometry_family_ablation.csv \
    results/geometry_family/tables/geometry_family_timing.csv; do
    if [ -f "$f" ]; then
        echo "  [OK]   $f ($(stat -c%s "$f" 2>/dev/null || echo '?') bytes)"
    else
        echo "  [MISS] $f"
    fi
done

echo ""
echo "[$(date -Iseconds)] Done."
echo "============================================================"
