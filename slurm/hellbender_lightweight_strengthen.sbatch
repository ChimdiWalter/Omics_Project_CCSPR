#!/usr/bin/env bash
#SBATCH --chdir=/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr
#SBATCH --job-name=ccspr_lightweight_strengthen
#SBATCH --output=results/logs/%x_%j.out
#SBATCH --error=results/logs/%x_%j.err
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

###############################################################################
# CC-SPR Lightweight Manuscript Strengthening
#
# Runs all 4 priority tasks in order:
#   1. LUAD multiseed paired analysis
#   2. GSE161711 bounded ablation sweep (small grid)
#   3. Diagnostics from existing outputs
#   4. Optional tiny scRNA k-sensitivity
#
# Then auto-updates manuscript and compiles PDFs.
#
# All progress is logged with timestamps. Each task skips if outputs exist,
# so this script is safe to re-run after interruption.
#
# Resource budget: 48 GB RAM, 8 CPUs, 18h wall clock
# Expected usage:  ~12-16 GB peak, ~6-10h total
###############################################################################

set -euo pipefail

# --- Setup ---
mkdir -p results/logs results/tables results/figures results_claude

echo "============================================================"
echo "[$(date -Iseconds)] CC-SPR Lightweight Strengthen starting"
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

# Use unbuffered Python for continuous logging
export PYTHONUNBUFFERED=1

# Log environment
echo "[$(date -Iseconds)] Python: $(which python3.11 2>/dev/null || which python3)"
echo "[$(date -Iseconds)] Working directory: $(pwd)"

# --- Phase 1: Run all lightweight experiments ---
echo ""
echo "[$(date -Iseconds)] === PHASE 1: Lightweight experiments ==="
echo ""

python3 -u experiments/run_lightweight_strengthen.py \
    --luad-config configs/tcga_luad.yaml \
    --gse161711-config configs/cll_venetoclax_bounded_ablation.yaml \
    --scrna-config configs/cll_rs_scrna_tiny_sensitivity.yaml \
    2>&1 | tee results/logs/lightweight_strengthen_experiments_${SLURM_JOB_ID:-local}.log

PHASE1_EXIT=$?
echo "[$(date -Iseconds)] Phase 1 exit code: ${PHASE1_EXIT}"

# --- Phase 2: Refresh manuscript and compile PDFs ---
echo ""
echo "[$(date -Iseconds)] === PHASE 2: Manuscript refresh and PDF compilation ==="
echo ""

python3 -u scripts/refresh_manuscript_after_lightweight.py \
    2>&1 | tee results/logs/lightweight_strengthen_manuscript_${SLURM_JOB_ID:-local}.log

PHASE2_EXIT=$?
echo "[$(date -Iseconds)] Phase 2 exit code: ${PHASE2_EXIT}"

# --- Summary ---
echo ""
echo "============================================================"
echo "[$(date -Iseconds)] CC-SPR Lightweight Strengthen complete"
echo "[$(date -Iseconds)] Phase 1 (experiments): exit ${PHASE1_EXIT}"
echo "[$(date -Iseconds)] Phase 2 (manuscript):  exit ${PHASE2_EXIT}"
echo ""

# Print final outputs
echo "Final outputs:"
for f in \
    results/paired_stats_luad_multiseed_v2.csv \
    results/ablation_summary_cll_venetoclax_gse161711_bounded.csv \
    results/tables/lightweight_diagnostics.csv \
    results/ablation_summary_cll_rs_scrna_gse165087_tiny_sensitivity.csv \
    results/checkpoint_lightweight_strengthen.json \
    results_claude/lightweight_strengthen_report.md \
    manuscript_claude/cc_spr_full_paper_claude.pdf \
    manuscript_claude/cc_spr_full_supplement_claude.pdf; do
    if [ -f "$f" ]; then
        echo "  [OK]   $f ($(stat -c%s "$f" 2>/dev/null || echo '?') bytes)"
    else
        echo "  [MISS] $f"
    fi
done

echo ""
echo "[$(date -Iseconds)] Done."
echo "============================================================"
