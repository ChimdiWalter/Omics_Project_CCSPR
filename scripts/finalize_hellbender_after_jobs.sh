#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

source /cluster/VAST/kazict-lab/e/lesion_phes/lesenv/bin/activate
export PYTHONPATH=src

JOBS="${1:-12833129,12833143}"

EXPECTED=(
  "results/metrics_cll_venetoclax_gse161711.csv"
  "results/main_results_cll_venetoclax_gse161711.csv"
  "results/ablation_summary_cll_venetoclax_gse161711.csv"
  "results/metrics_cll_rs_scrna_gse165087.csv"
  "results/main_results_cll_rs_scrna_gse165087.csv"
  "results/ablation_summary_cll_rs_scrna_gse165087.csv"
)

echo "Waiting for jobs: $JOBS"
while true; do
  if ! queue_rows="$(squeue -h -j "$JOBS" 2>&1)"; then
    echo "WARN: squeue query failed; retrying in 30s"
    echo "$queue_rows"
    sleep 30
    continue
  fi
  if [[ -z "${queue_rows//[[:space:]]/}" ]]; then
    break
  fi
  date
  squeue -j "$JOBS" -o '%.18i %.9P %.25j %.8T %.10M %.6D %R' || true
  sleep 30
done

echo "Jobs no longer in queue. Waiting for terminal accounting states..."
while true; do
  if ! acct_rows="$(sacct -j "$JOBS" --format=JobIDRaw,State,ExitCode,Elapsed -n -P 2>&1)"; then
    echo "WARN: sacct query failed; retrying in 30s"
    echo "$acct_rows"
    sleep 30
    continue
  fi
  echo "$acct_rows"
  if [[ -z "${acct_rows//[[:space:]]/}" ]]; then
    echo "WARN: sacct returned no rows; retrying in 30s"
    sleep 30
    continue
  fi
  NON_TERMINAL="$(echo "$acct_rows" | awk -F'|' '$1 ~ /^[0-9]+$/ && $2 !~ /^(COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL)/ {print $1":"$2}')"
  if [[ -z "$NON_TERMINAL" ]]; then
    break
  fi
  echo "Still waiting for terminal states:"
  echo "$NON_TERMINAL"
  sleep 30
done

BAD_STATES="$(echo "$acct_rows" | awk -F'|' '$1 ~ /^[0-9]+$/ && $2 !~ /^COMPLETED/ {print $1":"$2}')"
if [[ -n "$BAD_STATES" ]]; then
  echo "ERROR: non-completed jobs detected"
  echo "$BAD_STATES"
  exit 1
fi

echo "Validating expected outputs..."
missing=0
for f in "${EXPECTED[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "MISSING: $f"
    missing=1
  else
    echo "OK: $f"
  fi
done

if [[ "$missing" -ne 0 ]]; then
  echo "ERROR: at least one expected output is missing; not proceeding to merge"
  exit 1
fi

echo "Merging consolidated result tables..."
python3.11 scripts/merge_hellbender_results.py --results-dir results

if command -v pdflatex >/dev/null 2>&1 && [[ -f manuscript/cc_spr_public_substitutes.tex ]]; then
  echo "Compiling manuscript PDF..."
  pdflatex -interaction=nonstopmode -halt-on-error manuscript/cc_spr_public_substitutes.tex >/tmp/ccspr_pdflatex_1.log 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error manuscript/cc_spr_public_substitutes.tex >/tmp/ccspr_pdflatex_2.log 2>&1
  echo "Manuscript compile complete: manuscript/cc_spr_public_substitutes.pdf"
else
  echo "Skipping manuscript compile (pdflatex or tex file missing)."
fi

echo "Done. Next: update results/final_run_report_hellbender.md and results/resume_instructions_hellbender.md with completion status."
