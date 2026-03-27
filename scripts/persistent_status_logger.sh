#!/usr/bin/env bash
set -u -o pipefail

JOBS_CSV="${1:-12833143,12833568}"
OUT_FILE="${2:-results/logs/persistent_status_logger_$(date +%Y%m%d_%H%M%S).log}"
SLEEP_SEC="${3:-60}"

mkdir -p "$(dirname "$OUT_FILE")"

IFS=',' read -r -a ROOT_JOBS <<< "$JOBS_CSV"

now_iso() {
  date --iso-8601=seconds
}

is_terminal_state() {
  local s="${1:-UNKNOWN}"
  case "$s" in
    COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

{
  echo "# Persistent CC-SPR Status Log"
  echo "- started: $(now_iso)"
  echo "- jobs: $JOBS_CSV"
  echo "- sleep_sec: $SLEEP_SEC"
  echo "- host: $(hostname)"
  echo
} >> "$OUT_FILE"

while true; do
  {
    echo "## tick $(now_iso)"
    echo "### squeue"
  } >> "$OUT_FILE"

  if ! sq="$(squeue -j "$JOBS_CSV" -o '%.18i %.9P %.25j %.8T %.10M %.6D %R' 2>&1)"; then
    {
      echo "WARN: squeue query failed"
      echo "$sq"
    } >> "$OUT_FILE"
  else
    echo "$sq" >> "$OUT_FILE"
  fi

  {
    echo "### sacct"
  } >> "$OUT_FILE"

  acct_ok=1
  acct_rows=""
  if ! acct_rows="$(sacct -j "$JOBS_CSV" --format=JobID,JobName%30,State,ExitCode,Elapsed -n -P 2>&1)"; then
    acct_ok=0
    {
      echo "WARN: sacct query failed"
      echo "$acct_rows"
      echo
    } >> "$OUT_FILE"
  else
    echo "$acct_rows" >> "$OUT_FILE"
    echo >> "$OUT_FILE"
  fi

  if [[ "$acct_ok" -eq 1 ]]; then
    all_terminal=1
    for jid in "${ROOT_JOBS[@]}"; do
      row="$(echo "$acct_rows" | awk -F'|' -v j="$jid" '$1==j {print; exit}')"
      if [[ -z "$row" ]]; then
        all_terminal=0
        {
          echo "INFO: root job $jid not yet visible in sacct rows as terminal; continue waiting"
        } >> "$OUT_FILE"
        continue
      fi
      state="$(echo "$row" | awk -F'|' '{print $3}')"
      if ! is_terminal_state "$state"; then
        all_terminal=0
      fi
    done

    if [[ "$all_terminal" -eq 1 ]]; then
      {
        echo "All root jobs are in terminal states at $(now_iso)."
        echo "### key outputs snapshot"
        ls -lh \
          results/metrics_cll_rs_scrna_gse165087.csv \
          results/main_results_cll_rs_scrna_gse165087.csv \
          results/ablation_summary_cll_rs_scrna_gse165087.csv \
          results/metrics_all_datasets.csv \
          results/main_results_all_datasets.csv \
          results/ablation_summary_all_datasets.csv \
          manuscript/cc_spr_public_substitutes.pdf 2>&1 || true
        echo
        echo "# logger finished $(now_iso)"
      } >> "$OUT_FILE"
      break
    fi
  fi

  sleep "$SLEEP_SEC"
done
