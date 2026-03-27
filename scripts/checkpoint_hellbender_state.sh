#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

JOBS="${1:-12833129,12833143}"
SESSION="${2:-ccspr}"
mkdir -p results/logs

TS="$(date +%Y%m%d_%H%M%S)"
OUT="results/logs/checkpoint_hellbender_${TS}.md"

{
  echo "# Hellbender Checkpoint"
  echo
  echo "- timestamp: $(date --iso-8601=seconds)"
  echo "- repo: $REPO_DIR"
  echo "- jobs: $JOBS"
  echo "- tmux session: $SESSION"
  echo
  echo "## Queue"
  squeue -j "$JOBS" -o '%.18i %.9P %.25j %.8T %.10M %.6D %R' || true
  echo
  echo "## Accounting"
  sacct -j "$JOBS" --format=JobID,JobName%30,State,ExitCode,Elapsed -n -P || true
  echo
  echo "## tmux"
  if tmux has-session -t "$SESSION" 2>/dev/null; then
    tmux list-windows -t "$SESSION"
  else
    echo "session '$SESSION' not found"
  fi
  echo
  echo "## Current Result Files"
  find results -maxdepth 2 -type f | sort
} > "$OUT"

echo "$OUT"
