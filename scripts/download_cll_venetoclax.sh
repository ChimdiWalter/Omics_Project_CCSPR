#!/usr/bin/env bash
set -euo pipefail
source ~/.venvs/lesegenv/bin/activate
export PYTHONPATH=src
python3 scripts/download_cll_venetoclax.py
