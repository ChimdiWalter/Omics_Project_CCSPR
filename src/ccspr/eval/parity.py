from __future__ import annotations

import subprocess
from pathlib import Path


def run_matilda_parity_hook(
    r_script: str | Path,
    input_csv: str | Path,
    output_csv: str | Path,
) -> int:
    """Optional parity hook to compare outputs against an R/maTilDA script."""
    cmd = ["Rscript", str(r_script), str(input_csv), str(output_csv)]
    proc = subprocess.run(cmd, check=False)
    return int(proc.returncode)
