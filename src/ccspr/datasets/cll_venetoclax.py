from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_cll_venetoclax(
    processed_matrix_path: str | Path = "data/cll_venetoclax/processed_matrix.tsv",
    label_col: str = "label",
) -> dict:
    path = Path(processed_matrix_path)
    if not path.exists():
        raise FileNotFoundError(
            "CLL Venetoclax processed matrix not found. "
            "Expected file at data/cll_venetoclax/processed_matrix.tsv with columns: "
            "sample_id, label, feature_*. Add this file then rerun. "
            "Primary source study is request-based (Khalsa et al., Blood 2023, "
            "doi:10.1182/blood.2022016600)."
        )

    df = pd.read_csv(path, sep="\t")
    if "sample_id" not in df.columns or label_col not in df.columns:
        raise ValueError("processed_matrix.tsv must include 'sample_id' and label column.")

    feature_cols = [c for c in df.columns if c not in {"sample_id", label_col}]
    x = df[feature_cols].to_numpy(dtype=float)
    y = df[label_col].astype(str).to_numpy()

    return {
        "X": np.nan_to_num(x),
        "y": y,
        "sample_ids": df["sample_id"].astype(str).to_numpy(),
        "feature_names": np.array(feature_cols, dtype=object),
        "meta": {"placeholder": False, "source": str(path)},
    }
