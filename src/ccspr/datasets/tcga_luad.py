from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccspr.datasets.common import (
    download_first_available,
    normalize_sample_id,
    read_tsv_auto,
)
from ccspr.utils.io import ensure_dir


def download_tcga_luad(root: str | Path = "data/tcga_luad") -> dict:
    root = ensure_dir(root)
    expr_path = root / "expression.tsv.gz"
    labels_path = root / "labels.tsv.gz"

    expr_url = download_first_available(
        [
            "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
            "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/luad_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_mRNA_median_Zscores.txt",
        ],
        expr_path,
    )
    labels_url = download_first_available(
        [
            "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix",
            "https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
            "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/luad_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
            "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
        ],
        labels_path,
    )

    return {
        "expression_path": str(expr_path),
        "labels_path": str(labels_path),
        "expression_url": expr_url,
        "labels_url": labels_url,
    }


def _parse_expression(expr_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cols = [str(c) for c in expr_df.columns]

    if "sample" in expr_df.columns:
        # Xena sometimes uses "sample" as the gene-symbol column in gene x sample matrices.
        sample_like = expr_df["sample"].astype(str).str.startswith("TCGA-").mean()
        if sample_like > 0.5:
            sample_ids = expr_df["sample"].astype(str).map(normalize_sample_id).values
            x = expr_df.drop(columns=["sample"]).to_numpy(dtype=float)
            feature_names = expr_df.drop(columns=["sample"]).columns.astype(str).values
            return x, sample_ids, feature_names
        feature_names = expr_df["sample"].astype(str).values
        sample_ids = np.array([normalize_sample_id(c) for c in cols[1:]], dtype=object)
        x = expr_df.iloc[:, 1:].to_numpy(dtype=float).T
        return x, sample_ids, feature_names

    if len(cols) >= 3 and cols[0].lower().startswith("hugo"):
        sample_cols = cols[2:]
        feature_names = expr_df.iloc[:, 0].astype(str).values
        x = expr_df.iloc[:, 2:].to_numpy(dtype=float).T
        sample_ids = np.array([normalize_sample_id(c) for c in sample_cols], dtype=object)
        return x, sample_ids, feature_names

    first_col = cols[0]
    if first_col.lower() in {"gene", "genes", "symbol"}:
        feature_names = expr_df[first_col].astype(str).values
        x = expr_df.drop(columns=[first_col]).to_numpy(dtype=float).T
        sample_ids = np.array([normalize_sample_id(c) for c in cols[1:]], dtype=object)
        return x, sample_ids, feature_names

    raise ValueError("Could not parse LUAD expression format.")


def _parse_labels(labels_df: pd.DataFrame) -> pd.DataFrame:
    cands_id = [
        "sample",
        "Sample ID",
        "sampleID",
        "SAMPLE_ID",
        "PATIENT_ID",
        "bcr_patient_barcode",
        "submitter_id.samples",
    ]
    cands_lab = [
        "Subtype_mRNA",
        "subtype",
        "SUBTYPE",
        "Expression_Subtype",
        "CONSENSUSMCLUST",
        "mRNA_cluster",
        "PAM50",
        "Subtype",
        "expression_subtype",
        "sample_type",
        "SAMPLE_TYPE",
        "histological_type",
    ]

    id_col = next((c for c in cands_id if c in labels_df.columns), None)
    lab_col = next((c for c in cands_lab if c in labels_df.columns), None)

    if id_col is None:
        id_col = labels_df.columns[0]

    if lab_col is None:
        non_id = [c for c in labels_df.columns if c != id_col]
        for c in non_id:
            if labels_df[c].dtype == object and labels_df[c].nunique(dropna=True) >= 2:
                lab_col = c
                break

    if lab_col is None:
        raise ValueError(
            "Could not find subtype label column. Provide labels file with one of: "
            f"{cands_lab}."
        )

    out = labels_df[[id_col, lab_col]].copy()
    out.columns = ["sample_id", "label"]
    out["sample_id"] = out["sample_id"].astype(str).map(normalize_sample_id)
    out["label"] = out["label"].astype(str)
    out = out.replace({"label": {"nan": np.nan, "[Not Available]": np.nan, "": np.nan}})
    out = out.dropna(subset=["label"]).drop_duplicates(subset=["sample_id"])
    return out


def _fallback_label_from_sample_id(sample_ids: np.ndarray) -> pd.DataFrame:
    def _map_code(s: str) -> str:
        parts = str(s).split("-")
        if len(parts) < 4:
            return "unknown"
        code = parts[3][:2]
        if code == "01":
            return "Primary Tumor"
        if code == "11":
            return "Solid Tissue Normal"
        return f"SampleType_{code}"

    out = pd.DataFrame({"sample_id": sample_ids.astype(str)})
    out["label"] = out["sample_id"].map(_map_code)
    return out.drop_duplicates(subset=["sample_id"])


def load_tcga_luad(
    root: str | Path = "data/tcga_luad",
    expression_path: str | None = None,
    labels_path: str | None = None,
    min_class_size: int = 10,
) -> dict:
    root = Path(root)
    expression_path = root / "expression.tsv.gz" if expression_path is None else Path(expression_path)
    labels_path = root / "labels.tsv.gz" if labels_path is None else Path(labels_path)

    if not expression_path.exists() or not labels_path.exists():
        download_tcga_luad(root)

    expr_df = read_tsv_auto(expression_path)
    try:
        labels_df = read_tsv_auto(labels_path)
    except Exception:
        labels_df = pd.read_csv(labels_path, sep="\t", comment="#", low_memory=False, compression=None)
    if labels_df.empty:
        labels_df = pd.read_csv(labels_path, sep="\t", comment="#", low_memory=False, compression=None)

    x, sample_ids, feature_names = _parse_expression(expr_df)
    try:
        labels = _parse_labels(labels_df)
    except Exception:
        labels = _fallback_label_from_sample_id(sample_ids)
    if "_primary_disease" in labels_df.columns:
        labels_df = labels_df.copy()
        labels_df["_primary_disease"] = labels_df["_primary_disease"].astype(str).str.lower()
        luad_ids = set(
            labels_df.loc[
                labels_df["_primary_disease"].str.contains("lung adenocarcinoma", na=False),
                labels_df.columns[0],
            ]
            .astype(str)
            .map(normalize_sample_id)
        )
        if luad_ids:
            labels = labels[labels["sample_id"].isin(luad_ids)].copy()

    df_samples = pd.DataFrame({"sample_id": sample_ids, "idx": np.arange(len(sample_ids))})
    merged = df_samples.merge(labels, on="sample_id", how="inner")
    if merged["label"].nunique() < 2:
        merged = pd.DataFrame(columns=["sample_id", "idx", "label"])
    if merged.empty:
        merged = df_samples.merge(_fallback_label_from_sample_id(sample_ids), on="sample_id", how="inner")

    if merged.empty:
        raise ValueError(
            "No overlapping samples between expression and labels. "
            "Check source files or provide explicit expression/labels paths."
        )

    x_mat = x[merged["idx"].to_numpy()]
    y = merged["label"].astype(str).to_numpy()
    sids = merged["sample_id"].astype(str).to_numpy()

    classes, counts = np.unique(y, return_counts=True)
    keep_classes = set(classes[counts >= int(min_class_size)])
    keep = np.array([c in keep_classes for c in y])

    return {
        "X": x_mat[keep].astype(np.float64),
        "y": y[keep],
        "sample_ids": sids[keep],
        "feature_names": feature_names.astype(str),
        "meta": {
            "n_samples": int(np.sum(keep)),
            "n_features": int(x_mat.shape[1]),
            "classes": sorted(list(keep_classes)),
            "expression_path": str(expression_path),
            "labels_path": str(labels_path),
        },
    }
