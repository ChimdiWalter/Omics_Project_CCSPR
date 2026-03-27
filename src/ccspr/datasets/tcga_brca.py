from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from ccspr.datasets.common import download_first_available, normalize_sample_id, read_tsv_auto
from ccspr.utils.io import ensure_dir


def download_tcga_brca(root: str | Path = "data/tcga_brca") -> dict:
    root = ensure_dir(root)
    expr_path = root / "expression.tsv"
    meth_path = root / "methylation.tsv"
    clin_path = root / "clinical.tsv"

    expr_url = download_first_available(
        [
            "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
            "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_mRNA_median_Zscores.txt",
        ],
        expr_path,
    )
    meth_url = download_first_available(
        [
            "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HumanMethylation27.gz",
            "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HumanMethylation450.gz",
            "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_tcga_pan_can_atlas_2018/data_methylation_hm27_hm450_merged.txt",
        ],
        meth_path,
    )
    clin_url = download_first_available(
        [
            "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
            "https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
            "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
        ],
        clin_path,
    )

    return {
        "expression_path": str(expr_path),
        "methylation_path": str(meth_path),
        "clinical_path": str(clin_path),
        "expression_url": expr_url,
        "methylation_url": meth_url,
        "clinical_url": clin_url,
    }


def _parse_gene_by_sample(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cols = [str(c) for c in df.columns]
    if len(cols) >= 3 and cols[0].lower().startswith("hugo"):
        samples = np.array([normalize_sample_id(c) for c in cols[2:]], dtype=object)
        feats = df.iloc[:, 0].astype(str).to_numpy()
        x = df.iloc[:, 2:].to_numpy(dtype=float).T
        return x, samples, feats

    first = cols[0]
    feats = df[first].astype(str).to_numpy()
    x = df.drop(columns=[first]).to_numpy(dtype=float).T
    samples = np.array([normalize_sample_id(c) for c in cols[1:]], dtype=object)
    return x, samples, feats


def _load_labels(df: pd.DataFrame) -> pd.DataFrame:
    id_candidates = ["sampleID", "SAMPLE_ID", "sample", "PATIENT_ID", "patient_id"]
    id_col = next((c for c in id_candidates if c in df.columns), None)
    if id_col is None:
        id_col = df.columns[0]

    for c in [
        "PAM50Call_RNAseq",
        "BRCA_Subtype_PAM50",
        "PAM50",
        "PAM50 mRNA",
        "SUBTYPE",
        "Subtype",
        "Subtype_mRNA",
    ]:
        if c in df.columns:
            lab_col = c
            break
    else:
        non_id = [c for c in df.columns if c != id_col]
        lab_col = None
        for c in non_id:
            if df[c].dtype == object and 2 <= df[c].nunique(dropna=True) <= 20:
                lab_col = c
                break
        if lab_col is None:
            lab_col = non_id[0]

    out = df[[id_col, lab_col]].copy()
    out.columns = ["sample_id", "label"]
    out["sample_id"] = out["sample_id"].astype(str).map(normalize_sample_id)
    out = out.dropna().drop_duplicates(subset=["sample_id"])
    return out


def load_tcga_brca_multiomics(
    root: str | Path = "data/tcga_brca",
    pca_expr: int = 100,
    pca_meth: int = 100,
    min_class_size: int = 10,
    max_samples: int | None = None,
    seed: int = 42,
) -> dict:
    root = Path(root)
    expr_path = root / "expression.tsv"
    meth_path = root / "methylation.tsv"
    clin_path = root / "clinical.tsv"

    if not expr_path.exists() or not meth_path.exists() or not clin_path.exists():
        download_tcga_brca(root)

    expr_df = read_tsv_auto(expr_path)
    meth_df = read_tsv_auto(meth_path)
    clin_df = read_tsv_auto(clin_path)

    x_expr, s_expr, _ = _parse_gene_by_sample(expr_df)
    x_meth, s_meth, _ = _parse_gene_by_sample(meth_df)

    df_expr = pd.DataFrame({"sample_id": s_expr, "i": np.arange(len(s_expr))})
    df_meth = pd.DataFrame({"sample_id": s_meth, "j": np.arange(len(s_meth))})
    labels = _load_labels(clin_df)

    m = df_expr.merge(df_meth, on="sample_id", how="inner").merge(labels, on="sample_id", how="inner")
    if m.empty:
        raise ValueError("No overlapping BRCA samples across expression/methylation/clinical data.")

    xe = x_expr[m["i"].to_numpy()]
    xm = x_meth[m["j"].to_numpy()]

    xe = np.nan_to_num(xe)
    xm = np.nan_to_num(xm)

    pe = PCA(n_components=min(int(pca_expr), xe.shape[0] - 1, xe.shape[1]), random_state=0)
    pm = PCA(n_components=min(int(pca_meth), xm.shape[0] - 1, xm.shape[1]), random_state=0)
    ze = pe.fit_transform(xe)
    zm = pm.fit_transform(xm)

    x = np.concatenate([ze, zm], axis=1)
    y = m["label"].astype(str).to_numpy()

    classes, counts = np.unique(y, return_counts=True)
    keep_classes = set(classes[counts >= int(min_class_size)])
    keep = np.array([c in keep_classes for c in y])

    x = x[keep]
    y = y[keep]
    sample_ids = m["sample_id"].astype(str).to_numpy()[keep]

    if max_samples is not None and x.shape[0] > int(max_samples):
        rng = np.random.default_rng(int(seed))
        idx = np.arange(x.shape[0])
        classes, counts = np.unique(y, return_counts=True)
        take_idx = []
        for c, n in zip(classes, counts):
            cidx = idx[y == c]
            n_take = max(2, int(round((n / len(y)) * int(max_samples))))
            n_take = min(n_take, len(cidx))
            take_idx.extend(rng.choice(cidx, n_take, replace=False).tolist())
        take_idx = np.array(sorted(set(take_idx)), dtype=int)
        if len(take_idx) > int(max_samples):
            take_idx = rng.choice(take_idx, int(max_samples), replace=False)
        x = x[take_idx]
        y = y[take_idx]
        sample_ids = sample_ids[take_idx]

    return {
        "X": x,
        "y": y,
        "sample_ids": sample_ids,
        "feature_names": np.array([f"expr_pca_{i}" for i in range(ze.shape[1])] + [f"meth_pca_{i}" for i in range(zm.shape[1])]),
        "meta": {
            "n_samples": int(x.shape[0]),
            "n_features": int(x.shape[1]),
            "classes": sorted(list(keep_classes)),
            "expr_pca": int(ze.shape[1]),
            "meth_pca": int(zm.shape[1]),
        },
    }
