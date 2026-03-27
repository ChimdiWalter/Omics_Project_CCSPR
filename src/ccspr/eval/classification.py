from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedShuffleSplit

from ccspr.preprocess.basic import select_top_variable_genes, standardize
from ccspr.stability.tip import tip_bootstrap_topk


@dataclass
class FeatureParams:
    mode: str
    solver: str
    k: int
    iters: int
    alpha: float
    lambda_: float
    top_k: int
    normalize: str
    n_boot: int
    subsample_frac: float
    rescale_median: float


def _topk_indices(scores: np.ndarray, k: int) -> np.ndarray:
    k = min(int(k), scores.shape[0])
    return np.argsort(-scores)[:k]


def _fit_predict_logreg(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    seed: int,
) -> np.ndarray:
    clf = LogisticRegression(
        max_iter=5000,
        C=1.0,
        multi_class="auto",
        random_state=int(seed),
        n_jobs=1,
    )
    clf.fit(x_train, y_train)
    return clf.predict(x_test)


def _evaluate_feature_model(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    y_test: np.ndarray,
    params: FeatureParams,
    seed: int,
    cache_dir: str | None,
) -> tuple[float, dict]:
    tip_out = tip_bootstrap_topk(
        x_train,
        n_boot=params.n_boot,
        subsample_frac=params.subsample_frac,
        top_k=params.top_k,
        mode=params.mode,
        k=params.k,
        iters=params.iters,
        alpha=params.alpha,
        solver=params.solver,
        lambda_=params.lambda_,
        normalize=params.normalize,
        rescale_median=params.rescale_median,
        seed=int(seed),
        cache_dir=cache_dir,
    )

    scores = tip_out["tip"]
    idx = _topk_indices(scores, params.top_k)

    xtr = x_train[:, idx]
    xts = x_test[:, idx]
    xtr, xts = standardize(xtr, xts)

    pred = _fit_predict_logreg(xtr, y_train, xts, seed=seed)
    f1 = f1_score(y_test, pred, average="weighted")
    return f1, tip_out


def evaluate_luad_protocol(
    x: np.ndarray,
    y: np.ndarray,
    var_genes: int = 5000,
    n_splits: int = 5,
    test_size: float = 0.25,
    seed: int = 42,
    cache_dir: str | None = None,
    harmonic_params: FeatureParams | None = None,
    eu_sparse_params: FeatureParams | None = None,
    ccspr_params: FeatureParams | None = None,
) -> tuple[pd.DataFrame, dict]:
    if harmonic_params is None:
        harmonic_params = FeatureParams(
            mode="euclidean",
            solver="l2",
            k=10,
            iters=0,
            alpha=0.5,
            lambda_=1e-6,
            top_k=20,
            normalize="std",
            n_boot=20,
            subsample_frac=0.8,
            rescale_median=1.0,
        )
    if ccspr_params is None:
        ccspr_params = FeatureParams(
            mode="ricci",
            solver="elastic",
            k=10,
            iters=5,
            alpha=0.5,
            lambda_=1e-6,
            top_k=20,
            normalize="std",
            n_boot=20,
            subsample_frac=0.8,
            rescale_median=1.0,
        )
    if eu_sparse_params is None:
        eu_sparse_params = FeatureParams(
            mode="euclidean",
            solver="elastic",
            k=ccspr_params.k,
            iters=0,
            alpha=ccspr_params.alpha,
            lambda_=ccspr_params.lambda_,
            top_k=ccspr_params.top_k,
            normalize=ccspr_params.normalize,
            n_boot=ccspr_params.n_boot,
            subsample_frac=ccspr_params.subsample_frac,
            rescale_median=ccspr_params.rescale_median,
        )

    sss = StratifiedShuffleSplit(
        n_splits=int(n_splits),
        test_size=float(test_size),
        random_state=int(seed),
    )

    rows: list[dict] = []
    tips_h = []
    tips_eu_sparse = []
    tips_c = []
    lifetime_h = []
    lifetime_eu_sparse = []
    lifetime_c = []
    prom_h = []
    prom_eu_sparse = []
    prom_c = []

    for split_id, (itr, ite) in enumerate(sss.split(x, y)):
        x_train = x[itr]
        x_test = x[ite]
        y_train = y[itr]
        y_test = y[ite]

        x_train_var, idx_var = select_top_variable_genes(x_train, top_n=int(var_genes))
        x_test_var = x_test[:, idx_var]

        xtr_std, xts_std = standardize(x_train_var, x_test_var)
        y_pred_std = _fit_predict_logreg(xtr_std, y_train, xts_std, seed + split_id)
        f1_std = f1_score(y_test, y_pred_std, average="weighted")

        rows.append(
            {
                "split": split_id,
                "model": "standard",
                "f1_weighted": float(f1_std),
                "distance_mode": "euclidean",
                "solver": "none",
            }
        )

        f1_h, out_h = _evaluate_feature_model(
            x_train_var,
            y_train,
            x_test_var,
            y_test,
            params=harmonic_params,
            seed=seed + 1000 + split_id,
            cache_dir=cache_dir,
        )
        rows.append(
            {
                "split": split_id,
                "model": "harmonic",
                "f1_weighted": float(f1_h),
                "distance_mode": harmonic_params.mode,
                "solver": harmonic_params.solver,
                "k": harmonic_params.k,
                "iters": harmonic_params.iters,
                "lambda": harmonic_params.lambda_,
                "top_k": harmonic_params.top_k,
                "normalize": harmonic_params.normalize,
            }
        )

        f1_eu_sparse, out_eu_sparse = _evaluate_feature_model(
            x_train_var,
            y_train,
            x_test_var,
            y_test,
            params=eu_sparse_params,
            seed=seed + 1500 + split_id,
            cache_dir=cache_dir,
        )
        rows.append(
            {
                "split": split_id,
                "model": "eu_sparse",
                "f1_weighted": float(f1_eu_sparse),
                "distance_mode": eu_sparse_params.mode,
                "solver": eu_sparse_params.solver,
                "k": eu_sparse_params.k,
                "iters": eu_sparse_params.iters,
                "lambda": eu_sparse_params.lambda_,
                "top_k": eu_sparse_params.top_k,
                "normalize": eu_sparse_params.normalize,
            }
        )

        f1_c, out_c = _evaluate_feature_model(
            x_train_var,
            y_train,
            x_test_var,
            y_test,
            params=ccspr_params,
            seed=seed + 2000 + split_id,
            cache_dir=cache_dir,
        )
        rows.append(
            {
                "split": split_id,
                "model": "ccspr",
                "f1_weighted": float(f1_c),
                "distance_mode": ccspr_params.mode,
                "solver": ccspr_params.solver,
                "k": ccspr_params.k,
                "iters": ccspr_params.iters,
                "lambda": ccspr_params.lambda_,
                "top_k": ccspr_params.top_k,
                "normalize": ccspr_params.normalize,
            }
        )

        tips_h.append(out_h["tip"])
        tips_eu_sparse.append(out_eu_sparse["tip"])
        tips_c.append(out_c["tip"])
        lifetime_h.extend(out_h["lifetime"].tolist())
        lifetime_eu_sparse.extend(out_eu_sparse["lifetime"].tolist())
        lifetime_c.extend(out_c["lifetime"].tolist())
        prom_h.extend(out_h["prominence"].tolist())
        prom_eu_sparse.extend(out_eu_sparse["prominence"].tolist())
        prom_c.extend(out_c["prominence"].tolist())

    aux = {
        "tip_eu": np.mean(tips_eu_sparse, axis=0) if tips_eu_sparse else np.zeros(x.shape[1]),
        "tip_ricci": np.mean(tips_c, axis=0) if tips_c else np.zeros(x.shape[1]),
        "tip_harmonic": np.mean(tips_h, axis=0) if tips_h else np.zeros(x.shape[1]),
        "lifetime_eu": np.array(lifetime_eu_sparse, dtype=float),
        "lifetime_ricci": np.array(lifetime_c, dtype=float),
        "lifetime_harmonic": np.array(lifetime_h, dtype=float),
        "prominence_eu": np.array(prom_eu_sparse, dtype=float),
        "prominence_ricci": np.array(prom_c, dtype=float),
        "prominence_harmonic": np.array(prom_h, dtype=float),
    }
    return pd.DataFrame(rows), aux


def evaluate_ccspr_only(
    x: np.ndarray,
    y: np.ndarray,
    params: FeatureParams,
    var_genes: int = 5000,
    n_splits: int = 3,
    test_size: float = 0.25,
    seed: int = 42,
    cache_dir: str | None = None,
) -> pd.DataFrame:
    sss = StratifiedShuffleSplit(
        n_splits=int(n_splits),
        test_size=float(test_size),
        random_state=int(seed),
    )

    rows: list[dict] = []
    for split_id, (itr, ite) in enumerate(sss.split(x, y)):
        x_train = x[itr]
        x_test = x[ite]
        y_train = y[itr]
        y_test = y[ite]

        x_train_var, idx_var = select_top_variable_genes(x_train, top_n=int(var_genes))
        x_test_var = x_test[:, idx_var]

        f1_c, _ = _evaluate_feature_model(
            x_train_var,
            y_train,
            x_test_var,
            y_test,
            params=params,
            seed=seed + split_id,
            cache_dir=cache_dir,
        )
        rows.append(
            {
                "split": split_id,
                "model": "ccspr",
                "f1_weighted": float(f1_c),
                "distance_mode": params.mode,
                "solver": params.solver,
                "k": params.k,
                "iters": params.iters,
                "lambda": params.lambda_,
                "top_k": params.top_k,
                "normalize": params.normalize,
            }
        )
    return pd.DataFrame(rows)
