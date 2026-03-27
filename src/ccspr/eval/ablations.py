from __future__ import annotations

import pandas as pd

from ccspr.eval.classification import FeatureParams, evaluate_ccspr_only


def _single_param(
    base: FeatureParams,
    *,
    mode: str | None = None,
    k: int | None = None,
    iters: int | None = None,
    lambda_: float | None = None,
    top_k: int | None = None,
    normalize: str | None = None,
) -> FeatureParams:
    return FeatureParams(
        mode=base.mode if mode is None else mode,
        solver=base.solver,
        k=base.k if k is None else int(k),
        iters=base.iters if iters is None else int(iters),
        alpha=base.alpha,
        lambda_=base.lambda_ if lambda_ is None else float(lambda_),
        top_k=base.top_k if top_k is None else int(top_k),
        normalize=base.normalize if normalize is None else normalize,
        n_boot=base.n_boot,
        subsample_frac=base.subsample_frac,
        rescale_median=base.rescale_median,
    )


def run_luad_ablations(
    x,
    y,
    base_params: FeatureParams,
    var_genes: int,
    n_splits: int,
    test_size: float,
    seed: int,
    grids: dict | None = None,
    cache_dir: str | None = None,
) -> pd.DataFrame:
    all_rows = []

    default_grids = {
        "distance_mode": ["euclidean", "ricci"],
        "iters": [0, 2, 5, 10],
        "k": [5, 10, 15],
        "lambda": [1e-8, 1e-7, 1e-6, 1e-5, 1e-4],
        "top_k": [10, 20, 50],
        "normalize": ["none", "std"],
    }
    if grids is None:
        grids = default_grids
    else:
        for k, v in default_grids.items():
            grids.setdefault(k, v)

    run_id = 0

    for mode in grids["distance_mode"]:
        params = _single_param(base_params, mode=mode)
        m = evaluate_ccspr_only(
            x,
            y,
            params=params,
            var_genes=var_genes,
            n_splits=n_splits,
            test_size=test_size,
            seed=seed + run_id,
            cache_dir=cache_dir,
        )
        block = m.copy()
        block["ablation"] = "distance_mode"
        block["ablation_value"] = mode
        all_rows.append(block)
        run_id += 1

    for iters in grids["iters"]:
        params = _single_param(base_params, iters=iters)
        m = evaluate_ccspr_only(
            x,
            y,
            params=params,
            var_genes=var_genes,
            n_splits=n_splits,
            test_size=test_size,
            seed=seed + run_id,
            cache_dir=cache_dir,
        )
        block = m.copy()
        block["ablation"] = "iters"
        block["ablation_value"] = iters
        all_rows.append(block)
        run_id += 1

    for k in grids["k"]:
        params = _single_param(base_params, k=k)
        m = evaluate_ccspr_only(
            x,
            y,
            params=params,
            var_genes=var_genes,
            n_splits=n_splits,
            test_size=test_size,
            seed=seed + run_id,
            cache_dir=cache_dir,
        )
        block = m.copy()
        block["ablation"] = "k"
        block["ablation_value"] = k
        all_rows.append(block)
        run_id += 1

    for lam in grids["lambda"]:
        params = _single_param(base_params, lambda_=lam)
        m = evaluate_ccspr_only(
            x,
            y,
            params=params,
            var_genes=var_genes,
            n_splits=n_splits,
            test_size=test_size,
            seed=seed + run_id,
            cache_dir=cache_dir,
        )
        block = m.copy()
        block["ablation"] = "lambda"
        block["ablation_value"] = lam
        all_rows.append(block)
        run_id += 1

    for top_k in grids["top_k"]:
        params = _single_param(base_params, top_k=top_k)
        m = evaluate_ccspr_only(
            x,
            y,
            params=params,
            var_genes=var_genes,
            n_splits=n_splits,
            test_size=test_size,
            seed=seed + run_id,
            cache_dir=cache_dir,
        )
        block = m.copy()
        block["ablation"] = "top_k"
        block["ablation_value"] = top_k
        all_rows.append(block)
        run_id += 1

    for norm in grids["normalize"]:
        params = _single_param(base_params, normalize=norm)
        m = evaluate_ccspr_only(
            x,
            y,
            params=params,
            var_genes=var_genes,
            n_splits=n_splits,
            test_size=test_size,
            seed=seed + run_id,
            cache_dir=cache_dir,
        )
        block = m.copy()
        block["ablation"] = "normalize"
        block["ablation_value"] = norm
        all_rows.append(block)
        run_id += 1

    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
