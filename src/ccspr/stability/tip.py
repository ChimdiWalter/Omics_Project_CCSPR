from __future__ import annotations

import numpy as np

from ccspr.geometry.distance import build_distance
from ccspr.solver.cycle import compute_cycle_representative
from ccspr.solver.projection import project_to_features
from ccspr.topology.persistence import build_cycle_problem, compute_persistence, pick_dominant_h1


def _topk_mask(scores: np.ndarray, k: int) -> np.ndarray:
    k = int(k)
    if k <= 0:
        return np.zeros_like(scores, dtype=bool)
    idx = np.argsort(-scores)[:k]
    mask = np.zeros_like(scores, dtype=bool)
    mask[idx] = True
    return mask


def tip_bootstrap_topk(
    x: np.ndarray,
    n_boot: int = 30,
    subsample_frac: float = 0.8,
    top_k: int = 20,
    mode: str = "ricci",
    k: int = 10,
    iters: int = 5,
    alpha: float = 0.5,
    solver: str = "elastic",
    lambda_: float = 1e-6,
    normalize: str = "std",
    rescale_median: float = 1.0,
    seed: int = 0,
    cache_dir: str | None = None,
) -> dict:
    rng = np.random.default_rng(int(seed))
    n_samples, n_features = x.shape
    n_take = max(4, int(subsample_frac * n_samples))

    tip = np.zeros(n_features, dtype=float)
    lifetimes: list[float] = []
    prominences: list[float] = []
    n_ok = 0

    for _ in range(int(n_boot)):
        idx = rng.choice(n_samples, n_take, replace=False)
        x_sub = x[idx]

        try:
            dist = build_distance(
                x_sub,
                mode=mode,
                k=k,
                iters=iters,
                alpha=alpha,
                rescale_median=rescale_median,
                cache_dir=cache_dir,
            )
            ph = compute_persistence(dist, max_dimension=2, cache_dir=cache_dir)
            bars = ph.get("bars", [])
            if not bars:
                continue
            bar = pick_dominant_h1(ph)
            top1 = bars[0]["life"]
            top2 = bars[1]["life"] if len(bars) > 1 else 0.0
            lifetimes.append(float(top1))
            prominences.append(float(max(0.0, top1 - top2)))

            edges, d_mat, z0 = build_cycle_problem(dist, bar)
            if np.count_nonzero(z0) == 0:
                continue

            w = compute_cycle_representative(d_mat, z0, solver=solver, lambda_=lambda_)
            s = project_to_features(x_sub, edges, w, normalize=normalize)
            tip += _topk_mask(s, top_k).astype(float)
            n_ok += 1
        except Exception:
            continue

    den = max(n_ok, 1)
    return {
        "tip": tip / den,
        "n_ok": n_ok,
        "lifetime": np.array(lifetimes, dtype=float),
        "prominence": np.array(prominences, dtype=float),
    }
