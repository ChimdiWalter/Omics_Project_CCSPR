from __future__ import annotations

import numpy as np


def project_to_features(
    x: np.ndarray,
    edges,
    w: np.ndarray,
    normalize: str = "std",
) -> np.ndarray:
    n_features = x.shape[1]
    scores = np.zeros(n_features, dtype=float)

    normalize = str(normalize).lower()
    if normalize not in {"std", "none"}:
        raise ValueError("normalize must be one of {'std', 'none'}")

    if normalize == "std":
        feat_std = np.std(x, axis=0).astype(float)
        feat_std[feat_std < 1e-8] = 1.0
    else:
        feat_std = None

    active = np.where(np.abs(w) > 0)[0]
    for idx in active:
        u, v = edges[idx]
        diff = np.abs(x[u] - x[v]).astype(float)
        if feat_std is not None:
            diff = diff / feat_std
        scores += np.abs(w[idx]) * diff

    den = np.max(scores)
    if den > 0:
        scores = scores / den
    return scores
