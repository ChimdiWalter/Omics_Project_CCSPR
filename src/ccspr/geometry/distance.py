from __future__ import annotations

import hashlib
import logging

import networkx as nx
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.linalg import eigsh
from sklearn.neighbors import kneighbors_graph, NearestNeighbors

from ccspr.utils.cache import cache_path, load_npy, save_npy
from ccspr.utils.io import hash_array

LOGGER = logging.getLogger(__name__)

try:
    from GraphRicciCurvature.OllivierRicci import OllivierRicci

    HAS_RICCI = True
except Exception:
    HAS_RICCI = False


def _cache_key(
    x: np.ndarray,
    mode: str,
    k: int,
    iters: int,
    alpha: float,
    rescale_median: float,
) -> str:
    h = hashlib.sha1()
    h.update(hash_array(x).encode("utf-8"))
    h.update(f"{mode}|{k}|{iters}|{alpha}|{rescale_median}".encode("utf-8"))
    return h.hexdigest()


def _rescale_distance_matrix(dist: np.ndarray, target_median: float = 1.0) -> np.ndarray:
    if target_median <= 0:
        return dist
    tri = dist[np.triu_indices_from(dist, 1)]
    tri = tri[np.isfinite(tri)]
    tri = tri[tri > 0]
    if tri.size == 0:
        return dist
    med = np.median(tri)
    if med <= 0:
        return dist
    return dist * (target_median / med)


def _ensure_connected_knn_graph(g: nx.Graph, dist_euclid: np.ndarray) -> nx.Graph:
    comps = [list(c) for c in nx.connected_components(g)]
    if len(comps) <= 1:
        return g

    while len(comps) > 1:
        c0 = comps[0]
        best = (None, None, np.inf)
        for i in range(1, len(comps)):
            ci = comps[i]
            for u in c0:
                drow = dist_euclid[u, ci]
                j = int(np.argmin(drow))
                v = ci[j]
                d = float(drow[j])
                if d < best[2]:
                    best = (u, v, d)
        u, v, d = best
        if u is None:
            break
        g.add_edge(int(u), int(v), weight=max(d, 1e-9))
        comps = [list(c) for c in nx.connected_components(g)]
    return g


def _euclidean_distance(x: np.ndarray) -> np.ndarray:
    dist = squareform(pdist(x, metric="euclidean")).astype(np.float64)
    np.fill_diagonal(dist, 0.0)
    return dist


def _ricci_distance(
    x: np.ndarray,
    k: int,
    iters: int,
    alpha: float,
) -> np.ndarray:
    dist_euclid = _euclidean_distance(x)

    if not HAS_RICCI:
        LOGGER.warning("GraphRicciCurvature not installed; falling back to Euclidean distance.")
        return dist_euclid

    a = kneighbors_graph(
        x,
        n_neighbors=min(int(k), x.shape[0] - 1),
        mode="distance",
        include_self=False,
    )
    a = a.minimum(a.T)
    g = nx.from_scipy_sparse_array(a, edge_attribute="weight")

    for _, _, data in g.edges(data=True):
        data["weight"] = max(float(data.get("weight", 1.0)), 1e-9)

    g = _ensure_connected_knn_graph(g, dist_euclid)

    try:
        orc = OllivierRicci(g, alpha=float(alpha), verbose="ERROR", proc=1)
    except TypeError:
        orc = OllivierRicci(g, alpha=float(alpha), verbose="ERROR")

    if int(iters) > 0:
        orc.compute_ricci_flow(iterations=int(iters))

    dist_ricci = nx.floyd_warshall_numpy(orc.G, weight="weight").astype(np.float64)
    inf_mask = ~np.isfinite(dist_ricci)
    if np.any(inf_mask):
        max_finite = np.max(dist_ricci[np.isfinite(dist_ricci)])
        dist_ricci[inf_mask] = np.maximum(dist_euclid[inf_mask], max_finite)

    dist_ricci = 0.5 * (dist_ricci + dist_ricci.T)
    np.fill_diagonal(dist_ricci, 0.0)
    return dist_ricci


def _build_knn_affinity(x: np.ndarray, k: int) -> np.ndarray:
    """Build symmetric kNN affinity matrix with adaptive Gaussian kernel."""
    n = x.shape[0]
    k_use = min(int(k), n - 1)
    nn = NearestNeighbors(n_neighbors=k_use + 1, metric="euclidean").fit(x)
    dists, indices = nn.kneighbors(x)
    # sigma_i = distance to k-th neighbor (adaptive bandwidth)
    sigma = dists[:, -1].clip(min=1e-10)
    # Build sparse affinity
    W = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j_idx in range(1, k_use + 1):  # skip self
            j = indices[i, j_idx]
            d = dists[i, j_idx]
            w = np.exp(-d ** 2 / (sigma[i] * sigma[j] + 1e-12))
            W[i, j] = max(W[i, j], w)
            W[j, i] = max(W[j, i], w)
    return W


def _diffusion_distance(x: np.ndarray, k: int, iters: int) -> np.ndarray:
    """Diffusion distance via truncated eigendecomposition of Markov transition matrix.

    Parameters
    ----------
    x : data matrix (n, p)
    k : kNN neighbors for affinity graph
    iters : diffusion time parameter t (reuses 'iters' field)
    """
    t = max(1, int(iters))
    W = _build_knn_affinity(x, k)
    # Row-normalize to transition matrix
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    P = W / row_sums

    n = P.shape[0]
    n_eig = min(50, n - 2)

    # Eigendecomposition of P (real, but may not be symmetric)
    # Symmetrize: P_sym = D^{1/2} P D^{-1/2} has same eigenvalues
    d_sqrt = np.sqrt(row_sums.ravel())
    d_sqrt_inv = np.where(d_sqrt > 0, 1.0 / d_sqrt, 0.0)
    P_sym = (d_sqrt[:, None] * P) * d_sqrt_inv[None, :]
    P_sym = 0.5 * (P_sym + P_sym.T)  # ensure symmetric

    try:
        eigenvalues, eigenvectors = eigsh(P_sym, k=n_eig, which="LM")
    except Exception:
        LOGGER.warning("Eigendecomposition failed; falling back to Euclidean distance.")
        return _euclidean_distance(x)

    # Sort descending
    order = np.argsort(-eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    # Convert back to right eigenvectors of P: psi = D^{-1/2} * phi
    psi = eigenvectors * d_sqrt_inv[:, None]

    # Diffusion distance: D_t(i,j)^2 = sum_l (lambda_l^t)^2 * (psi_l(i) - psi_l(j))^2
    # Skip the trivial eigenvalue (lambda_0 ~ 1)
    lam_t = np.abs(eigenvalues[1:]) ** t
    coords = psi[:, 1:] * lam_t[None, :]

    dist = squareform(pdist(coords, metric="euclidean")).astype(np.float64)
    np.fill_diagonal(dist, 0.0)
    return dist


def _phate_like_distance(x: np.ndarray, k: int, iters: int) -> np.ndarray:
    """PHATE-like potential distance via information-geometric divergence.

    Computes potential U_t(i) = -log(P^t(i,:) + eps) and returns
    D_pot(i,j) = ||U_t(i) - U_t(j)||_2. This is a lightweight surrogate
    for full PHATE (avoids the MDS embedding step).

    Parameters
    ----------
    x : data matrix (n, p)
    k : kNN neighbors for affinity graph
    iters : diffusion time parameter t (reuses 'iters' field)
    """
    t = max(1, int(iters))
    W = _build_knn_affinity(x, k)
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    P = W / row_sums

    # Power P^t via repeated matrix multiplication (t is small)
    Pt = np.linalg.matrix_power(P, t)

    # Potential: U_t(i) = -log(P^t(i,:) + eps)
    eps = 1e-7 / P.shape[0]
    U = -np.log(Pt + eps)

    dist = squareform(pdist(U, metric="euclidean")).astype(np.float64)
    np.fill_diagonal(dist, 0.0)
    return dist


def _dtm_distance(x: np.ndarray, k: int) -> np.ndarray:
    """Distance-to-measure (DTM) weighted distance.

    For each point i, the DTM weight is:
        d_DTM(i) = sqrt(1/k * sum_{j in kNN(i)} d(i,j)^2)

    The DTM-weighted pairwise distance is:
        D_DTM(i,j) = sqrt(d(i,j)^2 + d_DTM(i)^2 + d_DTM(j)^2)

    This makes the metric robust to outliers.

    Parameters
    ----------
    x : data matrix (n, p)
    k : kNN neighbors for DTM computation
    """
    n = x.shape[0]
    k_use = min(int(k), n - 1)

    dist_euclid = _euclidean_distance(x)

    nn = NearestNeighbors(n_neighbors=k_use + 1, metric="euclidean").fit(x)
    dists_knn, _ = nn.kneighbors(x)
    # DTM weight per point (skip self-distance at index 0)
    dtm_sq = np.mean(dists_knn[:, 1:] ** 2, axis=1)

    # DTM-weighted distance matrix
    dist_dtm = np.sqrt(
        dist_euclid ** 2 + dtm_sq[:, None] + dtm_sq[None, :]
    ).astype(np.float64)
    np.fill_diagonal(dist_dtm, 0.0)
    return dist_dtm


VALID_MODES = {"euclidean", "ricci", "diffusion", "phate_like", "dtm"}


def build_distance(
    x: np.ndarray,
    mode: str = "euclidean",
    k: int = 10,
    iters: int = 10,
    alpha: float = 0.5,
    rescale_median: float = 1.0,
    cache_dir: str | None = None,
) -> np.ndarray:
    mode = str(mode).lower()
    if mode not in VALID_MODES:
        raise ValueError(f"mode must be one of {VALID_MODES}")

    key = _cache_key(
        x=x,
        mode=mode,
        k=int(k),
        iters=int(iters),
        alpha=float(alpha),
        rescale_median=float(rescale_median),
    )

    if cache_dir:
        cp = cache_path(cache_dir, "distance", key, "npy")
        cached = load_npy(cp)
        if cached is not None:
            return cached

    if mode == "euclidean":
        dist = _euclidean_distance(x)
    elif mode == "ricci":
        dist = _ricci_distance(x, k=int(k), iters=int(iters), alpha=float(alpha))
    elif mode == "diffusion":
        dist = _diffusion_distance(x, k=int(k), iters=int(iters))
    elif mode == "phate_like":
        dist = _phate_like_distance(x, k=int(k), iters=int(iters))
    elif mode == "dtm":
        dist = _dtm_distance(x, k=int(k))
    else:
        raise ValueError(f"Unknown mode: {mode}")

    dist = _rescale_distance_matrix(dist, target_median=float(rescale_median))

    if cache_dir:
        save_npy(cp, dist)
    return dist
