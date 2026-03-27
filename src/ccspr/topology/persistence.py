from __future__ import annotations

import hashlib
from typing import Any

import gudhi
import networkx as nx
import numpy as np
from scipy.sparse import dok_matrix

from ccspr.utils.cache import cache_path, load_pickle, save_pickle
from ccspr.utils.io import hash_array


def _cache_key(dist: np.ndarray, max_dimension: int) -> str:
    h = hashlib.sha1()
    h.update(hash_array(dist).encode("utf-8"))
    h.update(f"{max_dimension}".encode("utf-8"))
    return h.hexdigest()


def _extract_h1_bars(simplex_tree: gudhi.SimplexTree) -> list[dict[str, Any]]:
    simplex_tree.persistence()
    pairs = simplex_tree.persistence_pairs()
    bars: list[dict[str, Any]] = []
    for b, d in pairs:
        if len(b) == 2 and len(d) == 3:
            bt = float(simplex_tree.filtration(b))
            dt = float(simplex_tree.filtration(d))
            if np.isfinite(dt) and dt > bt:
                bars.append(
                    {
                        "birth_edge": tuple(sorted(map(int, b))),
                        "death_tri": tuple(sorted(map(int, d))),
                        "bt": bt,
                        "dt": dt,
                        "life": dt - bt,
                    }
                )
    bars.sort(key=lambda x: x["life"], reverse=True)
    return bars


def compute_persistence(
    dist: np.ndarray,
    max_dimension: int = 2,
    cache_dir: str | None = None,
) -> dict[str, Any]:
    key = _cache_key(dist=dist, max_dimension=int(max_dimension))
    if cache_dir:
        cp = cache_path(cache_dir, "persistence", key, "pkl")
        cached = load_pickle(cp)
        if cached is not None:
            return cached

    rips = gudhi.RipsComplex(distance_matrix=dist)
    st = rips.create_simplex_tree(max_dimension=int(max_dimension))
    bars = _extract_h1_bars(st)
    result = {"bars": bars, "max_dimension": int(max_dimension)}

    if cache_dir:
        save_pickle(cp, result)
    return result


def pick_dominant_h1(persistence: dict[str, Any]) -> dict[str, Any] | None:
    bars = persistence.get("bars", [])
    if not bars:
        return None
    return bars[0]


def top_h1_bars(persistence: dict[str, Any], top_k: int = 50) -> list[dict[str, Any]]:
    bars = persistence.get("bars", [])
    return bars[: int(top_k)]


def build_boundary_at_cutoff(simplex_tree: gudhi.SimplexTree, cutoff: float):
    edges: list[tuple[int, int]] = []
    triangles: list[tuple[int, int, int]] = []
    edge_index: dict[tuple[int, int], int] = {}

    for simplex, val in simplex_tree.get_filtration():
        if val >= cutoff:
            continue
        if len(simplex) == 2:
            u, v = sorted(simplex)
            e = (int(u), int(v))
            if e not in edge_index:
                edge_index[e] = len(edges)
                edges.append(e)
        elif len(simplex) == 3:
            triangles.append(tuple(sorted(map(int, simplex))))

    m, n = len(edges), len(triangles)
    d_mat = dok_matrix((m, n), dtype=np.float64)

    for j, (a, b, c) in enumerate(triangles):
        for (u, v), sgn in [((b, c), +1), ((a, c), -1), ((a, b), +1)]:
            if (u, v) in edge_index:
                d_mat[edge_index[(u, v)], j] = sgn

    return edges, triangles, edge_index, d_mat.tocsc()


def construct_z0_from_birth_mst(
    birth_edge,
    edge_index,
    edges,
    dist_matrix,
    bt,
    eps: float = 1e-12,
):
    u0, v0 = map(int, birth_edge)
    u, v = min(u0, v0), max(u0, v0)

    g = nx.Graph()
    for a, b in edges:
        w = float(dist_matrix[a, b])
        if w <= bt + eps:
            g.add_edge(int(a), int(b), weight=w)

    if g.has_edge(u, v):
        g.remove_edge(u, v)

    if (not g.has_node(u)) or (not g.has_node(v)) or (not nx.has_path(g, u, v)):
        g2 = nx.Graph()
        for a, b in edges:
            g2.add_edge(int(a), int(b), weight=float(dist_matrix[a, b]))
        if g2.has_edge(u, v):
            g2.remove_edge(u, v)
        try:
            path = nx.shortest_path(g2, u, v, weight="weight")
        except nx.NetworkXNoPath:
            return np.zeros(len(edges), dtype=float)
    else:
        t = nx.minimum_spanning_tree(g, weight="weight")
        path = nx.shortest_path(t, u, v, weight="weight")

    z0 = np.zeros(len(edges), dtype=float)

    for a, b in zip(path[:-1], path[1:]):
        e = (min(a, b), max(a, b))
        if e in edge_index:
            z0[edge_index[e]] += (+1.0 if a < b else -1.0)

    e_uv = (min(u0, v0), max(u0, v0))
    if e_uv in edge_index:
        z0[edge_index[e_uv]] += (-1.0 if u0 < v0 else +1.0)

    return z0


def build_cycle_problem(dist: np.ndarray, bar: dict[str, Any]):
    rips = gudhi.RipsComplex(distance_matrix=dist)
    st = rips.create_simplex_tree(max_dimension=2)

    cutoff = np.nextafter(float(bar["dt"]), -np.inf)
    edges, _, edge_idx, d_mat = build_boundary_at_cutoff(st, cutoff)
    z0 = construct_z0_from_birth_mst(
        bar["birth_edge"],
        edge_idx,
        edges,
        dist,
        float(bar["bt"]),
    )
    return edges, d_mat, z0
