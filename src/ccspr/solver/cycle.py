from __future__ import annotations

import cvxpy as cp
import numpy as np


def compute_cycle_representative(
    d_mat,
    z0: np.ndarray,
    solver: str = "l2",
    lambda_: float = 1e-6,
) -> np.ndarray:
    m, n = d_mat.shape
    if n == 0:
        return np.asarray(z0, dtype=float).copy()

    z0 = np.asarray(z0, dtype=np.float64)
    d64 = d_mat.astype(np.float64)

    y = cp.Variable(n)
    z = cp.Constant(z0) + cp.Constant(d64) @ y

    mode = str(solver).lower()
    if mode == "l2":
        prob = cp.Problem(cp.Minimize(cp.sum_squares(z)))
        solvers = [
            (
                cp.OSQP,
                dict(
                    eps_abs=1e-6,
                    eps_rel=1e-6,
                    max_iter=10000,
                    polish=False,
                    warm_start=True,
                    time_limit=5.0,
                ),
            ),
            (cp.ECOS, dict(abstol=1e-7, reltol=1e-7, feastol=1e-7, max_iters=2000)),
            (cp.SCS, dict(eps=1e-4, max_iters=5000)),
        ]
    elif mode == "elastic":
        t = cp.Variable(m)
        constraints = [t >= z, t >= -z]
        prob = cp.Problem(cp.Minimize(cp.sum(t) + float(lambda_) * cp.sum_squares(z)), constraints)
        solvers = [
            (
                cp.OSQP,
                dict(
                    eps_abs=1e-6,
                    eps_rel=1e-6,
                    max_iter=15000,
                    polish=False,
                    warm_start=True,
                    time_limit=7.0,
                ),
            ),
            (cp.ECOS, dict(abstol=1e-7, reltol=1e-7, feastol=1e-7, max_iters=3000)),
            (cp.SCS, dict(eps=1e-4, max_iters=6000)),
        ]
    else:
        raise ValueError("solver must be one of {'l2', 'elastic'}")

    z_val = None
    for solver_name, opts in solvers:
        try:
            prob.solve(solver=solver_name, verbose=False, **opts)
            if prob.status in ("optimal", "optimal_inaccurate") and z.value is not None:
                z_val = np.array(z.value).flatten()
                break
        except Exception:
            continue

    if z_val is None:
        return z0.copy()

    mx = np.max(np.abs(z_val)) + 1e-12
    z_val[np.abs(z_val) < 1e-6 * mx] = 0.0
    return z_val
