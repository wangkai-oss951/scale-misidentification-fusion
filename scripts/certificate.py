# -*- coding: utf-8 -*-
"""Prop. 2 / Cor. 2: ε_cert = r_max − Q and numerical near-max-r gap g(ε).

`numerical_frontier` is the paper's CPU search (grid + Dirichlet + SLSQP).
It is not an exact global optimum; the certificate itself is Q < r_max − ε.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from optimize_simplex import default_inits, max_r_weights
from scale_metrics import pack, q_heads


def r_q_of(P, w, y):
    p = np.asarray(w, np.float64) @ np.asarray(P, np.float64)
    st = pack(p, y)
    return st["r"], st["q"], st


def grid_simplex(k=5, step=5):
    import itertools

    nbin = 100 // step
    out = []
    for head in itertools.product(range(nbin + 1), repeat=k - 1):
        s = sum(head)
        if s > nbin:
            continue
        parts = head + (nbin - s,)
        out.append(np.array([p * step / 100.0 for p in parts], np.float64))
    return out


def sample_ws(k=5, n_dir=400, step=5, rng=0):
    rng = np.random.default_rng(rng)
    ws = grid_simplex(k, step)
    ws.extend(default_inits(k))
    for _ in range(n_dir):
        ws.append(rng.dirichlet(np.ones(k)))
    return ws


def certificate_quantities(P, y):
    """DEV component matrix P (k, n_dev) and labels y (n_dev,)."""
    w, r_max = max_r_weights(P, y)
    Q = float(np.max(q_heads(P, y)))
    eps_cert = float(r_max - Q)
    return {
        "w": w,
        "r_max": float(r_max),
        "Q": Q,
        "eps_cert": eps_cert,
        "certificate_at_0.01": bool(Q < r_max - 0.01),
    }


def numerical_frontier(P, y, r_max, eps, rng=0):
    """min |q-r| over sampled/SLSQP weights with r ≥ r_max − ε. Numerical, not exact."""
    P = np.asarray(P, np.float64)
    y = np.asarray(y, np.float64)
    k = P.shape[0]
    need = r_max - float(eps)
    cand = []
    for w in sample_ws(k, rng=rng):
        w = np.clip(np.asarray(w, np.float64), 0, 1)
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        r, q, _ = r_q_of(P, w, y)
        if r + 1e-12 >= need:
            cand.append((w, r, q, abs(q - r)))

    def rq(w):
        w = np.clip(w, 0, 1)
        if w.sum() <= 0:
            w = np.ones(k) / k
        else:
            w = w / w.sum()
        return r_q_of(P, w, y)

    cons = [
        {"type": "eq", "fun": lambda w: np.sum(w) - 1.0},
        {"type": "ineq", "fun": lambda w: rq(w)[0] - need},
    ]
    bounds = [(0.0, 1.0)] * k
    starts = [c[0] for c in cand[:8]] if cand else default_inits(k)[:8]
    starts.append(np.ones(k) / k)

    def min_gap(w):
        r, q, _ = rq(w)
        return abs(q - r)

    def neg_q(w):
        return -rq(w)[1]

    for fun in (min_gap, neg_q):
        for w0 in starts:
            try:
                res = minimize(
                    fun, w0, method="SLSQP", bounds=bounds, constraints=cons,
                    options={"maxiter": 200, "ftol": 1e-12},
                )
                w = np.clip(res.x, 0, 1)
                if w.sum() <= 0:
                    continue
                w = w / w.sum()
                r, q, _ = rq(w)
                if r + 1e-10 >= need:
                    cand.append((w, r, q, abs(q - r)))
            except Exception:
                continue
    if not cand:
        return {"g": float("nan"), "q_max": float("nan"), "n_feas": 0}
    g = min(c[3] for c in cand)
    q_max = max(c[2] for c in cand)
    return {"g": float(g), "q_max": float(q_max), "n_feas": int(len(cand))}


def cor2_floor(eps_cert, eps):
    """g(ε) ≥ ε_cert − ε on the near-max-r set (Corollary 2)."""
    return max(0.0, float(eps_cert) - float(eps))
