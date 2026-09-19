# -*- coding: utf-8 -*-
"""Continuous simplex fusion used in the paper (CEC2 max-r / min-MAE).

Copied from the analysis code (`p0_cpu_batch1.continuous_simplex`): SLSQP on the
probability simplex, multiple initializations. Not a new fusion algorithm.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

# Locked in the paper analysis (p0_cpu_batch1.py).
DIR_SEED = 20260911
N_DIRICHLET = 8


def pearson(pred, y) -> float:
    pred = np.asarray(pred, np.float64)
    y = np.asarray(y, np.float64)
    sp, sy = float(pred.std()), float(y.std())
    if sp < 1e-15 or sy < 1e-15:
        return 0.0
    pc, yc = pred - pred.mean(), y - y.mean()
    return float(np.dot(pc, yc) / (np.linalg.norm(pc) * np.linalg.norm(yc)))


def mae(pred, y) -> float:
    return float(np.mean(np.abs(np.asarray(pred, np.float64) - np.asarray(y, np.float64))))


def default_inits(k: int = 5, extra=None):
    rng = np.random.default_rng(DIR_SEED)
    inits = [np.ones(k) / k]
    for i in range(k):
        e = np.zeros(k)
        e[i] = 1.0
        inits.append(e)
    for _ in range(N_DIRICHLET):
        inits.append(rng.dirichlet(np.ones(k)))
    if extra:
        inits.extend(extra)
    return inits


def continuous_simplex(P, y, objective: str, inits=None):
    """Maximize Pearson (`max_r`) or minimize MAE (`min_mae`) over w in the simplex.

    P: (k, n) component predictions on the same rows as y (typically DEV).
    Returns (w_list, score) with w on the simplex.
    """
    P = np.asarray(P, np.float64)
    y = np.asarray(y, np.float64)
    k = P.shape[0]
    if inits is None:
        inits = default_inits(k)

    def fun(w):
        pred = w @ P
        if objective == "max_r":
            return -pearson(pred, y)
        if objective == "min_mae":
            return mae(pred, y)
        raise ValueError(objective)

    cons = {"type": "eq", "fun": lambda w: np.sum(w) - 1.0}
    bounds = [(0.0, 1.0)] * k
    best_w, best_v = None, None
    for w0 in inits:
        w0 = np.asarray(w0, np.float64)
        w0 = np.clip(w0, 0.0, 1.0)
        if w0.sum() <= 0:
            w0 = np.ones(k) / k
        else:
            w0 = w0 / w0.sum()
        res = minimize(
            fun, w0, method="SLSQP", bounds=bounds, constraints=cons,
            options={"maxiter": 400, "ftol": 1e-12},
        )
        w = np.clip(res.x, 0.0, 1.0)
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        v = float(fun(w))
        if best_v is None or v < best_v:
            best_v, best_w = v, w
    if best_w is None:
        best_w = np.ones(k) / k
        best_v = float(fun(best_w))
    score = -best_v if objective == "max_r" else best_v
    return best_w.tolist(), float(score)


def max_r_weights(P, y):
    w, r = continuous_simplex(P, y, "max_r", default_inits(P.shape[0]))
    return np.asarray(w, np.float64), float(r)
