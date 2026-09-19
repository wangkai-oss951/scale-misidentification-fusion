# -*- coding: utf-8 -*-
"""Scale / association metrics. Population moments ddof=0, as in the paper.

beta is the OLS prediction slope Cov(y, yhat) / Var(y), not a calibration slope.
"""
from __future__ import annotations

import numpy as np


def pearson(p, y) -> float:
    p = np.asarray(p, np.float64)
    y = np.asarray(y, np.float64)
    if p.std(ddof=0) < 1e-15 or y.std(ddof=0) < 1e-15:
        return 0.0
    return float(np.corrcoef(p, y)[0, 1])


def rmse(p, y) -> float:
    p = np.asarray(p, np.float64)
    y = np.asarray(y, np.float64)
    return float(np.sqrt(np.mean((p - y) ** 2)))


def beta_pred(p, y) -> float:
    y = np.asarray(y, np.float64)
    p = np.asarray(p, np.float64)
    sy2 = float(np.var(y, ddof=0))
    if sy2 < 1e-12:
        return float("nan")
    return float(np.cov(y, p, ddof=0)[0, 1] / sy2)


def pack(p, y) -> dict:
    p = np.asarray(p, np.float64)
    y = np.asarray(y, np.float64)
    sy = float(y.std(ddof=0))
    st = {
        "r": pearson(p, y),
        "rmse": rmse(p, y),
        "mae": float(np.mean(np.abs(p - y))),
        "beta": beta_pred(p, y),
        "q": float(p.std(ddof=0) / sy) if sy > 1e-12 else float("nan"),
        "mean_p": float(p.mean()),
        "mean_y": float(y.mean()),
        "dmu": float(p.mean() - y.mean()),
        "sp": float(p.std(ddof=0)),
        "sy": sy,
        "mse": float(np.mean((p - y) ** 2)),
    }
    return st


def mse_terms(st: dict) -> dict:
    """MSE = σ_y²(1-r²) + (Δμ)² + σ_y²(q-r)²  (ddof=0)."""
    sy2 = st["sy"] ** 2
    r, q, dmu = st["r"], st["q"], st["dmu"]
    assoc = sy2 * (1.0 - r * r)
    mean = dmu * dmu
    scale = sy2 * (q - r) ** 2
    recon = assoc + mean + scale
    rel = abs(recon - st["mse"]) / max(st["mse"], 1e-18)
    return {
        "mse_assoc": assoc,
        "mse_mean": mean,
        "mse_scale": scale,
        "mse_recon": recon,
        "mse_emp": st["mse"],
        "rel_err": rel,
    }


def q_heads(P, y) -> np.ndarray:
    P = np.asarray(P, np.float64)
    y = np.asarray(y, np.float64)
    sy = float(y.std(ddof=0))
    return P.std(axis=1, ddof=0) / max(sy, 1e-18)


def pi_from_w(w, sigma):
    w = np.asarray(w, np.float64)
    sigma = np.asarray(sigma, np.float64)
    num = w * sigma
    s = float(num.sum())
    if s <= 0:
        return np.ones_like(w) / w.size, 0.0
    return num / s, s


def w_from_pi(pi, sigma):
    pi = np.asarray(pi, np.float64)
    sigma = np.clip(np.asarray(sigma, np.float64), 1e-18, None)
    inv = pi / sigma
    z = float(inv.sum())
    return inv / z, 1.0 / z


def G_q(q_raw, r_raw, q_new, r_new) -> float:
    den = abs(q_raw - r_raw)
    if den < 1e-12:
        return float("nan")
    return 1.0 - abs(q_new - r_new) / den
