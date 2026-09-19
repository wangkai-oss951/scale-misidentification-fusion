# -*- coding: utf-8 -*-
"""DEV-only moment alignment controls (mechanism, not a submitted system).

Paper formula (std-only):
  ytilde_j = mu_j^DEV + (sigma_y^DEV / sigma_hat_j^DEV) * (hat_j - mu_hat_j^DEV)
Moments use ddof=0. Positive scaling does not change per-head Pearson.
"""
from __future__ import annotations

import numpy as np


def align_mean_only(hat, y, dv):
    hat = np.asarray(hat, np.float64)
    y = np.asarray(y, np.float64)
    return hat - float(hat[dv].mean()) + float(y[dv].mean())


def align_std_only(hat, y, dv):
    hat = np.asarray(hat, np.float64)
    y = np.asarray(y, np.float64)
    sig_y = float(y[dv].std(ddof=0))
    sig_h = float(hat[dv].std(ddof=0))
    mu_h = float(hat[dv].mean())
    if sig_h < 1e-12:
        return np.full_like(hat, mu_h)
    return mu_h + (sig_y / sig_h) * (hat - mu_h)


def align_full(hat, y, dv):
    hat = np.asarray(hat, np.float64)
    y = np.asarray(y, np.float64)
    mu_y = float(y[dv].mean())
    sig_y = float(y[dv].std(ddof=0))
    mu_h = float(hat[dv].mean())
    sig_h = float(hat[dv].std(ddof=0))
    if sig_h < 1e-12:
        return np.full_like(hat, mu_y)
    return mu_y + (sig_y / sig_h) * (hat - mu_h)


def align_heads(P, y, dv, mode: str = "std_only"):
    """P (k, n). Apply the same DEV-only map independently to each row."""
    fn = {"mean_only": align_mean_only, "std_only": align_std_only, "full": align_full}[mode]
    return np.stack([fn(P[j], y, dv) for j in range(P.shape[0])], 0)
