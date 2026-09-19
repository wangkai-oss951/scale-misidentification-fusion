# -*- coding: utf-8 -*-
"""DEV-only unclipped OLS affine map, frozen on TEST. G_aff = RMSE_raw − RMSE_ols."""
from __future__ import annotations

import numpy as np

from scale_metrics import rmse


def fit_ols(p, y):
    p = np.asarray(p, np.float64)
    y = np.asarray(y, np.float64)
    vx = float(np.var(p, ddof=0))
    if vx < 1e-18:
        return 0.0, float(y.mean())
    a = float(np.cov(p, y, ddof=0)[0, 1] / vx)
    b = float(y.mean() - a * p.mean())
    return a, b


def apply_ols(p, a, b):
    return a * np.asarray(p, np.float64) + b


def heldout_affine_headroom(p_dev, y_dev, p_test, y_test) -> dict:
    """Fit unclipped OLS on DEV, evaluate raw vs mapped RMSE on TEST."""
    a, b = fit_ols(p_dev, y_dev)
    raw = rmse(p_test, y_test)
    mapped = rmse(apply_ols(p_test, a, b), y_test)
    return {
        "a": a,
        "b": b,
        "rmse_raw": raw,
        "rmse_ols": mapped,
        "G_aff": float(raw - mapped),
    }
