# -*- coding: utf-8 -*-
"""Within-group residualization: subtract group means from y and every yhat_j."""
from __future__ import annotations

import numpy as np


def residualize(vals, keys):
    vals = np.asarray(vals, np.float64)
    keys = np.asarray(keys)
    out = np.empty_like(vals)
    for k in np.unique(keys):
        m = keys == k
        out[m] = vals[m] - vals[m].mean()
    return out


def residualize_heads(P, y, keys):
    yr = residualize(y, keys)
    Pr = np.stack([residualize(P[j], keys) for j in range(P.shape[0])], 0)
    return Pr, yr
