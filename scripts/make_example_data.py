# -*- coding: utf-8 -*-
"""Synthetic 5-head example so the optimizer/certificate path runs without Clarity data."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from affine_headroom import heldout_affine_headroom  # noqa: E402
from certificate import certificate_quantities, cor2_floor, numerical_frontier  # noqa: E402
from optimize_simplex import max_r_weights  # noqa: E402
from residualize import residualize_heads  # noqa: E402
from scale_metrics import mse_terms, pack, q_heads  # noqa: E402
from std_alignment import align_heads  # noqa: E402


def make_synthetic(n=400, k=5, seed=0):
    rng = np.random.default_rng(seed)
    y = rng.normal(50.0, 20.0, size=n)
    z = (y - y.mean()) / y.std()
    noise = rng.normal(0.0, 0.55, size=(k, n))
    # Shared direction, heterogeneous component scales (the paper's geometry).
    scales = np.array([2.0, 6.0, 3.0, 8.0, 1.2])
    P = scales[:, None] * (0.75 * z[None, :] + noise)
    P += rng.normal(0, 0.3, size=P.shape)
    listener = np.array([f"L{i % 8}" for i in range(n)])
    scene = np.array([f"S{i % 20}" for i in range(n)])
    dv = np.arange(0, n, 2)
    te = np.arange(1, n, 2)
    return {"P": P, "y": y, "listener": listener, "scene": scene, "dev": dv, "test": te}


def write_schema(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["column", "required", "notes"])
        w.writerow(["y", "yes", "intelligibility label on [0,100] (or any raw scale)"])
        w.writerow(["yhat_0..yhat_{k-1}", "yes", "component predictions, same rows as y"])
        w.writerow(["split", "yes", "train/dev/test"])
        w.writerow(["listener", "optional", "for within-listener residualization"])
        w.writerow(["scene", "optional", "for within-scene residualization"])
        w.writerow(["audio", "no", "not used; do not include participant audio"])
        w.writerow(["person_id", "no", "do not release identifiable IDs in public dumps"])


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=str(ROOT / "data" / "example"))
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    blob = make_synthetic()
    np.savez(
        out / "synthetic_heads.npz",
        P=blob["P"], y=blob["y"], listener=blob["listener"], scene=blob["scene"],
        dev=blob["dev"], test=blob["test"],
    )
    write_schema(out / "schema.csv")

    P, y, dv, te = blob["P"], blob["y"], blob["dev"], blob["test"]
    cert = certificate_quantities(P[:, dv], y[dv])
    num = numerical_frontier(P[:, dv], y[dv], cert["r_max"], 0.005, rng=0)
    w, _ = max_r_weights(P[:, dv], y[dv])
    st = pack(w @ P[:, te], y[te])
    terms = mse_terms(st)
    hr = heldout_affine_headroom(w @ P[:, dv], y[dv], w @ P[:, te], y[te])
    Pal = align_heads(P, y, dv, "std_only")
    cert_al = certificate_quantities(Pal[:, dv], y[dv])
    Pr, yr = residualize_heads(P, y, blob["listener"])
    cert_res = certificate_quantities(Pr[:, dv], yr[dv])

    print("example n,k", P.shape[1], P.shape[0])
    print("DEV r_max, Q, eps_cert", cert["r_max"], cert["Q"], cert["eps_cert"])
    print("g(.005) numerical", num["g"], "Cor.2 floor", cor2_floor(cert["eps_cert"], 0.005))
    print("TEST fused q,r,rmse", st["q"], st["r"], st["rmse"])
    print("MSE terms rel_err", terms["rel_err"])
    print("G_aff", hr["G_aff"])
    print("after std-only eps_cert", cert_al["eps_cert"])
    print("within-listener residual Q<r_max-0.01", cert_res["certificate_at_0.01"])
    print("WROTE", out / "synthetic_heads.npz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
