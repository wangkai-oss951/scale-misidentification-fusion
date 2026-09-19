# -*- coding: utf-8 -*-
"""Print Table 1 / Table 2 from locked aggregates (paper numbers)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from config_util import add_config_arg, load_config  # noqa: E402
from stdio_utf8 import configure_stdio  # noqa: E402


def main() -> int:
    configure_stdio()
    p = add_config_arg(argparse.ArgumentParser())
    cfg = load_config(p.parse_args().config)
    agg = cfg["aggregates_dir"]
    extra = json.loads((agg / "paper_tables.json").read_text(encoding="utf-8"))
    frontier = json.loads((agg / "p0_near_maxr_frontier.json").read_text(encoding="utf-8"))
    ols = json.loads((agg / "p0_matched_ols_headroom_table.json").read_text(encoding="utf-8"))
    clin = json.loads((agg / "p0_clinical_mse_decomp.json").read_text(encoding="utf-8"))
    within = json.loads((agg / "p0_within_group_barrier.json").read_text(encoding="utf-8"))

    lin = [r for r in frontier["rows"] if r["arm"] == "mse_linear"]
    sig = [r for r in frontier["rows"] if r["arm"] == "mse_sigmoid"]
    med_ec = float(np.median([r["eps_cert"] for r in lin]))
    g005 = []
    for r in lin:
        for f in r["frontier"]:
            if abs(f["eps"] - 0.005) < 1e-12:
                g005.append(f["g_numerical"])
    med_g = float(np.median(g005))
    n_lin_cert = int(sum(r["eps_cert"] >= 0.01 for r in lin))
    n_sig_cert = int(sum(r["eps_cert"] >= 0.01 for r in sig))

    print("=== Table 1 (structural summary) ===")
    t1 = extra["table1"]
    for row in t1:
        print(" | ".join(row))
    print()
    print("Recomputed from locked frontier JSON:")
    print(f"  CEC2 linear certificate @ ε=.01: {n_lin_cert}/15  median ε_cert={med_ec:.3f}")
    print(f"  CEC2 reparam certificate @ ε=.01: {n_sig_cert}/15")
    print(f"  median numerical g(.005)={med_g:.3f}  Cor.2 floor={med_ec-0.005:.3f}")
    print(f"  linear G_aff {ols['linear_G_aff_range'][0]:.2f}–{ols['linear_G_aff_range'][1]:.2f}")
    print(f"  reparam G_aff {ols['reparam_G_aff_range'][0]:.2f}–{ols['reparam_G_aff_range'][1]:.2f}")
    s = clin["summary"]
    print(
        f"  clinical mean ΔRMSE={s['mean_delta_rmse']:.3f}  "
        f"G_q from |q-r|: raw q={s['raw_mean_q']:.3f}<r={s['raw_mean_r']:.3f} "
        f"→ q={s['ac_mean_q']:.3f}, r={s['ac_mean_r']:.3f}; "
        f"scale share={s['mean_scale_share_of_mse_reduction']:.3f}"
    )
    print()
    print("=== Table 2 (stress tests) ===")
    for row in extra["table2"]:
        print(" | ".join(row))
    print()
    print("Within-listener / within-scene Q<r fractions (15 systems):")
    print("  listener", within["frac_Q_lt_r_within_listener"],
          "scene", within["frac_Q_lt_r_within_scene"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
