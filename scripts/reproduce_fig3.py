# -*- coding: utf-8 -*-
"""Fig. 3: matched DEV-OLS affine headroom and MSE-decomposition shares."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from config_util import add_config_arg, load_config  # noqa: E402


def _plot(cfg: dict) -> None:
    agg = cfg["aggregates_dir"]
    ols = json.loads((agg / "p0_matched_ols_headroom_table.json").read_text(encoding="utf-8"))
    extra = json.loads((agg / "paper_tables.json").read_text(encoding="utf-8"))
    rows = sorted(ols["rows"], key=lambda r: r["split"])

    fig, axes = plt.subplots(2, 1, figsize=(3.45, 4.55))
    fig.subplots_adjust(left=0.20, right=0.98, bottom=0.08, top=0.96, hspace=0.38)

    ax = axes[0]
    x = np.arange(5)
    w = 0.18
    series = [
        ("linear raw", [r["linear_raw_rmse"] for r in rows], "#1f4e79"),
        ("linear+OLS", [r["linear_plus_dev_ols_rmse"] for r in rows], "#7a93b0"),
        ("reparam raw", [r["reparam_raw_rmse"] for r in rows], "#c45c26"),
        ("reparam+OLS", [r["reparam_plus_dev_ols_rmse"] for r in rows], "#e0a078"),
    ]
    for i, (lab, vals, col) in enumerate(series):
        ax.bar(x + (i - 1.5) * w, vals, width=w, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([str(i) for i in range(5)])
    ax.set_xlabel("CEC2 split")
    ax.set_ylabel("TEST RMSE")
    ax.set_ylim(20, 42)
    ax.set_title("(a)", pad=2)
    ax.legend(frameon=False, ncol=2, loc="upper right", columnspacing=0.6, handlelength=1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    mse = extra["mse_share"]
    scale = np.array([mse["cec2_scale"], mse["clinical_scale"]])
    assoc = np.array([mse["cec2_assoc"], mse["clinical_assoc"]])
    mean_t = np.array([mse["cec2_mean"], mse["clinical_mean"]])
    xpos = np.arange(2)
    bw = 0.28
    ax.bar(xpos - bw, scale, width=bw, color="#1f4e79", label="scale")
    ax.bar(xpos, assoc, width=bw, color="#2e7d4f", label="association")
    ax.bar(xpos + bw, mean_t, width=bw, color="#c45c26", label="mean")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels(["CEC2", "Clinical"])
    ax.set_ylabel("share of MSE reduction")
    ax.set_ylim(-0.22, 1.18)
    ax.set_title("(b)", pad=2)
    for i, v in enumerate(scale):
        ax.text(i - bw, v + 0.03, f"{100 * v:.1f}%", ha="center", va="bottom", fontsize=7, color="#1f4e79")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = cfg["output_dir"] / "fig3_mechanism_transfer.png"
    fig.savefig(out, dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    print("WROTE", out)
    print("linear G_aff range", ols["linear_G_aff_range"])
    print("reparam G_aff range", ols["reparam_G_aff_range"])
    print("MSE shares CEC2/clinical scale", mse["cec2_scale"], mse["clinical_scale"])


def main() -> int:
    p = add_config_arg(argparse.ArgumentParser())
    _plot(load_config(p.parse_args().config))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
