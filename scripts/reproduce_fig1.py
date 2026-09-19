# -*- coding: utf-8 -*-
"""Fig. 1: selection-dependent slope, std-only alignment, spread ceiling."""
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
from stdio_utf8 import configure_stdio  # noqa: E402


def _plot(cfg: dict) -> Path:
    agg = cfg["aggregates_dir"]
    e2e = json.loads((agg / "p0_e2e_split.json").read_text(encoding="utf-8"))
    mech = json.loads((agg / "p0_mechanism_locked_protocol.json").read_text(encoding="utf-8"))
    cec1 = json.loads((agg / "p0_cec1_mechanism_replication.json").read_text(encoding="utf-8"))
    spread = json.loads((agg / "p0_spread_bound_table.json").read_text(encoding="utf-8"))

    rngs = [0, 1, 2, 3, 4]
    dbeta = [e2e["splits"][str(i)]["delta_maxr_minus_minmae_obs"]["beta"] for i in rngs]
    cec2_raw = [
        abs(mech["splits"][str(i)]["sanity1"]["unaligned_from_p0_e2e_split"]["delta"]["beta"])
        for i in rngs
    ]
    cec2_std = [
        abs(mech["splits"][str(i)]["sanity1"]["std_only"]["delta"]["beta"])
        for i in rngs
    ]
    g = cec1["summary_grid5"]
    cec1_raw_mean = g["mean_abs_dbeta_raw"]
    cec1_std_mean = g["mean_abs_dbeta_std_only"]
    cec1_keys = sorted(cec1["splits"].keys(), key=int)
    cec1_raw = [
        abs(cec1["splits"][k]["alignments"]["raw"]["grid5"]["delta_test"]["beta"])
        for k in cec1_keys
    ]
    cec1_std = [
        abs(cec1["splits"][k]["alignments"]["std_only"]["grid5"]["delta_test"]["beta"])
        for k in cec1_keys
    ]

    fig, axes = plt.subplots(1, 3, figsize=(7.16, 2.42))
    ax = axes[0]
    ax.axhline(0.0, color="#888", lw=0.6)
    ax.plot(rngs, dbeta, "o-", color="#1f4e79", ms=5, lw=1.2)
    ax.set_xticks(rngs)
    ax.set_xlabel("CEC2 split")
    ax.set_ylabel(r"$\Delta\beta$ (max-$r$ $-$ min-MAE)")
    ax.set_title("(a) Selection-dependent slope")

    ax = axes[1]
    x = np.arange(2)
    w = 0.32
    ax.bar(x - w / 2, [np.mean(cec2_raw), cec1_raw_mean], w, color="#c0392b", label="raw")
    ax.bar(x + w / 2, [np.mean(cec2_std), cec1_std_mean], w, color="#1f4e79", label="std-only")
    ax.scatter(np.full(5, x[0] - w / 2), cec2_raw, s=10, c="k", zorder=4)
    ax.scatter(np.full(5, x[0] + w / 2), cec2_std, s=10, c="k", zorder=4)
    ax.scatter(np.full(len(cec1_raw), x[1] - w / 2), cec1_raw, s=10, c="k", zorder=4)
    ax.scatter(np.full(len(cec1_std), x[1] + w / 2), cec1_std, s=10, c="k", zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(["CEC2 (5 splits)", "CEC1 (3 splits)"])
    ax.set_ylabel(r"mean $|\Delta\beta|$")
    ax.set_title("(b) Std-only alignment")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[2]
    heads_short = ["time", "fft", "L8", "wide", "L6"]
    q_mat, q_fus, q_ceil = [], [], []
    for i in rngs:
        rec = spread["splits"][str(i)]
        q_mat.append([h["q"] for h in rec["heads"]])
        q_fus.append(rec["fusions"]["max_r"]["q_fusion"])
        q_ceil.append(rec["max_component_q"])
    q_mat = np.array(q_mat, dtype=float)
    colors = ["#4c78a8", "#f58518", "#54a24b", "#e45756", "#b279a2"]
    jitter = np.linspace(-0.18, 0.18, q_mat.shape[1])
    for j, name in enumerate(heads_short):
        ax.scatter(np.array(rngs) + jitter[j], q_mat[:, j], s=14, c=colors[j],
                   label=name, zorder=3)
    ax.plot(rngs, q_ceil, "s--", color="#333", ms=4, lw=0.9, label=r"$Q$")
    ax.plot(rngs, q_fus, "D-", color="#c0392b", ms=4.5, lw=1.1, label=r"$q_{\mathrm{fus}}$")
    ax.axhline(1.0, color="#2c7a4b", lw=1.0, ls=":")
    ax.text(4.32, 1.02, r"$q=1$", color="#2c7a4b", fontsize=7.2, ha="right", va="bottom")
    ax.set_ylim(0, 1.18)
    ax.set_xlim(-0.35, 4.35)
    ax.set_xticks(rngs)
    ax.set_xlabel("CEC2 split")
    ax.set_ylabel(r"$q=\sigma_{\hat y}/\sigma_y$")
    ax.set_title("(c) Spread ceiling")
    ax.legend(loc="upper left", ncol=2, fontsize=7.0, frameon=True, fancybox=False,
              edgecolor="#cccccc", framealpha=0.92)
    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.tight_layout(w_pad=0.9)
    out = cfg["output_dir"] / "fig1_mechanism.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("WROTE", out)
    print("Fig.1a dbeta (max-r minus min-MAE):", [round(v, 4) for v in dbeta])
    print("Fig.1b CEC2 mean |dbeta| raw/std-only:", float(np.mean(cec2_raw)), float(np.mean(cec2_std)))
    return out


def main() -> int:
    configure_stdio()
    p = add_config_arg(argparse.ArgumentParser())
    cfg = load_config(p.parse_args().config)
    _plot(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
