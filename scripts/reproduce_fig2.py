# -*- coding: utf-8 -*-
"""Fig. 2: near-max-r frontier g(ε) vs Cor. 2 floor, Q vs r_max, certificate rates."""
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

EPS_LIST = [0.0, 0.0025, 0.005, 0.01, 0.02]


def collect(rows, arm: str):
    sub = [r for r in rows if r["arm"] == arm]
    eps_cert = [r["eps_cert"] for r in sub]
    rmax = [r["r_max_dev"] for r in sub]
    Q = [r["Q_dev"] for r in sub]
    g = {e: [] for e in EPS_LIST}
    for r in sub:
        for f in r["frontier"]:
            g[float(f["eps"])].append(f["g_numerical"])
    return eps_cert, rmax, Q, g


def _plot(cfg: dict) -> None:
    agg = cfg["aggregates_dir"]
    d = json.loads((agg / "p0_near_maxr_frontier.json").read_text(encoding="utf-8"))
    extra = json.loads((agg / "paper_tables.json").read_text(encoding="utf-8"))
    rows = d["rows"]
    lin_ec, lin_r, lin_Q, lin_g = collect(rows, "mse_linear")
    sig_ec, sig_r, sig_Q, sig_g = collect(rows, "mse_sigmoid")

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.15))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.22, top=0.88, wspace=0.38)

    ax = axes[0]
    for gdict, color, label in (
        (lin_g, "#1f4e79", "linear"),
        (sig_g, "#c45c26", "reparam"),
    ):
        xs, med, lo, hi = [], [], [], []
        for e in EPS_LIST:
            vals = np.asarray(gdict[e], float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                continue
            xs.append(e)
            med.append(float(np.median(vals)))
            lo.append(float(np.percentile(vals, 25)))
            hi.append(float(np.percentile(vals, 75)))
        ax.plot(xs, med, "-o", color=color, label=label, ms=3.5, lw=1.4)
        ax.fill_between(xs, lo, hi, color=color, alpha=0.16, linewidth=0)
    med_ec = float(np.median(lin_ec))
    xs_b = np.array(EPS_LIST, float)
    ax.plot(xs_b, np.maximum(0.0, med_ec - xs_b), "--", color="#1f4e79", lw=1.1,
            label=r"bound $\epsilon_{\mathrm{cert}}-\epsilon$")
    ax.axhline(0, color="0.55", lw=0.6)
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$g(\epsilon)$")
    ax.set_title("(a)", pad=4)
    ax.legend(frameon=False, loc="upper right", fontsize=6)
    ax.set_xlim(-0.001, 0.021)

    ax = axes[1]
    ax.scatter(lin_r, lin_Q, c="#1f4e79", s=22, label="linear", zorder=3)
    ax.scatter(sig_r, sig_Q, c="#c45c26", s=22, marker="s", label="reparam", zorder=3)
    lo = min(min(lin_Q), min(sig_Q), min(lin_r), min(sig_r)) - 0.05
    hi = max(max(lin_Q), max(sig_Q), max(lin_r), max(sig_r)) + 0.05
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.7)
    ax.text(0.04, 0.96, rf"med $\epsilon_{{\mathrm{{cert}}}}={med_ec:.2f}$",
            transform=ax.transAxes, va="top", ha="left", color="#1f4e79", fontsize=7)
    ax.set_xlabel(r"$r_{\max}$ (DEV)")
    ax.set_ylabel(r"$Q$ (DEV)")
    ax.set_title("(b)", pad=4)
    ax.legend(frameon=False, loc="lower right", fontsize=6.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    ax = axes[2]
    cert = extra["certificate_freq"]
    labs = ["CEC2\nlin.", "CEC2\nrep.", "CEC1", "WavLM"]
    vals = [cert["cec2_linear"], cert["cec2_reparam"], cert["cec1"], cert["wavlm"]]
    nlabs = [cert["cec2_linear_lab"], cert["cec2_reparam_lab"], cert["cec1_lab"], cert["wavlm_lab"]]
    cols = ["#1f4e79", "#7a93b0", "#2e7d4f", "#5c4d7a"]
    bars = ax.bar(np.arange(len(labs)), vals, color=cols, width=0.68)
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, fontsize=6.5)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("fraction")
    ax.set_title("(c)", pad=4)
    for b, nlab in zip(bars, nlabs):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.03, nlab,
                ha="center", va="bottom", fontsize=6.5)
    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)

    out = cfg["output_dir"] / "fig2_near_maxr_certificate.png"
    fig.savefig(out, dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)

    g005 = [v for v in lin_g[0.005] if v == v]
    print("WROTE", out)
    print("CEC2 linear n=", len(lin_ec))
    print("median eps_cert=", med_ec)
    print("median g(.005)=", float(np.median(g005)))
    print("Cor.2 floor at .005=", med_ec - 0.005)
    print("frac eps_cert>=.01=", float(np.mean([e >= 0.01 for e in lin_ec])))
    print("reparam frac eps_cert>=.01=", float(np.mean([e >= 0.01 for e in sig_ec])))


def main() -> int:
    p = add_config_arg(argparse.ArgumentParser())
    _plot(load_config(p.parse_args().config))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
