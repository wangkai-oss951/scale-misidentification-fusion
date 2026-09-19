# -*- coding: utf-8 -*-
"""One-command check that the released aggregates match the manuscript numbers.

This script does NOT retrain anything and does NOT re-derive the certificate from
raw prediction matrices (head predictions are not redistributed). It re-reads
`results/paper_locked/*.json` and asserts that each headline number in the paper
is the number those files actually contain, so a reader can confirm that the
released tables and the paper agree.

Usage:
    python scripts/verify_paper_numbers.py --config configs/cec2.yaml
Exit code is 0 only if every check passes.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from config_util import add_config_arg, load_config  # noqa: E402
from scale_metrics import G_q  # noqa: E402
from stdio_utf8 import configure_stdio  # noqa: E402

TOL = 1e-6
E2E_PAPER_DELTA_MAE = [1.94, 4.45, 1.46, 2.93, 2.18]


class Checker:
    def __init__(self) -> None:
        self.n = 0
        self.failed = 0

    def _emit(self, good: bool, name: str, detail: str) -> bool:
        self.n += 1
        if not good:
            self.failed += 1
        print(f"[{'PASS' if good else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        return bool(good)

    def eq(self, name: str, got, want, tol: float = TOL) -> bool:
        good = abs(float(got) - float(want)) <= tol
        return self._emit(good, name, f"got {float(got):.12g} expected {float(want):.12g}")

    def txt(self, name: str, got, want) -> bool:
        return self._emit(str(got) == str(want), name, f"got {got!r} expected {want!r}")

    def ge(self, name: str, got, bound) -> bool:
        good = float(got) >= float(bound) - TOL
        return self._emit(good, name, f"got {float(got):.12g} >= {float(bound):.12g}")

    def le(self, name: str, got, bound) -> bool:
        good = float(got) <= float(bound) + TOL
        return self._emit(good, name, f"got {float(got):.6g} <= {float(bound):.6g}")

    def ok(self, name: str, cond, detail: str = "") -> bool:
        return self._emit(bool(cond), name, detail)

    def report(self) -> int:
        print()
        print(f"{self.n - self.failed}/{self.n} checks passed")
        if self.failed:
            print(f"{self.failed} FAILED")
        return 1 if self.failed else 0


def load(agg: Path, name: str) -> dict:
    return json.loads((agg / name).read_text(encoding="utf-8"))


def main() -> int:
    configure_stdio()
    p = add_config_arg(argparse.ArgumentParser())
    cfg = load_config(p.parse_args().config)
    agg = Path(cfg["aggregates_dir"])
    c = Checker()
    print("Aggregates:", agg)
    print()

    frontier = load(agg, "p0_near_maxr_frontier.json")
    tables = load(agg, "paper_tables.json")
    ols = load(agg, "p0_matched_ols_headroom_table.json")
    e2e = load(agg, "p0_e2e_split.json")
    mech = load(agg, "p0_mechanism_locked_protocol.json")
    cec1 = load(agg, "p0_cec1_mechanism_replication.json")
    clin = load(agg, "p0_clinical_mse_decomp.json")
    within = load(agg, "p0_within_group_barrier.json")

    lin = [r for r in frontier["rows"] if r["arm"] == "mse_linear"]
    sig = [r for r in frontier["rows"] if r["arm"] == "mse_sigmoid"]

    print("-- CEC2 near-max-r certificate (Prop. 2 / Cor. 2) --")
    c.eq("CEC2 linear systems", len(lin), 15)
    c.eq("CEC2 reparam systems", len(sig), 15)
    eps_lin = [r["eps_cert"] for r in lin]
    med_ec = float(np.median(eps_lin))
    c.eq("median eps_cert (linear)", med_ec, 0.40435024073404935)
    c.eq("linear certificate at eps=.01", float(np.mean([e >= 0.01 for e in eps_lin])), 1.0)
    c.eq("reparam certificate at eps=.01", float(np.mean([r["eps_cert"] >= 0.01 for r in sig])), 1.0 / 15.0)
    g005 = [
        float(f["g_numerical"])
        for r in lin
        for f in r["frontier"]
        if abs(float(f["eps"]) - 0.005) < 1e-12
    ]
    med_g = float(np.median(g005))
    c.eq("median numerical g(.005)", med_g, 0.5646529106087451)
    c.eq("Cor.2 floor at eps=.005", med_ec - 0.005, 0.39935024073404934)
    c.ge("numerical g(.005) above the Cor.2 floor", med_g, med_ec - 0.005)
    per_sys = []
    for r in lin:
        gs = [float(f["g_numerical"]) for f in r["frontier"] if abs(float(f["eps"]) - 0.005) < 1e-12]
        if gs:
            per_sys.append(min(gs) - max(0.0, float(r["eps_cert"]) - 0.005))
    c.ok(
        "Cor.2 holds system-wise (15/15 at eps=.005)",
        len(per_sys) == 15 and all(v >= -1e-9 for v in per_sys),
        f"min slack {min(per_sys):.3g}",
    )
    freq = tables["certificate_freq"]
    c.eq("CEC1 certificate 3/3", freq["cec1"], 1.0)
    c.eq("CEC1 median eps_cert", freq["cec1_median_eps_cert"], 0.337)
    c.eq("WavLM certificate 3/3", freq["wavlm"], 1.0)
    c.eq("WavLM median eps_cert", freq["wavlm_median_eps_cert"], 0.259)

    print()
    print("-- Selection changes the slope (Fig. 1a/1b) --")
    maes = [float(e2e["splits"][str(i)]["delta_maxr_minus_minmae_obs"]["mae"]) for i in range(5)]
    for i, (got, want) in enumerate(zip(maes, E2E_PAPER_DELTA_MAE)):
        c.eq(f"split {i} e2e delta MAE rounds to paper value", round(got, 2), want, 1e-9)
    c.eq("split 0 delta MAE", maes[0], 1.9437325214695278)
    c.eq("locked protocol split 0 delta MAE", mech["paper_numbers_locked"]["e2e_split0_delta_mae"], round(maes[0], 2), 1e-9)
    direction = [
        e2e["splits"][str(i)]["delta_maxr_minus_minmae_obs"] for i in range(5)
    ]
    c.ok(
        "max-r minus min-MAE is directional 5/5 (r up, MAE up, beta down)",
        all(d["r"] > 0 and d["mae"] > 0 and d["beta"] < 0 for d in direction),
    )
    s1 = mech["sanity1_summary"]
    c.eq("CEC2 mean |dbeta| raw", s1["mean_abs_dbeta_unaligned"], 0.08652112727537493)
    c.eq("CEC2 mean |dbeta| std-only", s1["mean_abs_dbeta_std_only"], 0.003157497254845598)
    c.ok(
        "std-only alignment removes |dbeta|; mean-only does not",
        s1["mean_abs_dbeta_std_only"] < s1["mean_abs_dbeta_unaligned"]
        and s1["mean_abs_dbeta_mean_only"] >= s1["mean_abs_dbeta_unaligned"],
        f"mean-only {s1['mean_abs_dbeta_mean_only']:.4f}",
    )
    g5 = cec1["summary_grid5"]
    c.eq("CEC1 mean |dbeta| raw", g5["mean_abs_dbeta_raw"], 0.20228600705131314)
    c.eq("CEC1 mean |dbeta| std-only", g5["mean_abs_dbeta_std_only"], 0.00814723723617102)
    c.ok("CEC1 replication holds", g5["replication_holds"] and g5["std_only_kills_dbeta"] and g5["mean_only_does_not"])

    print()
    print("-- Affine headroom and MSE decomposition (Fig. 3) --")
    c.eq("linear G_aff range low", ols["linear_G_aff_range"][0], 8.90912241831615)
    c.eq("linear G_aff range high", ols["linear_G_aff_range"][1], 11.152351121944937)
    c.eq("reparam G_aff range low", ols["reparam_G_aff_range"][0], -0.02435024104601761)
    c.eq("reparam G_aff range high", ols["reparam_G_aff_range"][1], 0.34710440612820054)
    mse = tables["mse_share"]
    c.eq("CEC2 scale share of MSE reduction", mse["cec2_scale"], 0.9886)
    c.eq("clinical scale share of MSE reduction", mse["clinical_scale"], 0.8933)

    print()
    print("-- Clinical single head (TEST means; Q, r_max not applied) --")
    s = clin["summary"]
    c.eq("clinical raw q", s["raw_mean_q"], 0.10669175226021514)
    c.eq("clinical raw r", s["raw_mean_r"], 0.5273861307324736)
    c.ok("clinical raw q < r", s["raw_mean_q"] < s["raw_mean_r"], f"{s['raw_mean_q']:.4f} < {s['raw_mean_r']:.4f}")
    c.eq("clinical affine-coordinate q", s["ac_mean_q"], 0.4442298637177908)
    c.eq("clinical affine-coordinate r", s["ac_mean_r"], 0.5598851838217364)
    c.ok(
        "affine coordinate reduces |q-r|",
        abs(s["ac_mean_q"] - s["ac_mean_r"]) < abs(s["raw_mean_q"] - s["raw_mean_r"]),
        f"{abs(s['ac_mean_q'] - s['ac_mean_r']):.4f} < {abs(s['raw_mean_q'] - s['raw_mean_r']):.4f}",
    )
    c.eq("clinical mean delta RMSE", s["mean_delta_rmse"], -4.408686559675904)
    c.eq("clinical mean delta r", s["mean_delta_r"], 0.03249905308926271)
    gq = [G_q(r["raw"]["q"], r["raw"]["r"], r["affine_coordinate"]["q"], r["affine_coordinate"]["r"]) for r in clin["rows"]]
    c.eq("clinical seed-wise mean G_q", float(np.mean(gq)), 0.7268062569353408)
    c.le("clinical G_q below 1 (gap not closed)", float(np.mean(gq)), 1.0)

    print()
    print("-- Invariance and stress controls (Table 2) --")
    clo = tables["identifiability_closure"]
    c.le("T1 max relative MSE error", clo["t1_max_rel_mse_err"], 1e-14)
    c.le("T2 max |delta r|", clo["t2_max_abs_dr"], 1e-14)
    c.le("T3/T4 max |delta r|", clo["t3_t4_max_abs_dr"], 1e-14)
    uc = tables["unit_change_example"]
    c.eq("unit-change raw weights sum to 1", float(np.sum(uc["w"])), 1.0)
    c.eq("unit-change rescaled weights sum to 1", float(np.sum(uc["w_after"])), 1.0)
    c.ok("raw weights do change under affine units", uc["w"] != uc["w_after"], "pi and r unchanged to machine precision")
    c.eq("within-listener Q<r fraction", within["frac_Q_lt_r_within_listener"], 1.0)
    c.eq("within-scene Q<r fraction", within["frac_Q_lt_r_within_scene"], 1.0)
    c.eq("within-group systems", within["n"], 15)
    c.txt("archived census raw q<r", tables["census"]["explicitly_raw_uncalibrated_q_lt_r"], "61/74")
    c.ok("census is not an official leaderboard", tables["census"]["not_official_cpc2_leaderboard"] is True)
    census = load(agg, "p0_census_scale_landscape.json")
    comp = census["comparable"]
    vec = census["vector_G_aff"]
    c.txt(
        "census rank-agreement file pointer",
        tables["census"]["rank_agreement_file"],
        "p0_census_scale_landscape.json",
    )
    c.eq("census audited systems", census["n_rows"], 142)
    c.eq("census comparable systems", comp["n"], 141)
    c.eq("census systems with spread statistics", comp["n_with_q"], 77)
    c.eq("census Kendall tau (r vs RMSE)", comp["kendall_tau_r_rmse"], 0.018848804215646536)
    c.eq("census tau rounds to the paper value", round(comp["kendall_tau_r_rmse"], 3), 0.019, 1e-9)
    c.eq("census pairwise rank-inversion rate", comp["rank_inversion_rate"], 0.5094244021078232)
    c.eq("census inversion rounds to the paper value", round(comp["rank_inversion_rate"], 2), 0.51, 1e-9)
    c.eq("census rank pairs compared", comp["n_pairs"], 9868)
    c.eq("census affine-headroom subset size", vec["n"], 63)
    c.eq("census Spearman rho ((q-r)^2 vs G_aff)", vec["spearman_qr2_Gaff2"], 0.8141231817963609)
    c.eq("census rho rounds to the paper value", round(vec["spearman_qr2_Gaff2"], 2), 0.81, 1e-9)

    print()
    print("-- Guardrails carried by the locked protocol --")
    c.ok("not an official CPC2 evaluation", mech["not_official_cpc2_eval"] is True)
    c.ok("Job D is not a new algorithm", mech["job_d_is_not_a_new_algorithm"] is True)
    c.ok("CCC is not used as novelty", mech["do_not_use_ccc_as_novelty"] is True)
    c.ok("no variance-collapse framing", mech["do_not_use_variance_collapse_framing"] is True)
    c.ok("hits / n_words are not input features", mech["n_words"]["not_used_as_input_feature"] is True)

    print()
    return c.report()


if __name__ == "__main__":
    raise SystemExit(main())
