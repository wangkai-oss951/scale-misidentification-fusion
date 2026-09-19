# Evidence map: manuscript claim → released file

Every number below is reproduced by `scripts/verify_paper_numbers.py`, which reads
only the JSON files in `results/paper_locked/`. Run it first:

```bash
python scripts/verify_paper_numbers.py --config configs/cec2.yaml
```

To spot-check a single number by hand, open the file and search for the key path
(e.g. `sanity1_summary.mean_abs_dbeta_std_only`).

| Manuscript claim | File | Key path | Released value |
|---|---|---|---|
| CEC2 linear near-max-`r` certificate holds at `ε=.01` | `p0_near_maxr_frontier.json` | `rows[arm=mse_linear].eps_cert` | 15/15 |
| Median `ε_cert` over the 15 seed×split systems | `p0_near_maxr_frontier.json` | median of `rows[].eps_cert` | 0.404350 |
| Median numerical `g(.005)` | `p0_near_maxr_frontier.json` | median of `rows[].frontier[eps=.005].g_numerical` | 0.564653 |
| Cor. 2 floor at `ε=.005` | derived | `median ε_cert − .005` | 0.399350 |
| Certificate weakens under the matched reparameterization | `p0_near_maxr_frontier.json` | `rows[arm=mse_sigmoid].eps_cert` | 1/15 |
| CEC1 replication of the certificate | `p0_cec1_mechanism_replication.json` | `summary_grid5`, `splits` | 3/3, median `ε_cert` 0.337 |
| Cross-encoder (WavLM) replication | `paper_tables.json` | `certificate_freq.wavlm*` | 3/3, median `ε_cert` 0.259 |
| max-`r` vs min-MAE is directional in 5/5 splits | `p0_e2e_split.json` | `splits[i].delta_maxr_minus_minmae_obs` | ΔMAE 1.94 / 4.45 / 1.46 / 2.93 / 2.18 |
| DEV-only std alignment removes `Δβ`, mean-only does not | `p0_mechanism_locked_protocol.json` | `sanity1_summary` | `\|Δβ\|` 0.0865 → 0.00316 (mean-only 0.1377, no drop) |
| Same control on the CEC1 family | `p0_cec1_mechanism_replication.json` | `summary_grid5.mean_abs_dbeta_*` | 0.2023 → 0.00815 |
| Matched DEV-OLS affine headroom `G_aff` | `p0_matched_ols_headroom_table.json` | `linear_G_aff_range`, `reparam_G_aff_range` | 8.91–11.15 vs −0.02–0.35 |
| Scale dominates the MSE reduction (CEC2 / clinical) | `paper_tables.json` | `mse_share.cec2_scale`, `mse_share.clinical_scale` | 0.9886 / 0.8933 |
| Clinical single head: `q < r` before, closer after | `p0_clinical_mse_decomp.json` | `summary.raw_mean_q`, `ac_mean_q` | 0.1067 < 0.5274 → 0.4442 / 0.5599 |
| Clinical `G_q` closes 73% of the `\|q−r\|` gap | `p0_clinical_mse_decomp.json` | per-seed `G_q` over `rows[]` | mean 0.7268 |
| Positivity/units do not move `π` or `r` (machine precision) | `paper_tables.json` | `identifiability_closure`, `unit_change_example` | ≤2.1e-15, ≤7.8e-16 |
| Group-mean removal does not create the barrier | `p0_within_group_barrier.json` | `frac_Q_lt_r_within_listener`, `..._scene` | 1.0 / 1.0 over 15 systems |
| Archived census: explicitly raw/q<r | `p0_census_scale_landscape.json` | `by_calibration_status.explicitly_raw_uncalibrated` | 61/74 (82.4%) |
| Same file, unfiltered denominator | `p0_census_scale_landscape.json` | `comparable.frac_q_lt_r`, `comparable.n_with_q` | 61/77 (79.2%) — the 74/77 split is the calibration-status filter |
| Pearson vs raw-RMSE rankings nearly unrelated | `p0_census_scale_landscape.json` | `comparable.kendall_tau_r_rmse`, `.rank_inversion_rate`, `.n_pairs` | τ=0.01885, inversion 0.50942, 9868 pairs |
| Census coverage of the archived audit | `p0_census_scale_landscape.json` | `n_rows`, `n_comparable`, `comparable.n_with_q` | 142 audited / 141 comparable / 77 with spread statistics |
| `(q−r)²` tracks held-out affine headroom | `p0_census_scale_landscape.json` | `vector_G_aff.spearman_qr2_Gaff2`, `.n` | Spearman ρ=0.81412 over n=63 systems |

## What these files are not

- They are **aggregates of the CPU analysis**, not raw predictions. Head prediction
  matrices, waveforms, listener identifiers, and clinical vectors are not released.
- The near-max-`r` frontier is a **numerical** search (grid + Dirichlet + SLSQP), so
  `g(ε)` is an upper bound on what is achievable at that `ε`; the theorem direction is
  the certificate `Q < r_max − ε`. `g(ε)` is never used as a proof of the bound.
- The clinical corpus is a **single-head** person-holdout. `Q` and `r_max` are
  five-head feasible-set objects and are deliberately **not** applied to it.
- `p0_census_scale_landscape.json` is the **aggregate** census of archived Clarity-ecosystem
  systems: counts and rank statistics only. It is not an official CPC2 leaderboard census,
  and no prediction vector is released. It does carry one archival descriptor
  (`other_descriptive`) naming the single non-comparable system, which is a system
  published in the literature, not clinical data.
- `p0_mechanism_locked_protocol.json` also records the protocol guardrails (not an
  official CPC2 evaluation; Job D reparameterization is a control, not a new
  algorithm; CCC and variance-collapse framings are excluded; `hits`/`n_words` are
  likelihood labels only).
