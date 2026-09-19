# Scale misidentification in correlation-based intelligibility fusion

Reproduction code for the ICASSP analysis:

> *Scale Misidentification in Correlation-Based Intelligibility Fusion*

This repository redistributes **analysis code**, **synthetic examples**, **schema**, and **locked aggregate results**. It does **not** redistribute Clarity/CPC2 audio, listener identifiers, clinical recordings, or clinical prediction vectors.

The contribution implemented here is diagnostic, not a new architecture:

1. Component scale is **non-identifiable** under correlation-selected convex fusion (`w ↔ π` on the simplex).
2. A DEV margin `ε_cert = r_max − Q` is a **global near-optimal unattainability certificate** and a quantitative floor `g(ε) ≥ ε_cert − ε` (Prop. 2 / Cor. 2).
3. Raw simplex weights are **coordinate-dependent** and are not scale-free importance.

Pearson invariance itself is treated as a known fact, not a contribution.

## Paper → script map

| Paper | Script / module |
|---|---|
| Fig. 1a selection-dependent slope | `scripts/reproduce_fig1.py` |
| Fig. 1b DEV-only std alignment | `scripts/std_alignment.py` + Fig. 1 |
| Fig. 1c spread ceiling `Q`, `q_fusion` | `scripts/reproduce_fig1.py` |
| Simplex max-`r` / min-MAE | `scripts/optimize_simplex.py` (paper SLSQP, not a rewrite) |
| Prop. 2 / Cor. 2 `r_max`, `Q`, `ε_cert`, `g(ε)` | `scripts/certificate.py` |
| Fig. 2 near-max-`r` frontier | `scripts/reproduce_fig2.py` |
| Fig. 3a DEV-OLS affine headroom | `scripts/affine_headroom.py` + `reproduce_fig3.py` |
| Fig. 3b MSE decomposition shares | `scripts/scale_metrics.py` `mse_terms` + `reproduce_fig3.py` |
| Table 1 / Table 2 | `scripts/reproduce_tables.py` |
| Within-listener / within-scene residualization | `scripts/residualize.py` |

## Quick start

```bash
python -m pip install -r requirements.txt
python scripts/make_example_data.py
python scripts/reproduce_fig2.py --config configs/cec2.yaml
python scripts/reproduce_fig3.py --config configs/cec2.yaml
python scripts/reproduce_tables.py --config configs/cec2.yaml
python scripts/reproduce_fig1.py --config configs/cec2.yaml
```

`configs/cec2.yaml` defaults to **paper aggregates** in `results/paper_locked/`. Those JSON files are the locked CPU outputs used in the manuscript (CEC2 frontier, OLS headroom, Fig. 1 sources, clinical *seed-level* moments only).

The three `reproduce_fig2/3` + `reproduce_tables` commands reprint the core numbers:

- CEC2 linear certificate **15/15**, median `ε_cert ≈ 0.404`, median `g(.005) = 0.565`
- reparam certificate **1/15**
- `G_aff` **8.91–11.15 → −0.02–0.35**
- CEC2 / clinical scale shares **98.9% / 89.3%** (mean of seed-wise MSE-reduction shares)

## What is not in this repo

- HuBERT / WavLM training or checkpoints
- CEC2 / CEC1 waveforms or official CPC2 labels
- Clinical participant-level audio, names, or `full_pred.npz`
- hits / n_words as input features (labels only, in the original pipeline)

Clinical evidence in the paper is a **single-head** `raw768_time` person-holdout. The five-head quantities `Q` and `r_max` are fusion feasible-set objects and are **not** applied to that corpus. Released clinical numbers are aggregate `r,q,RMSE,β` and MSE-term shares.

To recompute fusion diagnostics from your own five-head matrices, store `P` as shape `(k, n)` with a DEV index and call:

```python
from certificate import certificate_quantities, numerical_frontier
cert = certificate_quantities(P[:, dev], y[dev])
g = numerical_frontier(P[:, dev], y[dev], cert["r_max"], eps=0.005)
```

Moments use **ddof=0**. `β = Cov(y,ŷ)/Var(y)` is an OLS prediction slope, not a calibration slope. Affine-coordinate training in the paper is a **matched training-coordinate control**, not a proposed new method.

## Data schema (public)

See `data/example/schema.csv`. `scripts/make_example_data.py` writes a tiny synthetic five-head set that exhibits `Q < r_max` so the certificate path is executable without Challenge data.

## License

MIT for the code in this repository. Clarity/CPC2 and clinical corpora remain under their original terms; obtain those data separately if you need to retrain heads.

## 中文说明

本仓库只开源分析代码、合成示例、字段说明和论文锁定的**汇总 JSON**。不公开临床原始数据、音频、可识别 ID 或预测向量。CEC2 图和表默认从 `results/paper_locked/` 复现；simplex 优化器是论文里实际使用的 SLSQP 实现，而非近似版。
