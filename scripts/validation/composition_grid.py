"""Composition grid: does the capstone's single-scenario finding generalize?

METHODS_NOTE_ROBUSTNESS.md Sec. 7 flagged that the sequential-composition
overcorrection (and the joint estimator's advantage) was "tested for exactly
one scenario". This script runs the full 5-arm comparison over a grid of
detection strengths x allocation regimes, all against the POPULATION truth:

  arms:  naive | detection-only | sparse-only | composed-sequential | joint

Writes figures/composition_grid.png + figures/composition_grid_results.csv.
Synthetic data only. Run from the repo root:

    python scripts/validation/composition_grid.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from imaihda_sim import (  # noqa: E402
    simulate_intersectional_data,
    fit_imaihda,
    correct_detection_bias,
    sparse_strata_vpc,
    joint_calibrated_vpc,
    vpc_latent,
)
from imaihda_sim.detection import _aggregate_strata, _stratum_score, _logit_inv  # noqa: E402
from imaihda_sim.uncertainty import _var_null_from_counts, _build_calibration_curve  # noqa: E402

N = 3500
INTERACTION_SD = 0.9
N_REP = 12
DETECTION_GRID = (0.4, 0.8, 1.2)
SPARSE_GRID = (False, True)
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "composition_grid.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "composition_grid_results.csv")


def pop_truth(df):
    strata = _aggregate_strata(df)
    add = (
        -2.10
        + 0.20 * pd.to_numeric(strata["sex"]).to_numpy(float)
        + 0.35 * pd.to_numeric(strata["education"]).to_numpy(float)
        + 0.30 * pd.to_numeric(strata["wealth"]).to_numpy(float)
        + 0.25 * pd.to_numeric(strata["rural"]).to_numpy(float)
    )
    resid = df.groupby("stratum", observed=True)["true_stratum_residual"].first().to_numpy(float)
    return vpc_latent(float(np.var(add + resid, ddof=1)))


def composed_sequential(df, delta, seed):
    """Naive composition from the capstone: detection-correct the counts,
    then push through the sparse-strata calibration built under a
    no-detection binomial model."""
    strata = _aggregate_strata(df)
    e_obs = strata["events"].to_numpy(float)
    n_j = strata["n"].to_numpy(float)
    s = _stratum_score(strata, None)
    d0 = float(_logit_inv(2.0))
    d_rel = np.clip(np.asarray(_logit_inv(2.0 - delta * s), float) / d0, 1e-6, 1.0)
    e_corr = np.minimum(n_j, e_obs / d_rel)
    p0 = float(np.clip(e_corr.sum() / n_j.sum(), 1e-6, 1 - 1e-6))
    var_naive = float(_var_null_from_counts(e_corr[None, :], n_j[None, :], p0)[0])
    rng = np.random.default_rng(seed)
    grid, means = _build_calibration_curve(n_j.astype(np.int64), p0, 2.0, 25, 200, rng)
    if var_naive <= means[0]:
        v = 0.0
    elif var_naive >= means[-1]:
        v = float(grid[-1])
    else:
        v = float(np.interp(var_naive, means, grid))
    return vpc_latent(v)


def main():
    rows = []
    for delta in DETECTION_GRID:
        for sparse in SPARSE_GRID:
            arms = {k: [] for k in ("naive", "det_only", "sparse_only", "composed", "joint")}
            for seed in range(N_REP):
                df = simulate_intersectional_data(
                    n=N, interaction_sd=INTERACTION_SD, detection_strength=delta,
                    sparse=sparse, seed=seed,
                )
                truth = pop_truth(df)
                arms["naive"].append(fit_imaihda(df)["vpc_null"] - truth)
                arms["det_only"].append(correct_detection_bias(df, delta=delta)["vpc_null"] - truth)
                arms["sparse_only"].append(
                    sparse_strata_vpc(df, seed=7000 + seed)["vpc_null_corrected"] - truth)
                arms["composed"].append(composed_sequential(df, delta, seed=8000 + seed) - truth)
                arms["joint"].append(
                    joint_calibrated_vpc(df, delta=delta, seed=9000 + seed)["vpc_corrected"] - truth)
            for arm, biases in arms.items():
                rows.append(dict(
                    detection_strength=delta, sparse=sparse, arm=arm, n_rep=N_REP,
                    bias=round(float(np.mean(biases)), 2),
                    rmse=round(float(np.sqrt(np.mean(np.square(biases)))), 2),
                ))

    summary = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(summary.pivot_table(index=["detection_strength", "sparse"], columns="arm",
                               values="bias").round(2).to_string())

    # Generalization checks
    cells = summary.pivot_table(index=["detection_strength", "sparse"], columns="arm", values="bias")
    overcorrect_cells = int((cells["composed"] > cells["det_only"] + 0.5).sum())
    joint_best2 = 0
    for idx in cells.index:
        ranks = cells.loc[idx].abs().rank()
        joint_best2 += ranks["joint"] <= 2
    print(f"\ncomposed overshoots det_only (+0.5pp margin) in {overcorrect_cells}/{len(cells)} cells")
    print(f"joint is best-or-second-best (|bias|) in {joint_best2}/{len(cells)} cells")

    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"wrote {os.path.abspath(CSV_OUT)}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    colors = {"naive": "gray", "det_only": "#21918c", "sparse_only": "#fde725",
              "composed": "#f89540", "joint": "#440154"}
    for ax, sparse in zip(axes, SPARSE_GRID):
        sub = summary[summary["sparse"] == sparse]
        width = 0.15
        x = np.arange(len(DETECTION_GRID))
        for k, arm in enumerate(colors):
            vals = [sub[(sub.detection_strength == d) & (sub.arm == arm)]["bias"].iloc[0]
                    for d in DETECTION_GRID]
            ax.bar(x + (k - 2) * width, vals, width=width, label=arm, color=colors[arm])
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"δ={d}" for d in DETECTION_GRID])
        ax.set_title(f"sparse={sparse}")
        ax.set_xlabel("true detection strength")
    axes[0].set_ylabel("Bias in VPC (pp) vs. population truth")
    axes[0].legend(fontsize=8)
    fig.suptitle("Composition grid: 5 estimators across detection x allocation regimes",
                  fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
