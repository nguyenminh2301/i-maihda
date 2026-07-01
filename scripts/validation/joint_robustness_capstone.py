"""Joint capstone: detection bias and sparsity co-occurring (M2-3).

The most realistic worst case for I-MAIHDA in practice: strata are both
genuinely sparse AND the outcome is under-detected in a SES-patterned way
(Scenario D+E combined). This script asks the question neither Study 1 nor
Study 2 alone can answer: do the two corrections COMPOSE? Runs four arms on
the *same* data:

  1. naive              -- no correction (fit_imaihda on observed y)
  2. detection-only      -- correct_detection_bias() alone (ignores sparsity)
  3. sparse-only         -- sparse_strata_vpc() applied directly to the
                            detection-biased observed y (its calibration
                            model assumes no detection process at all)
  4. both (composed)     -- detection-corrected stratum counts, THEN the same
                            small-sample calibration used by sparse_strata_vpc

Prints a 2x2-style decomposition table against the analytic ground truth
(VPC on y_true), writes figures/joint_robustness.png and
figures/joint_robustness_results.csv. Uses synthetic data only. Run from the
repo root:

    python scripts/validation/joint_robustness_capstone.py
"""
from __future__ import annotations

import math
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
    vpc_latent,
)
from imaihda_sim.detection import _aggregate_strata as _detection_aggregate  # noqa: E402
from imaihda_sim.detection import _stratum_score, _logit_inv  # noqa: E402
from imaihda_sim.uncertainty import _var_null_from_counts, _build_calibration_curve  # noqa: E402

N = 3500
INTERACTION_SD = 0.9
DETECTION_STRENGTH = 0.8
N_REP = 20
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "joint_robustness.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "joint_robustness_results.csv")

NONINFERIORITY_MARGIN = 0.5


def _vpc_on_true(df: pd.DataFrame) -> float:
    dft = df.copy()
    dft["y"] = df["y_true"]
    return fit_imaihda(dft)["vpc_null"]


def _detection_corrected_counts(df: pd.DataFrame, delta: float, baseline_logit: float = 2.0):
    """Reconstruct implied-true stratum event counts (same logic as
    correct_detection_bias's internal correction step), returned as raw
    arrays for composition with the sparse-strata calibration."""
    strata = _detection_aggregate(df)
    e_obs = strata["events"].to_numpy(dtype=float)
    n = strata["n"].to_numpy(dtype=float)
    s = _stratum_score(strata, None)
    d0 = float(_logit_inv(baseline_logit))
    d = np.asarray(_logit_inv(baseline_logit - delta * s), dtype=float) / d0
    d = np.clip(d, 1e-6, 1.0)
    e_true = np.minimum(n, e_obs / d)
    return e_true, n


def _composed_correction(df: pd.DataFrame, delta: float, seed: int) -> float:
    """Detection correction, then sparse-strata small-sample calibration
    applied to the detection-corrected counts."""
    e_true, n_j = _detection_corrected_counts(df, delta)
    p_overall = float(np.clip(e_true.sum() / n_j.sum(), 1e-6, 1 - 1e-6))
    var_naive = float(_var_null_from_counts(e_true[None, :], n_j[None, :], p_overall)[0])

    rng = np.random.default_rng(seed)
    grid, means = _build_calibration_curve(n_j.astype(np.int64), p_overall, sigma2_max=2.0,
                                            n_grid=25, n_sim=200, rng=rng)
    if var_naive <= means[0]:
        var_corrected = 0.0
    elif var_naive >= means[-1]:
        var_corrected = float(grid[-1])
    else:
        var_corrected = float(np.interp(var_naive, means, grid))
    return vpc_latent(var_corrected)


def run() -> pd.DataFrame:
    rows = []
    for seed in range(N_REP):
        df = simulate_intersectional_data(
            n=N, interaction_sd=INTERACTION_SD, detection_strength=DETECTION_STRENGTH,
            sparse=True, seed=seed,
        )
        vpc_true = _vpc_on_true(df)
        vpc_naive = fit_imaihda(df)["vpc_null"]
        vpc_detection_only = correct_detection_bias(df, delta=DETECTION_STRENGTH)["vpc_null"]
        vpc_sparse_only = sparse_strata_vpc(df, seed=3000 + seed)["vpc_null_corrected"]
        vpc_both = _composed_correction(df, DETECTION_STRENGTH, seed=4000 + seed)
        rows.append(dict(seed=seed, vpc_true=vpc_true, naive=vpc_naive,
                          detection_only=vpc_detection_only, sparse_only=vpc_sparse_only,
                          both=vpc_both))
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm in ("naive", "detection_only", "sparse_only", "both"):
        bias = float((df[arm] - df["vpc_true"]).mean())
        rmse = float(np.sqrt(((df[arm] - df["vpc_true"]) ** 2).mean()))
        rows.append(dict(arm=arm, mean_vpc_true=round(float(df["vpc_true"].mean()), 2),
                          bias=round(bias, 2), rmse=round(rmse, 2)))
    out = pd.DataFrame(rows)
    naive_bias = out.loc[out["arm"] == "naive", "bias"].iloc[0]
    out["noninferior_vs_naive"] = out["bias"].abs() <= abs(naive_bias) + NONINFERIORITY_MARGIN
    return out


def main() -> None:
    df = run()
    summary = summarize(df)
    pd.set_option("display.width", 140)
    print(summary.to_string(index=False))

    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    df.to_csv(CSV_OUT.replace("_results.csv", "_raw.csv"), index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    fig, ax = plt.subplots(figsize=(7, 5))
    labels = summary["arm"].tolist()
    biases = summary["bias"].tolist()
    colors = ["gray", "#21918c", "#fde725", "#440154"]
    ax.bar(labels, biases, color=colors)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Bias in VPC (pp) vs. analytic truth")
    ax.set_title("Joint capstone: detection bias + sparsity co-occurring\n"
                  "(does correction composition help?)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
