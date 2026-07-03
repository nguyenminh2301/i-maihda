"""Reproduce the headline numbers of docs/PHASE3_CAUSAL_IDENTIFICATION.md.

Runs the jointly-calibrated VPC estimator and the sharp partial-identification
bounds against the POPULATION estimand across three regimes, plus the
bounds-coverage check in the outcome-dependent-detection regime where point
identification provably fails. Writes figures/phase3_causal.png and
figures/phase3_causal_results.csv. Synthetic data only. Run from repo root:

    python scripts/validation/phase3_causal_demo.py
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
    joint_calibrated_vpc,
    vpc_partial_bounds,
    vpc_latent,
)
from imaihda_sim.detection import _aggregate_strata  # noqa: E402
from imaihda_sim.robustness import simulate_severity_dependent_detection  # noqa: E402

N_REP = 20
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "phase3_causal.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "phase3_causal_results.csv")


def pop_truth(df):
    strata = _aggregate_strata(df)
    add = (
        -2.10
        + 0.20 * strata["sex"].to_numpy(float)
        + 0.35 * strata["education"].to_numpy(float)
        + 0.30 * strata["wealth"].to_numpy(float)
        + 0.25 * strata["rural"].to_numpy(float)
    )
    resid = df.groupby("stratum", observed=True)["true_stratum_residual"].first().to_numpy(float)
    return vpc_latent(float(np.var(add + resid, ddof=1)))


def main():
    regimes = [
        ("capstone D+E (n=3500, sparse)", dict(n=3500, interaction_sd=0.9, detection_strength=0.8, sparse=True)),
        ("harsher (n=2000, shift=-3, sparse)", dict(n=2000, prevalence_shift=-3.0, interaction_sd=0.9,
                                                     detection_strength=0.8, sparse=True)),
        ("dense (n=12000)", dict(n=12000, interaction_sd=0.9, detection_strength=0.8, sparse=False)),
    ]
    rows = []
    for name, kw in regimes:
        res = dict(naive=[], det_only=[], joint=[])
        cov = 0
        for seed in range(N_REP):
            df = simulate_intersectional_data(seed=seed, **kw)
            truth = pop_truth(df)
            res["naive"].append(fit_imaihda(df)["vpc_null"] - truth)
            res["det_only"].append(correct_detection_bias(df, delta=kw["detection_strength"])["vpc_null"] - truth)
            r = joint_calibrated_vpc(df, delta=kw["detection_strength"], seed=5000 + seed)
            res["joint"].append(r["vpc_corrected"] - truth)
            cov += r["ci_lower"] <= truth <= r["ci_upper"]
        for arm, biases in res.items():
            rows.append(dict(regime=name, arm=arm,
                              bias=round(float(np.mean(biases)), 2),
                              rmse=round(float(np.sqrt(np.mean(np.square(biases)))), 2),
                              ci_coverage=round(cov / N_REP, 2) if arm == "joint" else None))

    # Bounds coverage where point identification fails (outcome-dependent detection)
    inside = 0
    for seed in range(15):
        df = simulate_severity_dependent_detection(
            n=12000, interaction_sd=0.9, detection_strength=0.8, severity_weight=1.0, seed=seed
        )
        truth = pop_truth(df)
        b = vpc_partial_bounds(df, delta_max=1.6)
        inside += b["vpc_lower"] <= truth <= b["vpc_upper"]
    rows.append(dict(regime="M1-2 severity=1.0 (point-ID broken)", arm="bounds delta_max=1.6",
                      bias=None, rmse=None, ci_coverage=round(inside / 15, 2)))

    summary = pd.DataFrame(rows)
    pd.set_option("display.width", 140)
    print(summary.to_string(index=False))
    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    fig, ax = plt.subplots(figsize=(9, 5))
    arms = ["naive", "det_only", "joint"]
    colors = {"naive": "gray", "det_only": "#21918c", "joint": "#440154"}
    width = 0.25
    xs = np.arange(len(regimes))
    for k, arm in enumerate(arms):
        vals = [summary[(summary.regime == name) & (summary.arm == arm)]["bias"].iloc[0]
                for name, _ in regimes]
        ax.bar(xs + (k - 1) * width, vals, width=width, label=arm, color=colors[arm])
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([n for n, _ in regimes], fontsize=8)
    ax.set_ylabel("Bias in VPC (pp) vs. POPULATION truth")
    ax.set_title("Jointly-calibrated VPC vs. naive and detection-only correction",
                  fontweight="bold")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
