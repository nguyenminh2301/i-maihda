"""Robustness/simulation study for correct_detection_bias() under violated
assumptions (Study 1 of the methods-note robustness study).

Runs three scenarios against an analytic-adjacent ground truth (VPC computed
on the simulator's own y_true column):

  M1-1  Score functional-form misspecification (analyst's score= formula
        differs from the true detection-driving covariates)
  M1-2  Detection depends on the true outcome directly (verification bias) --
        an identification-breaking violation, best-case score/delta otherwise
  M1-3  baseline_logit (detection curvature) misspecification

Prints a results table with pass/fail against the success criteria in
docs/METHODS_NOTE_ROBUSTNESS.md, writes figures/detection_robustness.png and
figures/detection_robustness_results.csv. Uses synthetic data only. Run from
the repo root:

    python scripts/validation/detection_bias_robustness.py
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

from imaihda_sim import simulate_intersectional_data, fit_imaihda, correct_detection_bias  # noqa: E402
from imaihda_sim.robustness import (  # noqa: E402
    misspecified_score,
    simulate_severity_dependent_detection,
)

N = 12000
INTERACTION_SD = 0.9
TRUE_DELTA = 0.8
N_REP = 20
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "detection_robustness.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "detection_robustness_results.csv")

NONINFERIORITY_MARGIN = 0.5  # pp
BIAS_REDUCTION_TARGET = 0.40
MILD_M1_1_KINDS = {"omit_rural", "quadratic"}   # mild-moderate misspecification
SEVERE_M1_1_KINDS = {"wrong_covariate"}          # severe misspecification


def _vpc_on_true(df: pd.DataFrame) -> float:
    dft = df.copy()
    dft["y"] = df["y_true"]
    return fit_imaihda(dft)["vpc_null"]


def _bias_reduction(bias_naive: float, bias_corrected: float):
    if abs(bias_naive) < 1.0:
        return None  # naive already ~unbiased; ratio uninformative
    return 1.0 - abs(bias_corrected) / abs(bias_naive)


def _classify(bias_naive: float, bias_corrected: float, require_reduction: bool) -> str:
    noninferior = abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN
    reduction = _bias_reduction(bias_naive, bias_corrected)
    if not noninferior:
        return "BREAKDOWN"
    if require_reduction:
        if reduction is not None and reduction >= BIAS_REDUCTION_TARGET:
            return "PASS"
        return "FAIL (below bias-reduction target)"
    if reduction is not None and reduction >= BIAS_REDUCTION_TARGET:
        return "attenuated-but-positive (exceeds mild-target anyway)"
    return "attenuated-but-net-neutral-or-positive"


def run_m1_1() -> pd.DataFrame:
    rows = []
    for seed in range(N_REP):
        df = simulate_intersectional_data(
            n=N, interaction_sd=INTERACTION_SD, detection_strength=TRUE_DELTA, seed=seed
        )
        vpc_true = _vpc_on_true(df)
        vpc_naive = fit_imaihda(df)["vpc_null"]
        vpc_correct_score = correct_detection_bias(df, delta=TRUE_DELTA)["vpc_null"]
        rows.append(dict(scenario="M1-1", arm="correctly_specified", seed=seed,
                          vpc_true=vpc_true, vpc_naive=vpc_naive, vpc_corrected=vpc_correct_score))
        for kind in ("omit_rural", "wrong_covariate", "quadratic"):
            score = misspecified_score(df, kind)
            vpc_corr = correct_detection_bias(df, delta=TRUE_DELTA, score=score)["vpc_null"]
            rows.append(dict(scenario="M1-1", arm=kind, seed=seed,
                              vpc_true=vpc_true, vpc_naive=vpc_naive, vpc_corrected=vpc_corr))
    return pd.DataFrame(rows)


def run_m1_2() -> pd.DataFrame:
    rows = []
    for severity_weight in (0.0, 0.3, 0.6, 0.8, 1.0):
        for seed in range(N_REP):
            df = simulate_severity_dependent_detection(
                n=N, interaction_sd=INTERACTION_SD, detection_strength=TRUE_DELTA,
                severity_weight=severity_weight, seed=seed,
            )
            vpc_true = _vpc_on_true(df)
            vpc_naive = fit_imaihda(df)["vpc_null"]
            vpc_corr = correct_detection_bias(df, delta=TRUE_DELTA)["vpc_null"]
            rows.append(dict(scenario="M1-2", arm=f"severity_weight={severity_weight}", seed=seed,
                              vpc_true=vpc_true, vpc_naive=vpc_naive, vpc_corrected=vpc_corr))
    return pd.DataFrame(rows)


def run_m1_3() -> pd.DataFrame:
    rows = []
    for baseline_logit in (1.0, 1.5, 2.0, 2.5, 3.0):
        for seed in range(N_REP):
            df = simulate_intersectional_data(
                n=N, interaction_sd=INTERACTION_SD, detection_strength=TRUE_DELTA, seed=seed
            )
            vpc_true = _vpc_on_true(df)
            vpc_naive = fit_imaihda(df)["vpc_null"]
            vpc_corr = correct_detection_bias(df, delta=TRUE_DELTA, baseline_logit=baseline_logit)["vpc_null"]
            rows.append(dict(scenario="M1-3", arm=f"baseline_logit={baseline_logit}", seed=seed,
                              vpc_true=vpc_true, vpc_naive=vpc_naive, vpc_corrected=vpc_corr))
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, require_reduction_arms: set) -> pd.DataFrame:
    out = []
    for (scenario, arm), g in df.groupby(["scenario", "arm"], sort=False):
        bias_naive = float((g["vpc_naive"] - g["vpc_true"]).mean())
        bias_corrected = float((g["vpc_corrected"] - g["vpc_true"]).mean())
        reduction = _bias_reduction(bias_naive, bias_corrected)
        classification = _classify(bias_naive, bias_corrected, arm in require_reduction_arms)
        out.append(dict(
            scenario=scenario, arm=arm, n_rep=len(g),
            mean_vpc_true=round(float(g["vpc_true"].mean()), 2),
            bias_naive=round(bias_naive, 2),
            bias_corrected=round(bias_corrected, 2),
            bias_reduction_pct=None if reduction is None else round(100 * reduction, 1),
            classification=classification,
        ))
    return pd.DataFrame(out)


def main() -> None:
    m1_1 = run_m1_1()
    m1_2 = run_m1_2()
    m1_3 = run_m1_3()

    require_reduction = {"omit_rural", "quadratic", "severity_weight=0.0", "severity_weight=0.3"}
    summary = pd.concat([
        summarize(m1_1, require_reduction),
        summarize(m1_2, require_reduction),
        summarize(m1_3, require_reduction),
    ], ignore_index=True)

    pd.set_option("display.width", 140)
    print(summary.to_string(index=False))
    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    # Figure: M1-1 (score misspecification) bar chart + M1-2 (severity) degradation line.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    m1_1_summary = summarize(m1_1, require_reduction)
    labels = m1_1_summary["arm"].tolist()
    biases_naive = m1_1_summary["bias_naive"].tolist()
    biases_corr = m1_1_summary["bias_corrected"].tolist()
    xpos = np.arange(len(labels))
    ax.bar(xpos - 0.18, biases_naive, width=0.36, label="naive", color="gray")
    ax.bar(xpos + 0.18, biases_corr, width=0.36, label="corrected", color="#440154")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Bias in VPC (pp)")
    ax.set_title("M1-1: score misspecification")
    ax.legend(fontsize=8)

    ax = axes[1]
    m1_2_summary = summarize(m1_2, require_reduction)
    sw = [float(a.split("=")[1]) for a in m1_2_summary["arm"]]
    ax.plot(sw, m1_2_summary["bias_naive"], "o-", color="gray", label="naive")
    ax.plot(sw, m1_2_summary["bias_corrected"], "o-", color="#440154", label="corrected")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("severity_weight (outcome-dependent detection strength)")
    ax.set_ylabel("Bias in VPC (pp)")
    ax.set_title("M1-2: detection depends on true outcome")
    ax.legend(fontsize=8)

    fig.suptitle("Detection-bias correction under assumption violations", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
