"""PCV-focused robustness pass (previously deferred in
METHODS_NOTE_ROBUSTNESS.md Sec. 7).

The robustness studies so far tracked the null-model VPC. The PCV — the
proportional change in variance from the null to the additive main-effects
model, the literature's standard evidence for "intersectional effects" — is a
*ratio* of two variances and can react differently. This script tracks the
PCV itself:

  Part A  naive vs detection-corrected PCV under Study-1 violations
          (score misspecification arms; outcome-dependent detection arms)
  Part B  naive PCV under pure sparsity (no PCV correction exists in the
          package — this quantifies the open gap honestly)

Population ground truth is analytic: V_null = Var(eta_true) and
V_main = Var(true_stratum_residual) (the DGP's residual is exactly the
component not explained by additive main effects), so
PCV_pop = 100 * (V_null - V_main) / V_null.

Writes figures/pcv_robustness.png + figures/pcv_robustness_results.csv.
Synthetic data only. Run from the repo root:

    python scripts/validation/pcv_robustness.py
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
from imaihda_sim.detection import _aggregate_strata  # noqa: E402
from imaihda_sim.robustness import misspecified_score, simulate_severity_dependent_detection  # noqa: E402

N = 12000
INTERACTION_SD = 0.9
TRUE_DELTA = 0.8
N_REP = 20
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "pcv_robustness.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "pcv_robustness_results.csv")

NONINFERIORITY_MARGIN = 5.0  # pp on the PCV scale (PCV is noisier than VPC)


def pop_truth_pcv(df: pd.DataFrame) -> float:
    strata = _aggregate_strata(df)
    add = (
        -2.10
        + 0.20 * pd.to_numeric(strata["sex"]).to_numpy(float)
        + 0.35 * pd.to_numeric(strata["education"]).to_numpy(float)
        + 0.30 * pd.to_numeric(strata["wealth"]).to_numpy(float)
        + 0.25 * pd.to_numeric(strata["rural"]).to_numpy(float)
    )
    resid = df.groupby("stratum", observed=True)["true_stratum_residual"].first().to_numpy(float)
    v_null = float(np.var(add + resid, ddof=1))
    v_main = float(np.var(resid, ddof=1))
    return 100.0 * (v_null - v_main) / v_null


def run_part_a() -> pd.DataFrame:
    rows = []
    # Score misspecification arms (true DGP fixed; analyst's score varies).
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=N, interaction_sd=INTERACTION_SD,
                                          detection_strength=TRUE_DELTA, seed=seed)
        truth = pop_truth_pcv(df)
        pcv_naive = fit_imaihda(df)["pcv"]
        rows.append(dict(part="A", arm="correctly_specified", seed=seed, pcv_true=truth,
                          pcv_naive=pcv_naive,
                          pcv_corrected=correct_detection_bias(df, delta=TRUE_DELTA)["pcv"]))
        for kind in ("omit_rural", "wrong_covariate", "quadratic"):
            score = misspecified_score(df, kind)
            rows.append(dict(part="A", arm=kind, seed=seed, pcv_true=truth,
                              pcv_naive=pcv_naive,
                              pcv_corrected=correct_detection_bias(df, delta=TRUE_DELTA, score=score)["pcv"]))
    # Outcome-dependent detection arms (correct score/delta; identification degrades).
    for severity_weight in (0.0, 0.6, 1.0):
        for seed in range(N_REP):
            df = simulate_severity_dependent_detection(
                n=N, interaction_sd=INTERACTION_SD, detection_strength=TRUE_DELTA,
                severity_weight=severity_weight, seed=seed)
            truth = pop_truth_pcv(df)
            rows.append(dict(part="A", arm=f"severity_weight={severity_weight}", seed=seed,
                              pcv_true=truth, pcv_naive=fit_imaihda(df)["pcv"],
                              pcv_corrected=correct_detection_bias(df, delta=TRUE_DELTA)["pcv"]))
    return pd.DataFrame(rows)


def run_part_b() -> pd.DataFrame:
    rows = []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=3500, interaction_sd=INTERACTION_SD,
                                          sparse=True, seed=seed)
        truth = pop_truth_pcv(df)
        rows.append(dict(part="B", arm="sparse_no_detection", seed=seed, pcv_true=truth,
                          pcv_naive=fit_imaihda(df)["pcv"], pcv_corrected=np.nan))
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (part, arm), g in df.groupby(["part", "arm"], sort=False):
        truth = float(g["pcv_true"].mean())
        bias_naive = float((g["pcv_naive"] - g["pcv_true"]).mean())
        has_corr = g["pcv_corrected"].notna().all()
        bias_corr = float((g["pcv_corrected"] - g["pcv_true"]).mean()) if has_corr else None
        noninferior = None
        if has_corr:
            noninferior = abs(bias_corr) <= abs(bias_naive) + NONINFERIORITY_MARGIN
        out.append(dict(part=part, arm=arm, n_rep=len(g),
                         mean_pcv_true=round(truth, 1),
                         bias_naive=round(bias_naive, 1),
                         bias_corrected=None if bias_corr is None else round(bias_corr, 1),
                         noninferior=noninferior))
    return pd.DataFrame(out)


def main() -> None:
    a = run_part_a()
    b = run_part_b()
    summary = pd.concat([summarize(a), summarize(b)], ignore_index=True)
    pd.set_option("display.width", 140)
    print(summary.to_string(index=False))

    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    sa = summarize(a)
    fig, ax = plt.subplots(figsize=(10, 4.8))
    x = np.arange(len(sa))
    ax.bar(x - 0.18, sa["bias_naive"], width=0.36, label="naive", color="gray")
    ax.bar(x + 0.18, [v if v is not None else 0 for v in sa["bias_corrected"]],
           width=0.36, label="corrected", color="#440154")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sa["arm"], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Bias in PCV (pp) vs. population truth")
    ax.set_title("PCV under detection-bias assumption violations", fontweight="bold")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
