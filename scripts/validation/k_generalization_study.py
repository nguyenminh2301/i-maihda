"""K-generalization study: does sparse_strata_vpc()'s finite-K calibration
hold up when the number of intersectional strata departs from the package's
fixed 36-stratum (2x3x3x2) design?

Deferred twice (METHODS_NOTE_ROBUSTNESS.md Sec. 7, PHASE3_CAUSAL_IDENTIFICATION.md
Sec. 6); completed here. Runs K in {8 (2x2x2), 36 (2x3x3x2), 108 (2x3x3x3x2)}
at the same total n=3500 under sparse allocation, so per-stratum sparsity
worsens mechanically as K grows. Ground truth is analytic (exact population
variance of the K true stratum logits). Prints pass/fail against the
established non-inferiority criterion, writes figures/k_generalization.png
and figures/k_generalization_results.csv. Synthetic data only. Run from the
repo root:

    python scripts/validation/k_generalization_study.py
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

from imaihda_sim import sparse_strata_vpc, vpc_latent  # noqa: E402
from imaihda_sim.robustness import simulate_k_strata  # noqa: E402

N = 3500
N_REP = 20
DIMS = {8: (2, 2, 2), 36: (2, 3, 3, 2), 108: (2, 3, 3, 3, 2)}
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "k_generalization.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "k_generalization_results.csv")

NONINFERIORITY_MARGIN = 0.5


def analytic_truth(df: pd.DataFrame) -> float:
    etas = df.groupby("stratum", observed=True)["true_stratum_logit"].first().to_numpy(float)
    return vpc_latent(float(np.var(etas, ddof=1)))


def main() -> None:
    rows = []
    for K, dims in DIMS.items():
        biases_naive, biases_corr, min_ns = [], [], []
        hits_boot, hits_ti, floor_count = 0, 0, 0
        for seed in range(N_REP):
            df = simulate_k_strata(dims=dims, n=N, sparse=True, seed=seed)
            truth = analytic_truth(df)
            r = sparse_strata_vpc(df, seed=6000 + seed)
            r_ti = sparse_strata_vpc(df, seed=6000 + seed, ci_method="test_inversion")
            biases_naive.append(r["vpc_null_naive"] - truth)
            biases_corr.append(r["vpc_null_corrected"] - truth)
            hits_boot += r["ci_lower"] <= truth <= r["ci_upper"]
            hits_ti += r_ti["ci_lower"] <= truth <= r_ti["ci_upper"]
            floor_count += r["at_floor"]
            min_ns.append(r["min_stratum_n"])
        bias_naive = float(np.mean(biases_naive))
        bias_corr = float(np.mean(biases_corr))
        noninferior = abs(bias_corr) <= abs(bias_naive) + NONINFERIORITY_MARGIN
        reduction = None if abs(bias_naive) < 1.0 else 1 - abs(bias_corr) / abs(bias_naive)
        rows.append(dict(
            K=K, dims=str(dims), n_rep=N_REP,
            mean_min_stratum_n=round(float(np.mean(min_ns)), 1),
            frac_at_floor=round(floor_count / N_REP, 2),
            bias_naive=round(bias_naive, 2),
            bias_corrected=round(bias_corr, 2),
            bias_reduction_pct=None if reduction is None else round(100 * reduction, 1),
            coverage_bootstrap=round(hits_boot / N_REP, 2),
            coverage_test_inversion=round(hits_ti / N_REP, 2),
            noninferior=noninferior,
        ))

    summary = pd.DataFrame(rows)
    pd.set_option("display.width", 140)
    print(summary.to_string(index=False))
    for _, row in summary.iterrows():
        verdict = "PASS" if row["noninferior"] else "BREAKDOWN"
        print(f"K={row['K']}: non-inferiority {verdict}")

    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    x = np.arange(len(summary))
    ax.bar(x - 0.18, summary["bias_naive"], width=0.36, label="naive", color="gray")
    ax.bar(x + 0.18, summary["bias_corrected"], width=0.36, label="corrected", color="#440154")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"K={k}\n{d}" for k, d in zip(summary["K"], summary["dims"])], fontsize=8)
    ax.set_ylabel("Bias in VPC (pp) vs. analytic truth")
    ax.set_title("K-generalization: sparse-strata calibration across stratum designs\n"
                  f"(fixed n={N}, sparse allocation; per-stratum sparsity worsens with K)",
                  fontweight="bold", fontsize=10)
    for i, cov in enumerate(summary["coverage_test_inversion"]):
        ax.annotate(f"cov {cov:.0%}", (x[i] + 0.18, summary["bias_corrected"][i]),
                    textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
