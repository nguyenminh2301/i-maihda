"""Generate figures/sparse_strata_ci.png — bias correction under sparse strata.

Demonstrates the imaihda-only sparse-strata calibration: compares naive vs.
bias-corrected VPC against the analytic ground truth across replications
under genuinely sparse stratum allocation (matching Scenario E). Uses
synthetic data only. Run from the repo root:

    python scripts/validation/sparse_strata_ci_demo.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from imaihda_sim import simulate_intersectional_data, sparse_strata_vpc  # noqa: E402
from imaihda_sim.simulate import _stratum_table  # noqa: E402
from imaihda_sim.fit import vpc_latent  # noqa: E402

INTERACTION_SD = 0.15
N_TRIALS = 12
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "sparse_strata_ci.png")


def analytic_truth(interaction_sd: float = INTERACTION_SD, seed: int = 999) -> float:
    strata = _stratum_table()
    eta_add = (
        -2.10
        + 0.20 * strata["sex"].to_numpy()
        + 0.35 * strata["education"].to_numpy()
        + 0.30 * strata["wealth"].to_numpy()
        + 0.25 * strata["rural"].to_numpy()
    )
    rng_re = np.random.default_rng(seed + 999)
    raw = rng_re.normal(0.0, interaction_sd, size=len(strata))
    raw += 0.90 * (
        (strata["education"].to_numpy() == 2)
        & (strata["wealth"].to_numpy() == 2)
        & (strata["rural"].to_numpy() == 1)
    )
    raw -= 0.60 * (
        (strata["education"].to_numpy() == 0)
        & (strata["wealth"].to_numpy() == 2)
        & (strata["rural"].to_numpy() == 0)
    )
    eta_true = eta_add + (raw - raw.mean())
    return vpc_latent(np.var(eta_true, ddof=1))


def main() -> None:
    truth = analytic_truth()
    naive, corrected, ci_lo, ci_hi = [], [], [], []
    for seed in range(N_TRIALS):
        df = simulate_intersectional_data(
            n=3500, interaction_sd=INTERACTION_SD, sparse=True, seed=seed
        )
        r = sparse_strata_vpc(df, seed=200 + seed)
        naive.append(r["vpc_null_naive"])
        corrected.append(r["vpc_null_corrected"])
        ci_lo.append(r["ci_lower"])
        ci_hi.append(r["ci_upper"])

    x = np.arange(N_TRIALS)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.axhline(truth, color="#21918c", ls="--", lw=1.5, label=f"Analytic truth = {truth:.1f}%")
    ax.errorbar(
        x, corrected,
        yerr=[np.array(corrected) - np.array(ci_lo), np.array(ci_hi) - np.array(corrected)],
        fmt="o", color="#440154", ecolor="#440154", alpha=0.85, capsize=3,
        label="Bias-corrected VPC (95% CI)",
    )
    ax.scatter(x, naive, color="gray", marker="x", s=45, label="Naive VPC (no correction)", zorder=5)

    ax.set_xlabel("Replicate (independent sparse-strata sample, Scenario E)")
    ax.set_ylabel("Null-model VPC (%)")
    ax.set_title("Sparse-strata bias correction: naive vs. corrected VPC", fontweight="bold")
    ax.text(0.5, -0.16, "imaihda-only calibration — no equivalent in CRAN MAIHDA",
            transform=ax.transAxes, ha="center", fontsize=8, color="gray")
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    ax.set_xticks(x)
    fig.tight_layout()
    fig.savefig(os.path.abspath(OUT), dpi=130, bbox_inches="tight")

    naive_gap = abs(np.mean(naive) - truth)
    corr_gap = abs(np.mean(corrected) - truth)
    covered = sum(lo <= truth <= hi for lo, hi in zip(ci_lo, ci_hi))
    print(f"wrote {os.path.abspath(OUT)}")
    print(f"mean naive gap={naive_gap:.2f}  mean corrected gap={corr_gap:.2f}  "
          f"CI coverage={covered}/{N_TRIALS}")


if __name__ == "__main__":
    main()
