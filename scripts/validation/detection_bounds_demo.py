"""Generate figures/detection_bounds.png — detection-bias sensitivity bounds.

Demonstrates the imaihda-only quantitative bias analysis: how the observed
null-model VPC could shift once SES-patterned under-detection is accounted for.
Uses synthetic data only. Run from the repo root:

    python scripts/validation/detection_bounds_demo.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Allow running from the repo root without installing the package.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from imaihda_sim import (  # noqa: E402
    simulate_intersectional_data,
    fit_imaihda,
    vpc_detection_bounds,
)

DELTA_TRUE = 0.8
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "detection_bounds.png")


def main() -> None:
    df = simulate_intersectional_data(
        n=12000, interaction_sd=0.9, detection_strength=DELTA_TRUE, seed=7
    )

    dft = df.copy()
    dft["y"] = df["y_true"]
    vpc_true = fit_imaihda(dft)["vpc_null"]
    vpc_obs = fit_imaihda(df)["vpc_null"]

    bounds = vpc_detection_bounds(df, delta_max=1.2, n_grid=25)
    lo, hi = bounds["vpc_null"].min(), bounds["vpc_null"].max()

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.fill_between(bounds["delta"], lo, hi, color="#440154", alpha=0.07,
                    label=f"Plausible range [{lo:.1f}, {hi:.1f}]%")
    ax.plot(bounds["delta"], bounds["vpc_null"], color="#440154", lw=2,
            label="Corrected null-model VPC")
    ax.axhline(vpc_obs, ls=":", color="gray",
               label=f"Observed VPC = {vpc_obs:.1f}% (delta = 0)")
    ax.axhline(vpc_true, ls="--", color="#21918c",
               label=f"True VPC (y_true) = {vpc_true:.1f}%")
    ax.scatter([0], [vpc_obs], color="#440154", zorder=5)

    ax.set_xlabel("Assumed SES-patterned under-detection strength (delta)")
    ax.set_ylabel("Null-model VPC (%)")
    ax.set_title("Detection-bias sensitivity bounds on VPC", fontweight="bold")
    ax.text(0.5, -0.16, "imaihda-only quantitative bias analysis — no equivalent in CRAN MAIHDA",
            transform=ax.transAxes, ha="center", fontsize=8, color="gray")
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.abspath(OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(OUT)}  (obs={vpc_obs:.2f}, true={vpc_true:.2f}, bounds=[{lo:.2f},{hi:.2f}])")


if __name__ == "__main__":
    main()
