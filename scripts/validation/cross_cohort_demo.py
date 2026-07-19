"""Cross-cohort VPC gap decomposition — reproduces docs/PHASE4_CROSS_COHORT_DECOMPOSITION.md.

Directly operationalizes the PhD exposé's central methodological question:
when can a raw HIC-MIC difference in VPC be trusted as a genuine structural
difference, and when is it an artefact of prevalence, sparse strata, selective
attrition, or differential detection?

Three self-validating scenarios (the true structural gap is known by
construction), plus a stacked-bar figure. Synthetic data only. Run from the
repo root:

    python scripts/validation/cross_cohort_demo.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from imaihda_sim import CohortSpec, cross_cohort_decomposition  # noqa: E402
from imaihda_sim.cohort import CHANNELS  # noqa: E402

FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "cross_cohort_decomposition.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "cross_cohort_decomposition_results.csv")

SCENARIOS = {
    "null_structure": (
        # Identical structure; MIC differs only in the four artefacts.
        CohortSpec(prevalence_shift=-2.10, interaction_sd=0.30, structure_seed=999),
        CohortSpec(prevalence_shift=-3.50, interaction_sd=0.30, structure_seed=999,
                   sparse=True, detection_strength=0.8, attrition_strength=0.8),
    ),
    "pure_structure": (
        # Genuine structural difference, no artefact difference.
        CohortSpec(interaction_sd=0.20, structure_seed=999),
        CohortSpec(interaction_sd=0.95, structure_seed=999),
    ),
    "masked_structure": (
        # Genuine structural difference HIDDEN by opposing artefacts: raw gap
        # is negative, decomposition unmasks a positive structural share.
        CohortSpec(prevalence_shift=-2.10, interaction_sd=0.20, structure_seed=999),
        CohortSpec(prevalence_shift=-2.60, interaction_sd=0.95, structure_seed=999,
                   sparse=True, detection_strength=0.5, attrition_strength=0.5),
    ),
}


def main() -> None:
    import pandas as pd

    rows = []
    results = {}
    for name, (a, b) in SCENARIOS.items():
        r = cross_cohort_decomposition(a, b, metric="vpc", n_rep=12, n_boot=40, seed=100)
        results[name] = r
        rows.append(dict(
            scenario=name,
            vpc_a=round(r["value_a"], 2), vpc_b=round(r["value_b"], 2),
            gap=round(r["gap"], 2),
            **{f"share_{c}": round(r["shares"][c], 2) for c in CHANNELS},
            structure_ci_lo=round(r["structure_ci"][0], 2),
            structure_ci_hi=round(r["structure_ci"][1], 2),
            distinguishable=r["structure_distinguishable_from_zero"],
            true_structural_gap=round(r["true_structural_gap"], 2),
        ))

    summary = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    print(summary.to_string(index=False))
    print("\nInterpretation:")
    print("  null_structure  -> structure share ~0, NOT distinguishable: the whole gap is artefact.")
    print("  pure_structure  -> structure share == gap, artefacts ~0, recovers the true gap.")
    print("  masked_structure-> raw gap negative, but structure share positive & distinguishable:")
    print("                      a genuinely more-intersectional cohort hidden by artefacts.")

    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    # Stacked contribution figure.
    colors = {"prevalence": "#bdbdbd", "sparsity": "#fde725", "detection": "#21918c",
              "attrition": "#f89540", "structure": "#440154"}
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(SCENARIOS))
    for i, (name, r) in enumerate(results.items()):
        pos = 0.0
        neg = 0.0
        for c in ["prevalence", "sparsity", "detection", "attrition", "structure"]:
            val = r["shares"][c]
            base = pos if val >= 0 else neg
            ax.bar(i, val, bottom=base, color=colors[c], edgecolor="white", linewidth=0.5,
                   label=c if i == 0 else None)
            if val >= 0:
                pos += val
            else:
                neg += val
        ax.plot(i, r["gap"], "D", color="black", markersize=8, zorder=5,
                label="observed gap" if i == 0 else None)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(list(SCENARIOS), fontsize=9)
    ax.set_ylabel("Contribution to VPC gap (VPC_MIC − VPC_HIC, pp)")
    ax.set_title("Cross-cohort VPC gap decomposition\n"
                  "(structural residual vs. four artefact channels; ♦ = observed gap)",
                  fontweight="bold", fontsize=11)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
