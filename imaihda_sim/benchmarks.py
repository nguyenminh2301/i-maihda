from __future__ import annotations

import pandas as pd


def evaluate_benchmarks(results: pd.DataFrame) -> pd.DataFrame:
    """Scenario-level checks that the synthetic demo answers its intended question."""
    r = results.set_index("scenario")
    checks = []

    def add(name: str, passed: bool, value: str, criterion: str):
        checks.append({"check": name, "passed": bool(passed), "value": value, "criterion": criterion})

    add(
        "A_additive_is_additive_dominant",
        r.loc["A", "pcv"] >= 80 and r.loc["A", "vpc_main"] < 1.0,
        f"PCV={r.loc['A','pcv']:.1f}; VPC_main={r.loc['A','vpc_main']:.2f}",
        "PCV >= 80 and VPC_main < 1",
    )
    add(
        "B_interaction_increases_vpc",
        r.loc["B", "vpc_null"] > r.loc["A", "vpc_null"] + 5,
        f"A={r.loc['A','vpc_null']:.2f}; B={r.loc['B','vpc_null']:.2f}",
        "B VPC_null > A VPC_null + 5 percentage points",
    )
    add(
        "B_interaction_leaves_residual_variance",
        r.loc["B", "pcv"] < 70,
        f"PCV={r.loc['B','pcv']:.1f}",
        "PCV < 70",
    )
    add(
        "C_detection_reduces_observed_prevalence",
        r.loc["C", "overall_prevalence"] < r.loc["A", "overall_prevalence"],
        f"A={100*r.loc['A','overall_prevalence']:.1f}%; C={100*r.loc['C','overall_prevalence']:.1f}%",
        "C observed prevalence < A observed prevalence",
    )
    add(
        "D_detection_can_mask_interaction_vpc",
        r.loc["D", "vpc_null"] < r.loc["B", "vpc_null"],
        f"B={r.loc['B','vpc_null']:.2f}; D={r.loc['D','vpc_null']:.2f}",
        "D VPC_null < B VPC_null despite same residual-interaction SD",
    )
    add(
        "E_sparse_strata_are_flagged",
        r.loc["E", "min_stratum_n"] < r.loc["B", "min_stratum_n"],
        f"B min_n={r.loc['B','min_stratum_n']:.0f}; E min_n={r.loc['E','min_stratum_n']:.0f}",
        "E minimum stratum size < B minimum stratum size",
    )
    return pd.DataFrame(checks)
