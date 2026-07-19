"""Tests for the cross-cohort VPC/PCV Shapley decomposition (Phase 4) and the
selective-attrition generator. Calibrated to the self-validating scenarios in
scripts/validation/cross_cohort_demo.py; the decomposition's correctness is
checkable because the true structural gap is known by construction.
"""
import math

import numpy as np
import pytest

from imaihda_sim import (
    simulate_intersectional_data,
    CohortSpec,
    cross_cohort_decomposition,
)
from imaihda_sim.robustness import simulate_selective_attrition
from imaihda_sim.cohort import CHANNELS, _ARTEFACT_CHANNELS


# ── attrition generator ───────────────────────────────────────────────────────

def test_attrition_zero_matches_original_mechanism():
    for kw in [dict(n=3000, interaction_sd=0.9, sparse=False, seed=7),
               dict(n=3500, interaction_sd=0.9, sparse=True, seed=42)]:
        a = simulate_intersectional_data(detection_strength=0.0, **kw)
        b = simulate_selective_attrition(attrition_strength=0.0, **kw)
        cols = ["sex", "education", "wealth", "rural", "stratum", "y", "y_true",
                "true_stratum_residual"]
        assert a[cols].equals(b[cols])


def test_attrition_drops_rows_and_is_ses_patterned():
    full = simulate_selective_attrition(n=8000, interaction_sd=0.9, attrition_strength=0.0, seed=1)
    thin = simulate_selective_attrition(n=8000, interaction_sd=0.9, attrition_strength=0.9, seed=1)
    assert len(thin) < len(full)
    assert thin["retention_probability"].mean() < 1.0
    # Retention is SES-patterned: the most-advantaged stratum keeps a larger
    # share of its members than the most-disadvantaged one.
    def kept_fraction(edu_level):
        f = (full["education"].astype(int) == edu_level).sum()
        t = (thin["education"].astype(int) == edu_level).sum()
        return t / f
    assert kept_fraction(0) > kept_fraction(2)


# ── decomposition: shares are exact and sum to the gap ────────────────────────

def test_shares_sum_exactly_to_gap():
    a = CohortSpec(prevalence_shift=-2.10, interaction_sd=0.30, structure_seed=999)
    b = CohortSpec(prevalence_shift=-3.50, interaction_sd=0.30, structure_seed=999,
                   sparse=True, detection_strength=0.8, attrition_strength=0.8)
    r = cross_cohort_decomposition(a, b, metric="vpc", n_rep=8, n_boot=10, seed=100)
    assert set(r["shares"]) == set(CHANNELS)
    assert abs(sum(r["shares"].values()) - r["gap"]) < 1e-9


def test_null_structure_gives_zero_structural_share():
    # Identical structure, differing nuisances: the whole gap must be artefact;
    # the structure share is exactly zero (flipping identical structure is a
    # no-op) and is not distinguishable from zero.
    a = CohortSpec(prevalence_shift=-2.10, interaction_sd=0.30, structure_seed=999)
    b = CohortSpec(prevalence_shift=-3.50, interaction_sd=0.30, structure_seed=999,
                   sparse=True, detection_strength=0.8, attrition_strength=0.8)
    r = cross_cohort_decomposition(a, b, metric="vpc", n_rep=10, n_boot=20, seed=100)
    assert abs(r["structure_share"]) < 1e-9
    assert r["structure_distinguishable_from_zero"] is False
    assert abs(r["true_structural_gap"]) < 1e-9
    assert abs(r["gap"]) > 3.0  # there is a real (artefact) gap to explain


def test_pure_structure_recovers_true_gap_with_zero_artefact():
    a = CohortSpec(interaction_sd=0.20, structure_seed=999)
    b = CohortSpec(interaction_sd=0.95, structure_seed=999)
    r = cross_cohort_decomposition(a, b, metric="vpc", n_rep=12, n_boot=20, seed=100)
    assert abs(r["artefact_total"]) < 0.5              # no nuisance difference
    assert abs(r["structure_share"] - r["true_structural_gap"]) < 2.0
    assert r["structure_distinguishable_from_zero"] is True


def test_masked_structure_is_unmasked():
    # Genuine structural difference hidden by opposing nuisance channels: the
    # raw gap is negative, but the decomposition recovers a positive,
    # distinguishable structural share.
    a = CohortSpec(prevalence_shift=-2.10, interaction_sd=0.20, structure_seed=999)
    b = CohortSpec(prevalence_shift=-2.60, interaction_sd=0.95, structure_seed=999,
                   sparse=True, detection_strength=0.5, attrition_strength=0.5)
    r = cross_cohort_decomposition(a, b, metric="vpc", n_rep=12, n_boot=20, seed=100)
    assert r["gap"] < 0.0                    # MIC looks LOWER at face value
    assert r["structure_share"] > 3.0        # but structure is genuinely higher
    assert r["structure_distinguishable_from_zero"] is True


def test_pcv_metric_flags_unreliable_under_heavy_nuisance():
    a = CohortSpec(prevalence_shift=-2.10, interaction_sd=0.20, structure_seed=999)
    b = CohortSpec(prevalence_shift=-2.60, interaction_sd=0.95, structure_seed=999,
                   sparse=True, detection_strength=0.5, attrition_strength=0.5)
    r = cross_cohort_decomposition(a, b, metric="pcv", n_rep=10, n_boot=20, seed=100)
    assert r["reliable"] is False
    assert r["structure_distinguishable_from_zero"] is None

    # In a mild regime PCV is well-defined and reliable.
    a2 = CohortSpec(interaction_sd=0.20, structure_seed=999)
    b2 = CohortSpec(interaction_sd=0.95, structure_seed=999)
    r2 = cross_cohort_decomposition(a2, b2, metric="pcv", n_rep=10, n_boot=20, seed=100)
    assert r2["reliable"] is True


def test_input_validation():
    a = CohortSpec()
    b = CohortSpec(interaction_sd=0.9)
    with pytest.raises(ValueError):
        cross_cohort_decomposition(a, b, metric="bogus")
    with pytest.raises(ValueError):
        cross_cohort_decomposition(a, b, level=1.5)
