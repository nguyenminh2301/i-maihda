"""Regression tests for causal.py: jointly-calibrated VPC and sharp
partial-identification bounds. Thin smoke tests calibrated to the empirical
findings in docs/PHASE3_CAUSAL_IDENTIFICATION.md; the validation runs recorded
there (n_rep=15-20) are the authoritative numbers.
"""
import numpy as np
from itertools import product

from imaihda_sim import (
    simulate_intersectional_data,
    fit_imaihda,
    correct_detection_bias,
    joint_calibrated_vpc,
    vpc_partial_bounds,
    vpc_latent,
)
from imaihda_sim.causal import _min_variance_over_box, _max_variance_over_box
from imaihda_sim.detection import _aggregate_strata
from imaihda_sim.robustness import simulate_severity_dependent_detection


def _pop_truth(df):
    """Population VPC of the true stratum logits (the estimand)."""
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


def test_box_optimizers_globally_correct_small_k():
    rng = np.random.default_rng(1)
    for _ in range(30):
        K = int(rng.integers(3, 9))
        a = rng.normal(0, 2, K)
        b = a + rng.gamma(1.0, 1.0, K)
        exhaustive = max(
            np.var(np.where(np.array(bits) == 1, b, a), ddof=1)
            for bits in product([0, 1], repeat=K)
        )
        v_hi, _ = _max_variance_over_box(a, b)
        assert abs(v_hi - exhaustive) < 1e-9
        v_lo, _ = _min_variance_over_box(a, b)
        assert v_lo <= v_hi + 1e-12
        # feasibility: min no less than 0, and if intervals all overlap at a
        # common point the min must be ~0
        if a.max() <= b.min():
            assert v_lo < 1e-10


def test_bounds_contain_population_truth_even_when_point_id_fails():
    # Outcome-dependent detection (severity_weight=1.0) is the regime where the
    # point correction provably breaks down; the bounds must still cover.
    for seed in range(4):
        df = simulate_severity_dependent_detection(
            n=12000, interaction_sd=0.9, detection_strength=0.8,
            severity_weight=1.0, seed=seed,
        )
        truth = _pop_truth(df)
        b = vpc_partial_bounds(df, delta_max=1.6)
        assert b["vpc_lower"] <= truth <= b["vpc_upper"]


def test_bounds_collapse_when_envelope_is_zero():
    df = simulate_intersectional_data(n=6000, interaction_sd=0.9, detection_strength=0.8, seed=0)
    b = vpc_partial_bounds(df, delta_max=0.0)
    assert abs(b["vpc_upper"] - b["vpc_lower"]) < 1e-9


def test_joint_calibration_beats_detection_only_against_population_truth():
    # Capstone regime (detection + sparsity). Empirically (n_rep=20):
    # naive -16.26, det_only -6.72, joint +0.87 vs the population estimand.
    naive, det_only, joint = [], [], []
    for seed in range(6):
        df = simulate_intersectional_data(
            n=3500, interaction_sd=0.9, detection_strength=0.8, sparse=True, seed=seed
        )
        truth = _pop_truth(df)
        naive.append(fit_imaihda(df)["vpc_null"] - truth)
        det_only.append(correct_detection_bias(df, delta=0.8)["vpc_null"] - truth)
        r = joint_calibrated_vpc(df, delta=0.8, seed=5000 + seed)
        joint.append(r["vpc_corrected"] - truth)
    assert abs(np.mean(joint)) < abs(np.mean(det_only))
    assert abs(np.mean(joint)) < abs(np.mean(naive))


def test_joint_calibration_ci_covers_truth_majority():
    hits = 0
    for seed in range(8):
        df = simulate_intersectional_data(
            n=12000, interaction_sd=0.9, detection_strength=0.8, seed=seed
        )
        truth = _pop_truth(df)
        r = joint_calibrated_vpc(df, delta=0.8, seed=5000 + seed)
        hits += r["ci_lower"] <= truth <= r["ci_upper"]
    assert hits >= 5  # loose smoke floor; driver-scale coverage is 19/20


def test_joint_calibration_reproducible_and_validates_inputs():
    df = simulate_intersectional_data(n=2000, interaction_sd=0.9, detection_strength=0.8, seed=1)
    r1 = joint_calibrated_vpc(df, delta=0.8, seed=7)
    r2 = joint_calibrated_vpc(df, delta=0.8, seed=7)
    assert r1 == r2
    try:
        joint_calibrated_vpc(df, delta=-1)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for negative delta")
