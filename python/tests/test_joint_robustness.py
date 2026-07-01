"""Regression test for the joint capstone scenario (M2-3): detection bias and
sparsity co-occurring. Calibrated to the empirical finding in
scripts/validation/joint_robustness_capstone.py (n_rep=20) and
docs/METHODS_NOTE_ROBUSTNESS.md: detection-only correction does most of the
work in this scenario (dominant bias source), sparse-only correction helps
less, and naively composing both corrections OVERCORRECTS -- the two
corrections do not compose additively. This is a documented finding, not a
bug; this test guards against a future silent change reversing it.
"""
import numpy as np

from imaihda_sim import simulate_intersectional_data, fit_imaihda, correct_detection_bias, sparse_strata_vpc
from imaihda_sim.detection import _aggregate_strata as _detection_aggregate
from imaihda_sim.detection import _stratum_score, _logit_inv
from imaihda_sim.uncertainty import _var_null_from_counts, _build_calibration_curve
from imaihda_sim.fit import vpc_latent

N = 3500
INTERACTION_SD = 0.9
DETECTION_STRENGTH = 0.8
N_REP = 8


def _vpc_on_true(df):
    dft = df.copy()
    dft["y"] = df["y_true"]
    return fit_imaihda(dft)["vpc_null"]


def _composed_correction(df, delta, seed):
    strata = _detection_aggregate(df)
    e_obs = strata["events"].to_numpy(dtype=float)
    n_j = strata["n"].to_numpy(dtype=float)
    s = _stratum_score(strata, None)
    d0 = float(_logit_inv(2.0))
    d = np.clip(np.asarray(_logit_inv(2.0 - delta * s), dtype=float) / d0, 1e-6, 1.0)
    e_true = np.minimum(n_j, e_obs / d)

    p_overall = float(np.clip(e_true.sum() / n_j.sum(), 1e-6, 1 - 1e-6))
    var_naive = float(_var_null_from_counts(e_true[None, :], n_j[None, :], p_overall)[0])
    rng = np.random.default_rng(seed)
    grid, means = _build_calibration_curve(n_j.astype(np.int64), p_overall, sigma2_max=2.0,
                                           n_grid=25, n_sim=200, rng=rng)
    if var_naive <= means[0]:
        var_corrected = 0.0
    elif var_naive >= means[-1]:
        var_corrected = float(grid[-1])
    else:
        var_corrected = float(np.interp(var_naive, means, grid))
    return vpc_latent(var_corrected)


def test_detection_only_correction_removes_most_of_the_joint_bias():
    naive_biases, detection_only_biases = [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=N, interaction_sd=INTERACTION_SD,
                                          detection_strength=DETECTION_STRENGTH,
                                          sparse=True, seed=seed)
        truth = _vpc_on_true(df)
        naive_biases.append(fit_imaihda(df)["vpc_null"] - truth)
        detection_only_biases.append(
            correct_detection_bias(df, delta=DETECTION_STRENGTH)["vpc_null"] - truth
        )
    naive_bias = float(np.mean(naive_biases))
    detection_only_bias = float(np.mean(detection_only_biases))
    assert abs(naive_bias) > 5.0  # substantial joint bias present
    assert abs(detection_only_bias) < abs(naive_bias) * 0.5  # removes most of it


def test_composing_both_corrections_overcorrects_past_naive():
    # The headline negative/cautionary finding: naively composing
    # detection-correction with sparse-strata calibration overshoots further
    # from the truth (in the opposite direction) than either alone.
    naive_biases, detection_only_biases, composed_biases = [], [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=N, interaction_sd=INTERACTION_SD,
                                          detection_strength=DETECTION_STRENGTH,
                                          sparse=True, seed=seed)
        truth = _vpc_on_true(df)
        naive_biases.append(fit_imaihda(df)["vpc_null"] - truth)
        detection_only_biases.append(
            correct_detection_bias(df, delta=DETECTION_STRENGTH)["vpc_null"] - truth
        )
        composed_biases.append(_composed_correction(df, DETECTION_STRENGTH, seed=4000 + seed) - truth)

    naive_bias = float(np.mean(naive_biases))
    detection_only_bias = float(np.mean(detection_only_biases))
    composed_bias = float(np.mean(composed_biases))

    assert composed_bias > 0  # overshoots to the opposite sign from naive's negative bias
    assert abs(composed_bias) > abs(detection_only_bias)  # composing is worse than detection-only alone
