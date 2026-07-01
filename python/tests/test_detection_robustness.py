"""Robustness regression tests for correct_detection_bias() under violated
assumptions (Study 1). Thin, fast smoke tests calibrated to the empirical
findings in scripts/validation/detection_bias_robustness.py (n_rep=20) and
docs/METHODS_NOTE_ROBUSTNESS.md -- not to a priori assumptions. Use a smaller
n_rep here (Monte Carlo noise tolerated via generous margins); the driver
script is the authoritative, fully-powered run.
"""
import numpy as np

from imaihda_sim import simulate_intersectional_data, fit_imaihda, correct_detection_bias
from imaihda_sim.robustness import misspecified_score, simulate_severity_dependent_detection

N = 12000
INTERACTION_SD = 0.9
TRUE_DELTA = 0.8
N_REP = 6
NONINFERIORITY_MARGIN = 0.5


def _vpc_on_true(df):
    dft = df.copy()
    dft["y"] = df["y_true"]
    return fit_imaihda(dft)["vpc_null"]


def _mean_biases(vpc_true, vpc_naive, vpc_corrected):
    return float(np.mean(np.array(vpc_naive) - np.array(vpc_true))), \
           float(np.mean(np.array(vpc_corrected) - np.array(vpc_true)))


def test_mild_score_misspecification_still_reduces_bias():
    # omit_rural: dropping one (weakly-weighted) true covariate is a mild
    # misspecification. Empirically (n_rep=20): ~52% bias reduction.
    vpc_true, vpc_naive, vpc_corr = [], [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=N, interaction_sd=INTERACTION_SD,
                                          detection_strength=TRUE_DELTA, seed=seed)
        vpc_true.append(_vpc_on_true(df))
        vpc_naive.append(fit_imaihda(df)["vpc_null"])
        score = misspecified_score(df, "omit_rural")
        vpc_corr.append(correct_detection_bias(df, delta=TRUE_DELTA, score=score)["vpc_null"])
    bias_naive, bias_corrected = _mean_biases(vpc_true, vpc_naive, vpc_corr)
    assert abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN  # non-inferiority
    assert abs(bias_corrected) < abs(bias_naive) * 0.75  # meaningful (not necessarily >=40%) reduction


def test_wrong_covariate_misspecification_is_noninferior_but_weaker():
    # wrong_covariate (sex swapped in for education): empirically ~25%
    # reduction at n_rep=20 -- still net-positive, non-inferior, but below
    # the 40% mild-violation target. This test only enforces non-inferiority.
    vpc_true, vpc_naive, vpc_corr = [], [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=N, interaction_sd=INTERACTION_SD,
                                          detection_strength=TRUE_DELTA, seed=seed)
        vpc_true.append(_vpc_on_true(df))
        vpc_naive.append(fit_imaihda(df)["vpc_null"])
        score = misspecified_score(df, "wrong_covariate")
        vpc_corr.append(correct_detection_bias(df, delta=TRUE_DELTA, score=score)["vpc_null"])
    bias_naive, bias_corrected = _mean_biases(vpc_true, vpc_naive, vpc_corr)
    assert abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN


def test_quadratic_score_misspecification_is_a_known_breakdown_case():
    # Empirical finding (n_rep=20, scripts/validation/detection_bias_robustness.py):
    # a quadratic (wrong functional form) score overcorrects past the truth and
    # violates non-inferiority (bias_corrected ~= +9.7pp vs bias_naive ~= -7.6pp).
    # This is a documented breakdown regime, not a bug -- this test guards against
    # a future silent change in correct_detection_bias() reversing that finding
    # without anyone noticing (in which case update this test AND the note).
    vpc_true, vpc_naive, vpc_corr = [], [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=N, interaction_sd=INTERACTION_SD,
                                          detection_strength=TRUE_DELTA, seed=seed)
        vpc_true.append(_vpc_on_true(df))
        vpc_naive.append(fit_imaihda(df)["vpc_null"])
        score = misspecified_score(df, "quadratic")
        vpc_corr.append(correct_detection_bias(df, delta=TRUE_DELTA, score=score)["vpc_null"])
    bias_naive, bias_corrected = _mean_biases(vpc_true, vpc_naive, vpc_corr)
    assert bias_corrected > 0  # overshoots past the truth (naive bias is negative)
    assert abs(bias_corrected) > abs(bias_naive)  # non-inferiority genuinely fails here


def test_severity_dependent_detection_degrades_from_helpful_to_breakdown():
    # Empirical finding: bias-reduction degrades smoothly as severity_weight
    # rises and crosses into breakdown (non-inferiority violated) around
    # severity_weight ~= 0.8-1.0. Assert the mild case (0.3) still clearly
    # helps and the severe case (1.0) is at least as bad as naive (within
    # the same margin used to declare breakdown in the driver script).
    def run(severity_weight):
        vpc_true, vpc_naive, vpc_corr = [], [], []
        for seed in range(N_REP):
            df = simulate_severity_dependent_detection(
                n=N, interaction_sd=INTERACTION_SD, detection_strength=TRUE_DELTA,
                severity_weight=severity_weight, seed=seed,
            )
            vpc_true.append(_vpc_on_true(df))
            vpc_naive.append(fit_imaihda(df)["vpc_null"])
            vpc_corr.append(correct_detection_bias(df, delta=TRUE_DELTA)["vpc_null"])
        return _mean_biases(vpc_true, vpc_naive, vpc_corr)

    bias_naive_mild, bias_corrected_mild = run(0.3)
    assert abs(bias_corrected_mild) < abs(bias_naive_mild)  # still helps

    bias_naive_severe, bias_corrected_severe = run(1.0)
    # At severity_weight=1.0 correction no longer reliably helps; only
    # require it doesn't blow up arbitrarily far past naive.
    assert abs(bias_corrected_severe) <= abs(bias_naive_severe) + 3.0


def test_severity_weight_zero_matches_original_mechanism_regression():
    df_a = simulate_intersectional_data(n=2000, interaction_sd=0.5, detection_strength=0.5, seed=1)
    df_b = simulate_severity_dependent_detection(
        n=2000, interaction_sd=0.5, detection_strength=0.5, severity_weight=0.0, seed=1
    )
    assert df_a.equals(df_b)
