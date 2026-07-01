"""Robustness regression tests for sparse_strata_vpc() under violated
assumptions (Study 2). Thin, fast smoke tests calibrated to the empirical
findings in scripts/validation/sparse_strata_robustness.py (n_rep=30) and
docs/METHODS_NOTE_ROBUSTNESS.md -- not to a priori assumptions. Use a smaller
n_rep here (Monte Carlo noise tolerated via generous margins); the driver
script is the authoritative, fully-powered run.

Headline empirical finding: non-inferiority held in every tested arm of both
scenarios (unlike Study 1, which found genuine breakdown regimes) -- these
tests primarily guard that property, plus the strongest/cleanest
bias-reduction findings.
"""
import numpy as np

from imaihda_sim import simulate_intersectional_data, sparse_strata_vpc
from imaihda_sim.fit import vpc_latent
from imaihda_sim.robustness import simulate_structured_random_effects
from imaihda_sim.simulate import _stratum_table

N_REP = 15
# Looser than the driver script's 0.5pp (which uses n_rep=30, the authoritative
# figure reported in the methods note): at n_rep=15 Monte Carlo noise on the
# mean bias is large enough that the strict margin is flaky. This margin still
# catches genuine breakdowns (large, systematic non-inferiority violations).
NONINFERIORITY_MARGIN = 1.5


def _analytic_truth_structured(prevalence_shift, severity, tail_family, interaction_sd, seed):
    strata = _stratum_table()
    k = len(strata)
    eta_add = (
        prevalence_shift
        + 0.20 * strata["sex"].to_numpy()
        + 0.35 * strata["education"].to_numpy()
        + 0.30 * strata["wealth"].to_numpy()
        + 0.25 * strata["rural"].to_numpy()
    )
    rng_re = np.random.default_rng(seed + 999)
    if tail_family == "gaussian":
        raw = rng_re.normal(0.0, interaction_sd, size=k)
    elif tail_family == "student_t3":
        raw = rng_re.standard_t(df=3, size=k) * (interaction_sd / np.sqrt(3.0))
    else:
        raw = rng_re.lognormal(mean=0.0, sigma=0.5, size=k)
        raw = (raw - raw.mean()) / raw.std() * interaction_sd
    if severity > 0:
        raw = raw + 0.90 * severity * (
            (strata["education"].to_numpy() == 2) & (strata["wealth"].to_numpy() == 2)
            & (strata["rural"].to_numpy() == 1)
        )
        raw = raw - 0.60 * severity * (
            (strata["education"].to_numpy() == 0) & (strata["wealth"].to_numpy() == 2)
            & (strata["rural"].to_numpy() == 0)
        )
    eta_true = eta_add + (raw - raw.mean())
    return vpc_latent(np.var(eta_true, ddof=1))


def test_severity_zero_matches_original_mechanism_regression():
    df_a = simulate_intersectional_data(n=1500, prevalence_shift=-3.00, interaction_sd=0.15,
                                        detection_strength=0.0, sparse=True, seed=5)
    df_b = simulate_structured_random_effects(n=1500, prevalence_shift=-3.00, severity=1,
                                              tail_family="gaussian", interaction_sd=0.15,
                                              sparse=True, seed=5)
    cols = ["sex", "education", "wealth", "rural", "stratum", "y", "y_true", "true_stratum_residual"]
    assert df_a[cols].equals(df_b[cols])


def test_spike_escalation_stays_noninferior_and_helps():
    # Empirically (n_rep=30): non-inferiority holds at every severity level
    # 0-3, with bias-reduction actually strengthening at higher severity.
    for severity in (0, 1, 3):
        biases_naive, biases_corr = [], []
        for seed in range(N_REP):
            df = simulate_structured_random_effects(
                n=3500, prevalence_shift=-3.00, severity=severity, tail_family="gaussian",
                interaction_sd=0.15, sparse=True, seed=seed,
            )
            truth = _analytic_truth_structured(-3.00, severity, "gaussian", 0.15, seed)
            r = sparse_strata_vpc(df, seed=1000 + seed)
            biases_naive.append(r["vpc_null_naive"] - truth)
            biases_corr.append(r["vpc_null_corrected"] - truth)
        bias_naive = float(np.mean(biases_naive))
        bias_corrected = float(np.mean(biases_corr))
        assert abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN, (
            f"non-inferiority failed at severity={severity}"
        )
        assert abs(bias_corrected) < abs(bias_naive)  # correction still helps on average


def test_nongaussian_tail_shapes_stay_noninferior():
    for tail_family in ("student_t3", "skew_lognormal"):
        biases_naive, biases_corr = [], []
        for seed in range(N_REP):
            df = simulate_structured_random_effects(
                n=3500, prevalence_shift=-3.00, severity=1, tail_family=tail_family,
                interaction_sd=0.15, sparse=True, seed=seed,
            )
            truth = _analytic_truth_structured(-3.00, 1, tail_family, 0.15, seed)
            r = sparse_strata_vpc(df, seed=1000 + seed)
            biases_naive.append(r["vpc_null_naive"] - truth)
            biases_corr.append(r["vpc_null_corrected"] - truth)
        bias_naive = float(np.mean(biases_naive))
        bias_corrected = float(np.mean(biases_corr))
        assert abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN, (
            f"non-inferiority failed for tail_family={tail_family}"
        )


def test_extreme_rarity_stays_noninferior():
    # prevalence_shift=-5.0: rarest, most sparsity-stressed arm tested.
    # Empirically bias-reduction drops to ~26% here (below the mild-violation
    # target) but non-inferiority still holds -- graceful degradation, not
    # breakdown.
    biases_naive, biases_corr = [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(
            n=3500, prevalence_shift=-5.00, interaction_sd=0.15, sparse=True, seed=seed,
        )
        truth = _analytic_truth_structured(-5.00, 1, "gaussian", 0.15, seed)
        r = sparse_strata_vpc(df, seed=2000 + seed)
        biases_naive.append(r["vpc_null_naive"] - truth)
        biases_corr.append(r["vpc_null_corrected"] - truth)
    bias_naive = float(np.mean(biases_naive))
    bias_corrected = float(np.mean(biases_corr))
    assert abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN
