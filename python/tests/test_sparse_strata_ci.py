import numpy as np
import pytest

from imaihda_sim import simulate_intersectional_data, sparse_strata_vpc
from imaihda_sim.simulate import _stratum_table
from imaihda_sim.fit import vpc_latent

INTERACTION_SD = 0.15
INTERACTION_SEED = 999  # matches simulate_intersectional_data(seed=999)


def _analytic_truth(interaction_sd=INTERACTION_SD, seed=INTERACTION_SEED):
    """Exact population VPC of the true stratum logits, computed directly from
    the generative formula in simulate.py -- no simulation, no estimator bias.
    This is the ground truth that vpc_null (of any estimator) targets."""
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


TRUTH = _analytic_truth()  # ~10.06, computed once for the module


def test_naive_underestimates_under_genuine_sparsity():
    # Scenario E-style allocation (sparse=True): most strata get very few
    # individuals. The naive fast estimator should be substantially biased
    # downward relative to the analytic ground truth.
    naive = []
    for seed in range(8):
        df = simulate_intersectional_data(n=3500, interaction_sd=INTERACTION_SD,
                                          sparse=True, seed=seed)
        naive.append(sparse_strata_vpc(df, seed=100 + seed)["vpc_null_naive"])
    assert np.mean(naive) < TRUTH - 2.0


def test_correction_reduces_average_bias():
    naive, corrected = [], []
    for seed in range(8):
        df = simulate_intersectional_data(n=3500, interaction_sd=INTERACTION_SD,
                                          sparse=True, seed=seed)
        r = sparse_strata_vpc(df, seed=100 + seed)
        naive.append(r["vpc_null_naive"])
        corrected.append(r["vpc_null_corrected"])
    naive_gap = abs(np.mean(naive) - TRUTH)
    corrected_gap = abs(np.mean(corrected) - TRUTH)
    assert corrected_gap < naive_gap


def test_ci_well_formed_and_sparse_flag():
    df = simulate_intersectional_data(n=3500, interaction_sd=INTERACTION_SD,
                                      sparse=True, seed=1)
    r = sparse_strata_vpc(df, seed=7)
    assert r["ci_lower"] <= r["ci_upper"]
    assert 0.0 <= r["ci_lower"] <= 100.0
    assert 0.0 <= r["ci_upper"] <= 100.0
    assert r["sparse"] is True
    assert r["min_stratum_n"] < 20


def test_ci_covers_truth_at_reasonable_rate():
    # Loose coverage check (not a strict 95% test, to avoid flakiness with a
    # small number of trials): the interval should contain the analytic truth
    # most of the time under this sparse regime.
    hits = 0
    trials = 10
    for seed in range(trials):
        df = simulate_intersectional_data(n=3500, interaction_sd=INTERACTION_SD,
                                          sparse=True, seed=seed)
        r = sparse_strata_vpc(df, seed=200 + seed)
        if r["ci_lower"] <= TRUTH <= r["ci_upper"]:
            hits += 1
    assert hits >= 6  # >= 60% of nominal 95%, generous margin for n=10 trials


def test_dense_well_powered_data_needs_little_correction():
    # With large, evenly-allocated strata the naive estimator is already close
    # to unbiased, so calibration should not move the estimate much.
    df = simulate_intersectional_data(n=60000, interaction_sd=INTERACTION_SD, seed=3)
    r = sparse_strata_vpc(df, seed=42)
    assert r["sparse"] is False
    assert abs(r["vpc_null_corrected"] - r["vpc_null_naive"]) < 1.5


def test_reproducible_with_same_seed():
    df = simulate_intersectional_data(n=1500, interaction_sd=INTERACTION_SD,
                                      sparse=True, seed=5)
    r1 = sparse_strata_vpc(df, seed=11)
    r2 = sparse_strata_vpc(df, seed=11)
    assert r1 == r2


def test_input_validation():
    df = simulate_intersectional_data(n=2000, seed=1)
    with pytest.raises(ValueError):
        sparse_strata_vpc(df, level=1.5)
    with pytest.raises(ValueError):
        sparse_strata_vpc(df, level=0.0)

    tiny = df[df["stratum"] == df["stratum"].iloc[0]]
    with pytest.raises(ValueError):
        sparse_strata_vpc(tiny)
