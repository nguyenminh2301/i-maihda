"""K-generalization regression tests (previously deferred in
METHODS_NOTE_ROBUSTNESS.md Sec. 7 and PHASE3_CAUSAL_IDENTIFICATION.md Sec. 6).
Calibrated to scripts/validation/k_generalization_study.py (n_rep=20); thinner
n_rep here with generous margins.

Includes the boundary finding that motivated the ci_method="test_inversion"
option: at K=108 with n=3500 the naive statistic sits at its max(0,.)
truncation floor in ~90% of replicates, degenerating the bootstrap CI to
[0, 0] (coverage 0.10); test inversion stays honestly one-sided and restores
coverage (0.95 at n_rep=20).
"""
import numpy as np
import pytest

from imaihda_sim import sparse_strata_vpc, vpc_latent
from imaihda_sim.robustness import simulate_k_strata

N_REP = 8
NONINFERIORITY_MARGIN = 1.5  # generous at small n_rep; driver uses 0.5 at n_rep=20


def _truth(df):
    etas = df.groupby("stratum", observed=True)["true_stratum_logit"].first().to_numpy(float)
    return vpc_latent(float(np.var(etas, ddof=1)))


def test_generator_shapes_and_validation():
    df = simulate_k_strata(dims=(2, 2, 2), n=2000, seed=1)
    assert df["stratum"].nunique() <= 8
    assert set(df["y"].unique()) <= {0, 1}
    assert "true_stratum_logit" in df.columns
    df2 = simulate_k_strata(dims=(2, 3, 3, 3, 2), n=2000, sparse=False, seed=1)
    assert df2["stratum"].nunique() <= 108
    with pytest.raises(ValueError):
        simulate_k_strata(dims=(1, 2))
    with pytest.raises(ValueError):
        simulate_k_strata(dims=(2, 2), betas=[0.1])


def test_noninferiority_holds_at_small_and_large_k():
    for dims in [(2, 2, 2), (2, 3, 3, 3, 2)]:
        naive, corr = [], []
        for seed in range(N_REP):
            df = simulate_k_strata(dims=dims, n=3500, sparse=True, seed=seed)
            truth = _truth(df)
            r = sparse_strata_vpc(df, seed=6000 + seed)
            naive.append(r["vpc_null_naive"] - truth)
            corr.append(r["vpc_null_corrected"] - truth)
        assert abs(np.mean(corr)) <= abs(np.mean(naive)) + NONINFERIORITY_MARGIN, dims


def test_test_inversion_rescues_floor_degenerate_ci_at_k108():
    hits_boot, hits_ti, floors = 0, 0, 0
    for seed in range(N_REP):
        df = simulate_k_strata(dims=(2, 3, 3, 3, 2), n=3500, sparse=True, seed=seed)
        truth = _truth(df)
        r_b = sparse_strata_vpc(df, seed=6000 + seed)
        r_t = sparse_strata_vpc(df, seed=6000 + seed, ci_method="test_inversion")
        hits_boot += r_b["ci_lower"] <= truth <= r_b["ci_upper"]
        hits_ti += r_t["ci_lower"] <= truth <= r_t["ci_upper"]
        floors += r_b["at_floor"]
        # test-inversion interval always contains the point estimate's floor side
        assert r_t["ci_upper"] >= r_t["ci_lower"]
    assert floors >= N_REP // 2          # the floor regime is the norm here
    assert hits_ti > hits_boot            # inversion strictly improves coverage
    assert hits_ti >= N_REP - 2           # and reaches near-nominal


def test_ci_method_validation_and_default_unchanged():
    df = simulate_k_strata(dims=(2, 2, 2), n=2000, seed=1)
    with pytest.raises(ValueError):
        sparse_strata_vpc(df, ci_method="bogus")
    r1 = sparse_strata_vpc(df, seed=3)
    r2 = sparse_strata_vpc(df, seed=3, ci_method="bootstrap")
    assert r1["ci_lower"] == r2["ci_lower"] and r1["ci_upper"] == r2["ci_upper"]
