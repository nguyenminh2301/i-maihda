"""Composition-grid regression tests (generalizing the single-scenario
capstone finding, previously flagged in METHODS_NOTE_ROBUSTNESS.md Sec. 7).
Calibrated to scripts/validation/composition_grid.py (n_rep=12, full 3x2
grid); thinner slice here for speed.

Grid findings guarded: (1) the sequential composition's bias is *unstable in
delta* — it climbs monotonically from negative through zero to strongly
positive as detection strengthens, so its occasional accuracy is two errors
cancelling, not reliability; (2) the joint estimator stays stable (|bias|
within a few pp) across regimes and ranked best-or-second in 6/6 grid cells.
"""
import numpy as np
import pandas as pd

from imaihda_sim import simulate_intersectional_data, fit_imaihda, joint_calibrated_vpc, vpc_latent
from imaihda_sim.detection import _aggregate_strata, _stratum_score, _logit_inv
from imaihda_sim.uncertainty import _var_null_from_counts, _build_calibration_curve

N_REP = 6


def _pop_truth(df):
    strata = _aggregate_strata(df)
    add = (
        -2.10
        + 0.20 * pd.to_numeric(strata["sex"]).to_numpy(float)
        + 0.35 * pd.to_numeric(strata["education"]).to_numpy(float)
        + 0.30 * pd.to_numeric(strata["wealth"]).to_numpy(float)
        + 0.25 * pd.to_numeric(strata["rural"]).to_numpy(float)
    )
    resid = df.groupby("stratum", observed=True)["true_stratum_residual"].first().to_numpy(float)
    return vpc_latent(float(np.var(add + resid, ddof=1)))


def _composed(df, delta, seed):
    strata = _aggregate_strata(df)
    e_obs = strata["events"].to_numpy(float)
    n_j = strata["n"].to_numpy(float)
    s = _stratum_score(strata, None)
    d0 = float(_logit_inv(2.0))
    d_rel = np.clip(np.asarray(_logit_inv(2.0 - delta * s), float) / d0, 1e-6, 1.0)
    e_corr = np.minimum(n_j, e_obs / d_rel)
    p0 = float(np.clip(e_corr.sum() / n_j.sum(), 1e-6, 1 - 1e-6))
    var_naive = float(_var_null_from_counts(e_corr[None, :], n_j[None, :], p0)[0])
    rng = np.random.default_rng(seed)
    grid, means = _build_calibration_curve(n_j.astype(np.int64), p0, 2.0, 25, 200, rng)
    if var_naive <= means[0]:
        v = 0.0
    elif var_naive >= means[-1]:
        v = float(grid[-1])
    else:
        v = float(np.interp(var_naive, means, grid))
    return vpc_latent(v)


def _mean_biases(delta):
    composed, joint = [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=3500, interaction_sd=0.9,
                                          detection_strength=delta, sparse=True, seed=seed)
        truth = _pop_truth(df)
        composed.append(_composed(df, delta, seed=8000 + seed) - truth)
        joint.append(joint_calibrated_vpc(df, delta=delta, seed=9000 + seed)["vpc_corrected"] - truth)
    return float(np.mean(composed)), float(np.mean(joint))


def test_composed_bias_climbs_with_delta_while_joint_stays_stable():
    comp_lo, joint_lo = _mean_biases(0.4)
    comp_hi, joint_hi = _mean_biases(1.2)
    # Composed drifts upward strongly as delta grows (driver: -0.4 -> +7.3).
    assert comp_hi - comp_lo > 3.0
    # Joint stays bounded in both regimes (driver: within +-3pp in all cells).
    assert abs(joint_lo) < 6.0
    assert abs(joint_hi) < 6.0


def test_joint_not_worse_than_naive_anywhere_on_slice():
    for delta in (0.4, 1.2):
        naives, joints = [], []
        for seed in range(N_REP):
            df = simulate_intersectional_data(n=3500, interaction_sd=0.9,
                                              detection_strength=delta, sparse=True, seed=seed)
            truth = _pop_truth(df)
            naives.append(fit_imaihda(df)["vpc_null"] - truth)
            joints.append(joint_calibrated_vpc(df, delta=delta, seed=9000 + seed)["vpc_corrected"] - truth)
        assert abs(np.mean(joints)) < abs(np.mean(naives))
