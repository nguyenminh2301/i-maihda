"""PCV-focused robustness regression tests (previously deferred in
METHODS_NOTE_ROBUSTNESS.md Sec. 7). Calibrated to the empirical findings of
scripts/validation/pcv_robustness.py (n_rep=20) — which are largely honest
NEGATIVE findings, guarded here so future changes can't silently erase them:

1. Under pure sparsity the naive PCV is catastrophically inflated (~+45pp:
   var_main truncates to zero, so PCV reads ~100% "purely additive").
2. The detection correction's stratum-level WLS main-effects model degrades
   the PCV relative to naive even under CORRECT specification (+15.5 vs
   +3.7pp) — a real, documented limitation; a PCV-specific correction
   remains open.
3. Under strong outcome-dependent detection the sign flips: naive PCV bias
   grows large (+15.3pp) and the corrected PCV is closer to truth (+3.1pp).
"""
import numpy as np
import pandas as pd

from imaihda_sim import simulate_intersectional_data, fit_imaihda, correct_detection_bias
from imaihda_sim.detection import _aggregate_strata
from imaihda_sim.robustness import simulate_severity_dependent_detection

N_REP = 8


def _pop_truth_pcv(df):
    strata = _aggregate_strata(df)
    add = (
        -2.10
        + 0.20 * pd.to_numeric(strata["sex"]).to_numpy(float)
        + 0.35 * pd.to_numeric(strata["education"]).to_numpy(float)
        + 0.30 * pd.to_numeric(strata["wealth"]).to_numpy(float)
        + 0.25 * pd.to_numeric(strata["rural"]).to_numpy(float)
    )
    resid = df.groupby("stratum", observed=True)["true_stratum_residual"].first().to_numpy(float)
    v_null = float(np.var(add + resid, ddof=1))
    v_main = float(np.var(resid, ddof=1))
    return 100.0 * (v_null - v_main) / v_null


def test_sparse_pcv_catastrophically_inflated():
    biases = []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=3500, interaction_sd=0.9, sparse=True, seed=seed)
        biases.append(fit_imaihda(df)["pcv"] - _pop_truth_pcv(df))
    assert np.mean(biases) > 30.0  # PCV reads spuriously "additive" under sparsity


def test_detection_corrected_pcv_degrades_under_correct_specification():
    # Documented negative finding: the correction's WLS main-model inflates
    # corrected PCV past the naive value even at the true delta and score.
    naive, corr = [], []
    for seed in range(N_REP):
        df = simulate_intersectional_data(n=12000, interaction_sd=0.9,
                                          detection_strength=0.8, seed=seed)
        truth = _pop_truth_pcv(df)
        naive.append(fit_imaihda(df)["pcv"] - truth)
        corr.append(correct_detection_bias(df, delta=0.8)["pcv"] - truth)
    assert abs(np.mean(corr)) > abs(np.mean(naive))


def test_corrected_pcv_helps_under_strong_outcome_dependent_detection():
    # PCV is very seed-noisy (per-seed corrected biases range -12 to +14pp),
    # so this test uses the driver script's full n_rep=20 rather than the
    # thinner N_REP used above; at n_rep=8 the naive mean happens to sit near
    # zero and the comparison is unstable.
    naive, corr = [], []
    for seed in range(20):
        df = simulate_severity_dependent_detection(
            n=12000, interaction_sd=0.9, detection_strength=0.8,
            severity_weight=1.0, seed=seed)
        truth = _pop_truth_pcv(df)
        naive.append(fit_imaihda(df)["pcv"] - truth)
        corr.append(correct_detection_bias(df, delta=0.8)["pcv"] - truth)
    assert abs(np.mean(corr)) < abs(np.mean(naive))
