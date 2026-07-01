import numpy as np
import pandas as pd

from imaihda_sim import (
    simulate_intersectional_data,
    fit_imaihda,
    correct_detection_bias,
    vpc_detection_bounds,
    detection_tipping_point,
)


def _vpc_on_true(df):
    """VPC computed on the true (pre-detection) outcome — the recovery target."""
    dft = df.copy()
    dft["y"] = df["y_true"]
    return fit_imaihda(dft)["vpc_null"]


def test_delta_zero_is_identity():
    # With no correction, the sensitivity analysis must return exactly the
    # observed VPC/PCV (the anchor point of the sweep).
    df = simulate_intersectional_data(n=6000, interaction_sd=0.9,
                                      detection_strength=0.8, seed=7)
    obs = fit_imaihda(df)
    corr = correct_detection_bias(df, delta=0.0)
    # Null-model VPC (the headline discriminatory-accuracy quantity) is an exact
    # identity at delta = 0.
    assert abs(corr["vpc_null"] - obs["vpc_null"]) < 1e-9
    # Main-effects VPC uses a stratum-level WLS fit (the correction path has no
    # individual-level data to refit a GLM), so it matches the individual-level
    # GLM to within small aggregation error rather than exactly.
    assert abs(corr["vpc_main"] - obs["vpc_main"]) < 0.5


def test_correction_recovers_true_vpc():
    # When delta matches the generating strength, correcting the observed data
    # should land much closer to the true VPC than the (biased) observed VPC.
    delta_true = 0.8
    df = simulate_intersectional_data(n=12000, interaction_sd=0.9,
                                      detection_strength=delta_true, seed=7)
    vpc_true = _vpc_on_true(df)
    vpc_obs = fit_imaihda(df)["vpc_null"]
    vpc_corr = correct_detection_bias(df, delta=delta_true)["vpc_null"]

    # Detection bias masks between-stratum variance -> observed < true.
    assert vpc_obs < vpc_true
    # Correction removes the differential distortion.
    assert abs(vpc_corr - vpc_true) < abs(vpc_obs - vpc_true)
    assert abs(vpc_corr - vpc_true) < 2.0  # within ~2 percentage points


def test_bounds_contain_true_vpc():
    delta_true = 0.8
    df = simulate_intersectional_data(n=12000, interaction_sd=0.9,
                                      detection_strength=delta_true, seed=7)
    vpc_true = _vpc_on_true(df)
    bounds = vpc_detection_bounds(df, delta_max=1.2, n_grid=25)

    assert len(bounds) == 25
    assert bounds["delta"].iloc[0] == 0.0
    # delta = 0 row equals the observed VPC.
    assert abs(bounds["vpc_null"].iloc[0] - fit_imaihda(df)["vpc_null"]) < 1e-9
    # The plausible range brackets the true VPC.
    assert bounds["vpc_null"].min() <= vpc_true <= bounds["vpc_null"].max()
    # Implied true prevalence rises monotonically as assumed under-detection grows.
    assert np.all(np.diff(bounds["implied_true_prevalence"].to_numpy()) >= -1e-9)


def test_tipping_point_finds_crossing_and_is_monotone():
    df = simulate_intersectional_data(n=12000, interaction_sd=0.9,
                                      detection_strength=0.8, seed=7)
    vpc_obs = fit_imaihda(df)["vpc_null"]
    # A threshold above the observed VPC is reached at some positive delta,
    # since correcting for under-detection raises the VPC.
    target = vpc_obs + 2.0
    tp = detection_tipping_point(df, threshold=target, quantity="vpc_null",
                                 delta_max=3.0)
    assert not np.isnan(tp)
    assert 0.0 < tp <= 3.0
    # At the tipping point the corrected quantity equals the threshold.
    assert abs(correct_detection_bias(df, tp)["vpc_null"] - target) < 1e-3

    # A threshold that is never reached within [0, delta_max] returns nan.
    assert np.isnan(detection_tipping_point(df, threshold=vpc_obs - 5.0,
                                            quantity="vpc_null", delta_max=1.0))


def test_score_length_validation():
    df = simulate_intersectional_data(n=2000, seed=1)
    try:
        correct_detection_bias(df, delta=0.5, score=[1.0, 2.0])  # wrong length
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for mismatched score length")
