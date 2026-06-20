import numpy as np
from imaihda_sim import simulate_intersectional_data, fit_imaihda, vpc_latent, pcv


def test_vpc_pcv_formulas():
    assert round(vpc_latent(0), 8) == 0
    assert 0 < vpc_latent(0.5) < 100
    assert pcv(1.0, 0.25) == 75.0


def test_simulation_has_expected_structure():
    df = simulate_intersectional_data(n=1000, seed=123)
    assert df["stratum"].nunique() == 36
    assert set(df["y"].unique()).issubset({0, 1})
    assert 0 < df["y"].mean() < 1


def test_additive_scenario_fits_and_is_additive_dominant():
    df = simulate_intersectional_data(n=1800, seed=123)
    res = fit_imaihda(df)
    assert res["pcv"] > 70
    assert res["vpc_main"] < res["vpc_null"]


def test_detection_bias_reduces_observed_prevalence():
    d0 = simulate_intersectional_data(n=1200, seed=123, detection_strength=0)
    d1 = simulate_intersectional_data(n=1200, seed=123, detection_strength=1.0)
    assert d1["y"].mean() < d0["y"].mean()
    assert d1["y_true"].mean() == d0["y_true"].mean()
