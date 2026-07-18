"""Tests for the fit_imaihda(weighting=) option added to close the documented
Python/R variance-formula discrepancy (see uncertainty.py module docstring):
the default stays precision-weighted for benchmark continuity, while
weighting="unweighted" matches R's v0.2.1-corrected formula.
"""
import numpy as np
import pytest

from imaihda_sim import simulate_intersectional_data, fit_imaihda
from imaihda_sim.detection import _aggregate_strata
from imaihda_sim.uncertainty import _var_null_from_counts


def test_default_weighting_unchanged():
    df = simulate_intersectional_data(n=3000, interaction_sd=0.9, seed=3)
    r_default = fit_imaihda(df)
    r_explicit = fit_imaihda(df, weighting="precision")
    assert r_default["var_null"] == r_explicit["var_null"]
    assert r_default["var_main"] == r_explicit["var_main"]
    assert r_default["weighting"] == "precision"


def test_unweighted_matches_uncertainty_module_formula():
    # weighting="unweighted" must reproduce exactly the null-model variance
    # formula used by sparse_strata_vpc's calibration (uncertainty.py), which
    # in turn matches R's diagnostics.R.
    df = simulate_intersectional_data(n=3000, interaction_sd=0.9, seed=3)
    r = fit_imaihda(df, weighting="unweighted")

    strata = _aggregate_strata(df)
    e = strata["events"].to_numpy(float)
    n = strata["n"].to_numpy(float)
    p0 = float(np.clip(e.sum() / n.sum(), 1e-6, 1 - 1e-6))
    expected = float(_var_null_from_counts(e[None, :], n[None, :], p0)[0])
    assert abs(r["var_null"] - expected) < 1e-12
    assert r["weighting"] == "unweighted"


def test_formulas_differ_on_sparse_data():
    # The two formulas are genuinely different estimators; on sparse data they
    # should not coincide (guards against a silent no-op implementation).
    df = simulate_intersectional_data(n=3500, interaction_sd=0.9, sparse=True, seed=7)
    v_w = fit_imaihda(df, weighting="precision")["var_null"]
    v_u = fit_imaihda(df, weighting="unweighted")["var_null"]
    assert abs(v_w - v_u) > 1e-6


def test_invalid_weighting_rejected():
    df = simulate_intersectional_data(n=1000, seed=1)
    with pytest.raises(ValueError):
        fit_imaihda(df, weighting="bogus")
