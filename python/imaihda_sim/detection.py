"""Detection-bias sensitivity analysis for intersectional MAIHDA.

Quantitative bias analysis (QBA) for the between-stratum VPC and PCV under
SES-patterned under-detection of the outcome. Unlike CRAN ``MAIHDA`` — which
assumes the outcome is measured equally well across strata — these tools ask:
*if disadvantaged strata under-record true cases, how far could the observed
VPC/PCV be from the truth?*

The observed outcome is ``y = y_true * detected``. Detection removes only true
positives, so within a stratum ``j`` the observed prevalence deflates by the
stratum detection probability ``d_j``:

    E[p_obs_j] = p_true_j * d_j   =>   p_true_j = p_obs_j / d_j(delta)

The analyst does not know the detection strength, so ``delta`` is a sensitivity
parameter swept over a plausible range. ``delta = 0`` applies no correction and
reproduces the observed VPC/PCV exactly. As ``delta`` approaches the true
generating strength, the correction removes the SES-patterned *differential*
distortion and moves the estimate back toward the true VPC (see
``tests/test_detection_correction.py``). Only the differential pattern is
identifiable from observed data; a uniform ascertainment level that hits every
stratum equally is not, and is left as an explicit assumption.

These functions operate entirely at the stratum level (aggregated counts), so
they are fast and need no individual-level refit.
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from .fit import (
    LOGISTIC_L1_VARIANCE,
    _aggregate_strata,
    _between_stratum_variance,
    pcv,
    vpc_latent,
)

__all__ = [
    "correct_detection_bias",
    "vpc_detection_bounds",
    "detection_tipping_point",
]

# Pattern of SES-patterned under-detection, matching the simulator's
# detection model (simulate.py): logit(P(detect)) = baseline - delta * score,
# with score = education + wealth + 0.4 * rural.
_RURAL_WEIGHT = 0.40


def _logit_inv(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-x))


def _numeric(col) -> np.ndarray:
    """Coerce a (possibly categorical) stratum column to a float array."""
    return pd.to_numeric(pd.Series(np.asarray(col)), errors="coerce").to_numpy(dtype=float)


def _stratum_score(strata: pd.DataFrame, score: Optional[Sequence[float]]) -> np.ndarray:
    if score is not None:
        s = np.asarray(score, dtype=float)
        if s.shape[0] != len(strata):
            raise ValueError(
                f"score has length {s.shape[0]} but there are {len(strata)} strata"
            )
        return s
    return (
        _numeric(strata["education"])
        + _numeric(strata["wealth"])
        + _RURAL_WEIGHT * _numeric(strata["rural"])
    )


def _wls_main_pred(strata: pd.DataFrame, logit: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted-least-squares additive fit of stratum logits.

    Mirrors the individual-level additive GLM used by ``fit_imaihda`` but
    operates on corrected stratum logits, which is all the QBA needs.
    """
    cols = []
    for name in ("sex", "education", "wealth", "rural"):
        vals = _numeric(strata[name])
        levels = np.unique(vals)
        # Dummy-code with the lowest level as reference (drop first).
        for lv in levels[1:]:
            cols.append((vals == lv).astype(float))
    X = np.column_stack([np.ones(len(strata))] + cols) if cols else np.ones((len(strata), 1))

    sw = np.sqrt(weights)
    Xw = X * sw[:, None]
    yw = logit * sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return X @ beta


def correct_detection_bias(
    df: pd.DataFrame,
    delta: float,
    baseline_logit: float = 2.0,
    score: Optional[Sequence[float]] = None,
) -> Dict[str, float]:
    """Correct VPC/PCV for a hypothesized SES-patterned under-detection strength.

    Parameters
    ----------
    df : DataFrame
        Individual-level data with columns ``stratum, sex, education, wealth,
        rural, y`` (as produced by :func:`simulate_intersectional_data`).
    delta : float
        Sensitivity strength of SES-patterned under-detection. ``0`` applies no
        correction (returns the observed VPC/PCV).
    baseline_logit : float, default 2.0
        Controls how steeply detection declines with disadvantage. Detection is
        taken relative to the most-advantaged stratum (``score = 0``), so
        ``delta = 0`` is the identity regardless of this value. The default
        matches the simulator's detection curvature.
    score : sequence of float, optional
        Per-stratum disadvantage score (one value per stratum, in the order
        returned by the internal aggregation). Defaults to
        ``education + wealth + 0.4 * rural``.

    Returns
    -------
    dict
        ``delta, baseline_logit, var_null, vpc_null, var_main, vpc_main, pcv,
        observed_prevalence, implied_true_prevalence, min_detection_probability``.
    """
    if delta < 0:
        raise ValueError("delta must be non-negative")

    strata = _aggregate_strata(df)
    e_obs = strata["events"].to_numpy(dtype=float)
    n = strata["n"].to_numpy(dtype=float)

    s = _stratum_score(strata, score)
    # Detection probability *relative* to the most-advantaged stratum (score = 0).
    # Normalizing this way makes delta = 0 the identity (only the SES-patterned
    # *differential* under-detection is identifiable from observed data; a uniform
    # ascertainment level that hits every stratum equally is not). baseline_logit
    # controls how steeply detection declines with disadvantage.
    d0 = float(_logit_inv(baseline_logit))
    d = np.asarray(_logit_inv(baseline_logit - delta * s), dtype=float) / d0
    d = np.clip(d, 1e-6, 1.0)

    # Invert the detection process: reconstruct implied true events, capped at n.
    e_true = np.minimum(n, e_obs / d)

    # Recompute Laplace-smoothed logits and delta-method sampling variance
    # from the corrected counts (same formulas as _aggregate_strata).
    p_corr = (e_true + 0.5) / (n + 1.0)
    logit_corr = np.log(p_corr / (1.0 - p_corr))
    sampling_var = 1.0 / (e_true + 0.5) + 1.0 / (n - e_true + 0.5)
    weights = 1.0 / sampling_var

    p_overall = float(np.clip(e_true.sum() / n.sum(), 1e-6, 1.0 - 1e-6))
    null_pred = np.repeat(math.log(p_overall / (1.0 - p_overall)), len(strata))
    var0 = _between_stratum_variance(logit_corr, null_pred, sampling_var, weights)

    main_pred = _wls_main_pred(strata, logit_corr, weights)
    var1 = _between_stratum_variance(logit_corr, main_pred, sampling_var, weights)

    return {
        "delta": float(delta),
        "baseline_logit": float(baseline_logit),
        "var_null": var0,
        "vpc_null": vpc_latent(var0),
        "var_main": var1,
        "vpc_main": vpc_latent(var1),
        "pcv": pcv(var0, var1),
        "observed_prevalence": float(e_obs.sum() / n.sum()),
        "implied_true_prevalence": p_overall,
        "min_relative_detection": float(d.min()),
    }


def vpc_detection_bounds(
    df: pd.DataFrame,
    delta_max: float = 1.0,
    n_grid: int = 21,
    baseline_logit: float = 2.0,
    score: Optional[Sequence[float]] = None,
) -> pd.DataFrame:
    """Sweep the under-detection strength and tabulate corrected VPC/PCV.

    Returns one row per ``delta`` in ``linspace(0, delta_max, n_grid)``. The
    range of ``vpc_null`` across rows is the sensitivity bound on the true VPC.
    """
    if delta_max < 0:
        raise ValueError("delta_max must be non-negative")
    if n_grid < 2:
        raise ValueError("n_grid must be at least 2")

    deltas = np.linspace(0.0, delta_max, n_grid)
    rows = [
        correct_detection_bias(df, float(d), baseline_logit=baseline_logit, score=score)
        for d in deltas
    ]
    return pd.DataFrame(rows)


def detection_tipping_point(
    df: pd.DataFrame,
    threshold: float,
    quantity: str = "vpc_null",
    delta_max: float = 3.0,
    baseline_logit: float = 2.0,
    score: Optional[Sequence[float]] = None,
) -> float:
    """Minimum under-detection strength that moves ``quantity`` to ``threshold``.

    An E-value analogue for MAIHDA: the smallest ``delta`` at which the
    corrected quantity (default null-model VPC) first crosses ``threshold``.
    Returns ``nan`` if the threshold is not crossed within ``[0, delta_max]``.
    """

    def f(d: float) -> float:
        return correct_detection_bias(
            df, d, baseline_logit=baseline_logit, score=score
        )[quantity] - threshold

    f0 = f(0.0)
    if f0 == 0.0:
        return 0.0
    f_max = f(delta_max)
    if np.isnan(f0) or np.isnan(f_max) or np.sign(f0) == np.sign(f_max):
        return float("nan")
    return float(brentq(f, 0.0, delta_max))
