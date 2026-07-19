"""Assumption-violation data-generating processes for the robustness/simulation
study of `detection.py` and `uncertainty.py`.

Both `correct_detection_bias()` and `sparse_strata_vpc()` were validated under
the exact assumptions their generative model encodes (the analyst's detection
`score=` formula matches the true mechanism; true random effects are i.i.d.
Gaussian). This module generates data under DEVIATIONS from those assumptions,
so the two methods' behavior can be stress-tested and honestly characterized
-- including the regime where they break down -- rather than only validated
on the easy case.

Every generator here returns the same column contract as
`simulate_intersectional_data()` (`stratum, sex, education, wealth, rural, y,
y_true, detection_probability, true_stratum_residual`), so `fit_imaihda()`,
`correct_detection_bias()`, and `sparse_strata_vpc()` all consume its output
unmodified. `simulate_intersectional_data()` itself is left untouched (it is
relied on by the existing, already-validated test suite).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from .fit import _aggregate_strata
from .simulate import _stratum_table

__all__ = [
    "simulate_severity_dependent_detection",
    "simulate_structured_random_effects",
    "misspecified_score",
    "simulate_k_strata",
    "simulate_selective_attrition",
]

# Same asymmetric spike ratio as simulate.py's baseline (+0.90 / -0.60),
# scaled by `severity` so severity=1 reproduces the original mechanism exactly.
_SPIKE_HIGH_UNIT = 0.90
_SPIKE_LOW_UNIT = -0.60


def _logit_inv(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def _numeric(col) -> np.ndarray:
    return pd.to_numeric(pd.Series(np.asarray(col)), errors="coerce").to_numpy(dtype=float)


def simulate_severity_dependent_detection(
    n: int = 6000,
    prevalence_shift: float = -2.10,
    interaction_sd: float = 0.0,
    detection_strength: float = 0.0,
    severity_weight: float = 0.0,
    sparse: bool = False,
    seed: int = 42,
) -> pd.DataFrame:
    """Like `simulate_intersectional_data`, but detection can also depend
    directly on the individual's TRUE outcome status (classic verification /
    outcome-dependent ascertainment bias), not only on stratum covariates:

        logit(P(detect)) = 2.0 - detection_strength*(education + wealth
                            + 0.4*rural) - severity_weight * y_true

    `severity_weight=0` reproduces `simulate_intersectional_data`'s detection
    mechanism exactly (regression-checked in
    `python/tests/test_detection_robustness.py`). `severity_weight > 0` means
    true cases are additionally under-detected *regardless of SES* -- this is
    not identifiable from stratum-level covariates alone, so no `score=`
    choice passed to `correct_detection_bias()` can fully correct for it. It
    is included specifically to characterize where the detection-bias
    correction's identifying assumption (differential detection acts only
    through observed covariates) breaks down.
    """
    rng = np.random.default_rng(seed)
    strata = _stratum_table()
    k = len(strata)

    if sparse:
        weights = rng.gamma(shape=0.35, scale=1.0, size=k)
        weights = weights / weights.sum()
    else:
        weights = np.repeat(1.0 / k, k)

    stratum_id = rng.choice(k, size=n, p=weights)
    x = strata.iloc[stratum_id].reset_index(drop=True)

    linear_predictor = (
        prevalence_shift
        + 0.20 * x["sex"].to_numpy()
        + 0.35 * x["education"].to_numpy()
        + 0.30 * x["wealth"].to_numpy()
        + 0.25 * x["rural"].to_numpy()
    )

    true_stratum_residual = np.zeros(k)
    if interaction_sd > 0:
        rng_re = np.random.default_rng(seed + 999)
        raw = rng_re.normal(loc=0.0, scale=interaction_sd, size=k)
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
        true_stratum_residual = raw - raw.mean()
        linear_predictor = linear_predictor + true_stratum_residual[stratum_id]

    p_true = _logit_inv(linear_predictor)
    y_true = rng.binomial(1, p_true)

    if detection_strength > 0 or severity_weight > 0:
        detect_logit = (
            2.0
            - detection_strength * x["education"].to_numpy()
            - detection_strength * x["wealth"].to_numpy()
            - 0.40 * detection_strength * x["rural"].to_numpy()
            - severity_weight * y_true
        )
        detect_p = _logit_inv(detect_logit)
        detected = rng.binomial(1, detect_p)
        y_observed = y_true * detected
    else:
        detect_p = np.ones(n)
        y_observed = y_true

    out = x.copy()
    out["stratum"] = pd.Categorical(stratum_id.astype(str))
    for col in ["sex", "education", "wealth", "rural"]:
        out[col] = pd.Categorical(out[col])
    out["y"] = y_observed.astype(int)
    out["y_true"] = y_true.astype(int)
    out["detection_probability"] = detect_p
    out["true_stratum_residual"] = true_stratum_residual[stratum_id]
    return out


def simulate_structured_random_effects(
    n: int = 3500,
    prevalence_shift: float = -3.00,
    severity: int = 1,
    tail_family: str = "gaussian",
    interaction_sd: float = 0.15,
    sparse: bool = True,
    seed: int = 42,
) -> pd.DataFrame:
    """Like `simulate_intersectional_data(interaction_sd=..., sparse=...)`,
    but the true residual stratum effects can deviate from i.i.d. Gaussian in
    two independent ways, to stress-test `sparse_strata_vpc()`'s Gaussian
    calibration model:

    - `tail_family`: shape of the base (unspiked) noise, standardized to the
      same dispersion (`interaction_sd`) regardless of family so only shape
      varies: {"gaussian" (matches the original mechanism), "student_t3"
      (heavy-tailed), "skew_lognormal" (right-skewed)}.
    - `severity`: magnitude of the two structured outlier strata (matching
      simulate.py's education=2/wealth=2/rural=1 -> +0.90 and
      education=0/wealth=2/rural=0 -> -0.60 pattern, replicated across both
      sex strata). `severity=0` disables spikes (pure `tail_family` noise);
      `severity=1` reproduces the original mechanism exactly (regression-
      checked in `python/tests/test_sparse_strata_robustness.py`);
      `severity=2`/`3` scale the spike magnitude 2x/3x, holding location and
      the +0.90:-0.60 asymmetry ratio fixed.

    No detection bias is modeled here (`y == y_true` always) -- this
    generator isolates the "non-Gaussian random effects" assumption from
    detection-bias violations, which are covered separately (see
    `simulate_severity_dependent_detection` and the joint capstone scenario).
    """
    if tail_family not in ("gaussian", "student_t3", "skew_lognormal"):
        raise ValueError(f"Unknown tail_family: {tail_family!r}")
    if severity not in (0, 1, 2, 3):
        raise ValueError("severity must be one of 0, 1, 2, 3")

    rng = np.random.default_rng(seed)
    strata = _stratum_table()
    k = len(strata)

    if sparse:
        weights = rng.gamma(shape=0.35, scale=1.0, size=k)
        weights = weights / weights.sum()
    else:
        weights = np.repeat(1.0 / k, k)

    stratum_id = rng.choice(k, size=n, p=weights)
    x = strata.iloc[stratum_id].reset_index(drop=True)

    linear_predictor = (
        prevalence_shift
        + 0.20 * x["sex"].to_numpy()
        + 0.35 * x["education"].to_numpy()
        + 0.30 * x["wealth"].to_numpy()
        + 0.25 * x["rural"].to_numpy()
    )

    rng_re = np.random.default_rng(seed + 999)
    if tail_family == "gaussian":
        raw = rng_re.normal(loc=0.0, scale=interaction_sd, size=k)
    elif tail_family == "student_t3":
        # Student-t(df=3) has variance df/(df-2) = 3; rescale to match
        # interaction_sd so only kurtosis (heavy tails), not dispersion, varies.
        raw = rng_re.standard_t(df=3, size=k) * (interaction_sd / np.sqrt(3.0))
    else:  # skew_lognormal
        raw = rng_re.lognormal(mean=0.0, sigma=0.5, size=k)
        raw = (raw - raw.mean()) / raw.std() * interaction_sd

    if severity > 0:
        high_magnitude = _SPIKE_HIGH_UNIT * severity
        low_magnitude = _SPIKE_LOW_UNIT * severity
        raw = raw + high_magnitude * (
            (strata["education"].to_numpy() == 2)
            & (strata["wealth"].to_numpy() == 2)
            & (strata["rural"].to_numpy() == 1)
        )
        raw = raw + low_magnitude * (
            (strata["education"].to_numpy() == 0)
            & (strata["wealth"].to_numpy() == 2)
            & (strata["rural"].to_numpy() == 0)
        )

    true_stratum_residual = raw - raw.mean()
    linear_predictor = linear_predictor + true_stratum_residual[stratum_id]

    p_true = _logit_inv(linear_predictor)
    y_true = rng.binomial(1, p_true)

    out = x.copy()
    out["stratum"] = pd.Categorical(stratum_id.astype(str))
    for col in ["sex", "education", "wealth", "rural"]:
        out[col] = pd.Categorical(out[col])
    out["y"] = y_true.astype(int)
    out["y_true"] = y_true.astype(int)
    out["detection_probability"] = np.ones(n)
    out["true_stratum_residual"] = true_stratum_residual[stratum_id]
    return out


def misspecified_score(df: pd.DataFrame, kind: str) -> np.ndarray:
    """A stratum-level detection `score` deliberately misspecified relative to
    the true generating mechanism (`education + wealth + 0.4*rural`), for
    testing `correct_detection_bias()`'s robustness to the analyst not
    knowing the true detection-driving covariates.

    Computed from the SAME aggregation `correct_detection_bias()` uses
    internally (`_aggregate_strata`), not a static 36-row stratum table --
    stratum IDs are stored as strings and sort lexicographically ("10" <
    "2"), and only strata actually present in `df` appear, so a
    precomputed/static score vector would silently misalign. Pass the same
    `df` given to `correct_detection_bias`/`vpc_detection_bounds`.

    Parameters
    ----------
    df : DataFrame
        Individual-level data, as passed to `correct_detection_bias()`.
    kind : {"omit_rural", "wrong_covariate", "quadratic"}
        - "omit_rural": drops the rural term (`education + wealth`).
        - "wrong_covariate": swaps education for sex (`sex + wealth`), a
          covariate that plays no role in the true detection mechanism.
        - "quadratic": wrong functional form, same covariates
          (`education**2 + wealth**2 + 0.4*rural`).
    """
    strata = _aggregate_strata(df)
    education = _numeric(strata["education"])
    wealth = _numeric(strata["wealth"])
    rural = _numeric(strata["rural"])
    sex = _numeric(strata["sex"])

    if kind == "omit_rural":
        return education + wealth
    elif kind == "wrong_covariate":
        return sex + wealth
    elif kind == "quadratic":
        return education**2 + wealth**2 + 0.4 * rural
    else:
        raise ValueError(f"Unknown kind: {kind!r}")


def simulate_k_strata(
    dims: Sequence[int] = (2, 3, 3, 2),
    betas: Sequence[float] | None = None,
    n: int = 3500,
    prevalence_shift: float = -3.00,
    interaction_sd: float = 0.15,
    sparse: bool = True,
    seed: int = 42,
) -> pd.DataFrame:
    """Generalized intersectional design with an arbitrary number of strata.

    Builds the full factorial stratum table over generic covariates
    ``x1..xd`` with the level counts in ``dims`` (K = prod(dims); e.g.
    ``(2,2,2)`` -> K=8, ``(2,3,3,2)`` -> K=36 matching the package's fixed
    design, ``(2,3,3,3,2)`` -> K=108). The linear predictor is
    ``prevalence_shift + sum_i beta_i * x_i`` plus i.i.d. Gaussian stratum
    residuals ``N(0, interaction_sd)`` (no structured spikes — kept clean so
    the K-study isolates the number-of-strata effect). ``betas`` defaults to
    cycling the canonical gradient magnitudes ``[0.20, 0.35, 0.30, 0.25]``.

    Intended for studying how ``sparse_strata_vpc()``'s finite-K calibration
    behaves as K varies, so the output carries only the columns that function
    needs (``stratum, y, y_true``) plus ``true_stratum_logit`` — the exact
    per-stratum eta, from which the analytic population truth is
    ``vpc_latent(var(unique etas, ddof=1))``. It does NOT carry the named
    ``sex/education/wealth/rural`` columns, so it is not consumable by
    ``fit_imaihda`` or ``correct_detection_bias`` (deliberately out of the
    K-study's scope).
    """
    dims = tuple(int(d) for d in dims)
    if len(dims) < 1 or any(d < 2 for d in dims):
        raise ValueError("dims must have at least one dimension, each with >= 2 levels")
    canonical = [0.20, 0.35, 0.30, 0.25]
    if betas is None:
        betas = [canonical[i % len(canonical)] for i in range(len(dims))]
    betas = np.asarray(betas, dtype=float)
    if len(betas) != len(dims):
        raise ValueError(f"betas has length {len(betas)} but dims has {len(dims)} dimensions")

    grids = np.meshgrid(*[np.arange(d) for d in dims], indexing="ij")
    table = np.column_stack([g.ravel() for g in grids]).astype(float)
    k = table.shape[0]

    rng = np.random.default_rng(seed)
    if sparse:
        weights = rng.gamma(shape=0.35, scale=1.0, size=k)
        weights = weights / weights.sum()
    else:
        weights = np.repeat(1.0 / k, k)

    eta_add = prevalence_shift + table @ betas
    rng_re = np.random.default_rng(seed + 999)
    raw = rng_re.normal(0.0, interaction_sd, size=k)
    residual = raw - raw.mean()
    eta = eta_add + residual

    stratum_id = rng.choice(k, size=n, p=weights)
    p_true = _logit_inv(eta[stratum_id])
    y = rng.binomial(1, p_true)

    out = pd.DataFrame({
        "stratum": pd.Categorical(stratum_id.astype(str)),
        "y": y.astype(int),
        "y_true": y.astype(int),
        "true_stratum_logit": eta[stratum_id],
    })
    return out


def simulate_selective_attrition(
    n: int = 6000,
    prevalence_shift: float = -2.10,
    interaction_sd: float = 0.0,
    attrition_strength: float = 0.0,
    attrition_outcome_weight: float = 0.30,
    sparse: bool = False,
    seed: int = 42,
) -> pd.DataFrame:
    """Like ``simulate_intersectional_data`` but individuals are selectively
    LOST before observation (the fourth artefact channel named in the PhD
    exposé, alongside prevalence, sparsity, and detection). Retention is
    SES-patterned and mildly outcome-dependent:

        logit(P(retained)) = 2.0 - attrition_strength * (education + wealth
                             + 0.4*rural) - attrition_outcome_weight
                             * attrition_strength * y_true

    so disadvantaged strata — and, within them, true cases — are
    differentially dropped. Unlike detection (which zeroes the recorded
    outcome but keeps the individual in the denominator), attrition removes
    the individual entirely, shrinking BOTH the stratum's event count and its
    size. ``attrition_strength=0`` reproduces ``simulate_intersectional_data``
    exactly (no rows dropped; regression-checked).

    Returns the standard column contract (``stratum, sex, education, wealth,
    rural, y, y_true, detection_probability, true_stratum_residual``) plus
    ``retention_probability``, restricted to retained rows — so
    ``fit_imaihda``, ``sparse_strata_vpc``, etc. consume it unmodified.
    """
    rng = np.random.default_rng(seed)
    strata = _stratum_table()
    k = len(strata)

    if sparse:
        weights = rng.gamma(shape=0.35, scale=1.0, size=k)
        weights = weights / weights.sum()
    else:
        weights = np.repeat(1.0 / k, k)

    stratum_id = rng.choice(k, size=n, p=weights)
    x = strata.iloc[stratum_id].reset_index(drop=True)

    linear_predictor = (
        prevalence_shift
        + 0.20 * x["sex"].to_numpy()
        + 0.35 * x["education"].to_numpy()
        + 0.30 * x["wealth"].to_numpy()
        + 0.25 * x["rural"].to_numpy()
    )

    true_stratum_residual = np.zeros(k)
    if interaction_sd > 0:
        rng_re = np.random.default_rng(seed + 999)
        raw = rng_re.normal(loc=0.0, scale=interaction_sd, size=k)
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
        true_stratum_residual = raw - raw.mean()
        linear_predictor = linear_predictor + true_stratum_residual[stratum_id]

    p_true = _logit_inv(linear_predictor)
    y_true = rng.binomial(1, p_true)

    if attrition_strength > 0:
        score = (
            x["education"].to_numpy()
            + x["wealth"].to_numpy()
            + 0.40 * x["rural"].to_numpy()
        )
        retain_logit = (
            2.0
            - attrition_strength * score
            - attrition_outcome_weight * attrition_strength * y_true
        )
        retain_p = _logit_inv(retain_logit)
        retained = rng.binomial(1, retain_p).astype(bool)
    else:
        retain_p = np.ones(n)
        retained = np.ones(n, dtype=bool)

    out = x.copy()
    out["stratum"] = pd.Categorical(stratum_id.astype(str))
    for col in ["sex", "education", "wealth", "rural"]:
        out[col] = pd.Categorical(out[col])
    out["y"] = y_true.astype(int)
    out["y_true"] = y_true.astype(int)
    out["detection_probability"] = np.ones(n)
    out["true_stratum_residual"] = true_stratum_residual[stratum_id]
    out["retention_probability"] = retain_p
    return out.loc[retained].reset_index(drop=True)
