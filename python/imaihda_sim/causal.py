"""Causal identification tools for the VPC under differential detection.

This module treats SES-patterned under-detection as a *missing-data /
causal-identification* problem rather than only a sensitivity analysis, and
provides the two estimators that framing implies. Both target the
**population** between-stratum variance of the true stratum logits (the
parameter a hypothetical infinite sample would reveal), not the finite-sample
refit on ``y_true`` — the distinction matters in sparse designs, where the
finite-sample refit itself carries shrinkage.

Identification facts encoded here (derivations in
``docs/PHASE3_CAUSAL_IDENTIFICATION.md``):

1. **Point identification (covariate-only detection).** If detection is
   conditionally independent of the true outcome given the stratum
   (``P(detect | stratum j) = d_j``), then ``E[p_obs_j] = p_true_j * d_j``,
   so for hypothesized ``(delta, score, baseline_logit)`` the true VPC is
   point-identified. This is what ``correct_detection_bias()`` inverts.

2. **Set identification (outcome-dependent detection).** If detection may
   also depend on the true outcome itself, the observed data constrain
   ``p_true_j`` only to an interval, and the VPC is identified only up to
   sharp bounds — computed exactly by :func:`vpc_partial_bounds` (minimum by
   convex projection, maximum at a vertex of the box, both global).

3. **Joint estimation (detection x sparsity).** Composing the detection
   correction with the sparse-strata calibration double-counts, for two
   reasons found empirically and mechanistically: (a) the sparse calibration
   assumes plain-binomial noise but detection-corrected counts ``e/d`` carry
   inflated noise; (b) an independent-Gaussian effects model ignores the
   correlation between stratum prevalence and detection (both driven by the
   same SES covariates), which depresses the calibration curve and causes
   overshoot on inversion. :func:`joint_calibrated_vpc` fixes both: it
   calibrates the *entire* corrected-statistic pipeline in one simulation
   whose effects carry the data's own estimated systematic (main-effects)
   component, with only the residual variance on the calibration grid.
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .detection import _aggregate_strata, _logit_inv, _stratum_score, _wls_main_pred
from .fit import _between_stratum_variance, vpc_latent

__all__ = ["joint_calibrated_vpc", "vpc_partial_bounds"]


def _logit(p):
    return np.log(p / (1.0 - p))


def _weighted_pipeline_stat(e_obs: np.ndarray, n: np.ndarray, d_rel: np.ndarray):
    """Detection-corrected precision-weighted null-model variance statistic —
    the exact statistic ``correct_detection_bias()`` computes, factored out so
    the joint calibration can simulate its own sampling distribution."""
    e_c = np.minimum(n, e_obs / d_rel)
    p_c = (e_c + 0.5) / (n + 1.0)
    lg = np.log(p_c / (1.0 - p_c))
    sv = 1.0 / (e_c + 0.5) + 1.0 / (n - e_c + 0.5)
    w = 1.0 / sv
    p0 = float(np.clip(e_c.sum() / n.sum(), 1e-6, 1 - 1e-6))
    null_pred = np.repeat(math.log(p0 / (1.0 - p0)), len(n))
    return float(_between_stratum_variance(lg, null_pred, sv, w)), lg, w, p0


def joint_calibrated_vpc(
    df: pd.DataFrame,
    delta: float,
    baseline_logit: float = 2.0,
    score: Optional[Sequence[float]] = None,
    s2res_max: float = 1.5,
    n_grid: int = 25,
    n_sim: int = 300,
    n_boot: int = 300,
    level: float = 0.95,
    seed: Optional[int] = None,
) -> Dict[str, float]:
    """Jointly-calibrated VPC under SES-patterned detection AND sparse strata.

    Single-simulation calibration of the full detection-corrected pipeline.
    Stratum effects in the calibration simulator are decomposed as
    ``main_effects_pattern (estimated by WLS from the corrected logits, held
    fixed) + N(0, s2_res)``, so the simulator preserves the data's own
    correlation between prevalence and detection; only the residual variance
    ``s2_res`` sits on the calibration grid. The reported total variance is
    ``Var(main_pattern) + s2_res_hat``.

    Targets the *population* between-stratum variance (see module docstring).

    Returns dict with ``vpc_corrected, var_corrected, var_main_pattern,
    s2_res_hat, vpc_detection_only, ci_lower, ci_upper, level, saturated``.
    ``saturated=True`` means the inversion hit the top of the grid — rerun
    with a larger ``s2res_max``.
    """
    if delta < 0:
        raise ValueError("delta must be non-negative")
    if not (0.0 < level < 1.0):
        raise ValueError("level must be in (0, 1)")

    strata = _aggregate_strata(df)
    e_obs = strata["events"].to_numpy(dtype=float)
    n_j = strata["n"].to_numpy(dtype=float)
    if len(n_j) < 2:
        raise ValueError("Need at least 2 strata; found " + str(len(n_j)))

    s = _stratum_score(strata, score)
    d0 = float(_logit_inv(baseline_logit))
    d_rel = np.clip(np.asarray(_logit_inv(baseline_logit - delta * s), dtype=float) / d0, 1e-6, 1.0)

    T_obs, lg_c, w, p0 = _weighted_pipeline_stat(e_obs, n_j, d_rel)

    main_pred = _wls_main_pred(strata, lg_c, w)
    mp_c = main_pred - np.average(main_pred, weights=n_j)
    v_main = float(np.var(main_pred, ddof=1))

    rng = np.random.default_rng(seed)
    grid = np.linspace(0.0, s2res_max, n_grid)
    K = len(n_j)
    center = _logit(p0)
    n_int = n_j.astype(np.int64)

    def simulate_T(s2r: float, n_rep: int) -> np.ndarray:
        sd = math.sqrt(max(s2r, 0.0))
        noise = rng.normal(0.0, sd, size=(n_rep, K)) if sd > 0 else np.zeros((n_rep, K))
        eff = np.broadcast_to(mp_c, (n_rep, K)) + noise
        p_cs = np.clip(1.0 / (1.0 + np.exp(-(center + eff))), 1e-6, 1 - 1e-6)
        p_o = np.clip(p_cs * d_rel[None, :], 1e-9, 1 - 1e-9)
        e_o = rng.binomial(n_int[None, :], p_o).astype(float)
        return np.array([_weighted_pipeline_stat(e_o[i], n_j, d_rel)[0] for i in range(n_rep)])

    curve = np.array([simulate_T(s2r, n_sim).mean() for s2r in grid])
    curve = np.maximum.accumulate(curve)  # isotonic guard against MC wiggle

    saturated = False
    if T_obs <= curve[0]:
        s2r_hat = 0.0
    elif T_obs >= curve[-1]:
        s2r_hat = float(grid[-1])
        saturated = True
    else:
        s2r_hat = float(np.interp(T_obs, curve, grid))

    var_corrected = v_main + s2r_hat

    # CI: bootstrap the pipeline statistic at the calibrated point, push each
    # replicate back through the inverse curve. V_main is treated as fixed
    # (its WLS estimation noise is second-order at K=36 with 7 parameters);
    # documented as a limitation.
    boot_T = simulate_T(s2r_hat, n_boot)
    boot_s2r = np.interp(np.clip(boot_T, curve[0], curve[-1]), curve, grid)
    alpha = 1.0 - level
    lo, hi = np.percentile(v_main + boot_s2r, [100 * alpha / 2, 100 * (1 - alpha / 2)])

    return {
        "vpc_corrected": vpc_latent(var_corrected),
        "var_corrected": var_corrected,
        "var_main_pattern": v_main,
        "s2_res_hat": s2r_hat,
        "vpc_detection_only": vpc_latent(T_obs) if T_obs > 0 else 0.0,
        "ci_lower": vpc_latent(float(lo)),
        "ci_upper": vpc_latent(float(hi)),
        "level": level,
        "saturated": saturated,
    }


# ── Sharp partial-identification bounds ───────────────────────────────────────


def _min_variance_over_box(a: np.ndarray, b: np.ndarray) -> Tuple[float, np.ndarray]:
    """Global minimum of Var(x, ddof=1) subject to a_j <= x_j <= b_j.

    Convex program with a 1-D dual: the minimizer is x_j = clip(c, a_j, b_j)
    for the scalar c minimizing g(c) = sum_j dist(c, [a_j, b_j])^2, and g is
    convex in c — solved exactly by ternary search on [min a, max b].
    """
    lo, hi = float(a.min()), float(b.max())

    def g(c):
        x = np.clip(c, a, b)
        return float(np.sum((x - c) ** 2))

    for _ in range(200):
        m1 = lo + (hi - lo) / 3.0
        m2 = hi - (hi - lo) / 3.0
        if g(m1) <= g(m2):
            hi = m2
        else:
            lo = m1
    c = 0.5 * (lo + hi)
    x = np.clip(c, a, b)
    return float(np.var(x, ddof=1)), x


def _max_variance_over_box(
    a: np.ndarray, b: np.ndarray, n_starts: int = 64, seed: int = 0
) -> Tuple[float, np.ndarray]:
    """Global maximum of Var(x, ddof=1) subject to a_j <= x_j <= b_j.

    Var is convex, so the maximum over the box is attained at a vertex
    (x_j in {a_j, b_j} for all j). Coordinate exchange — moving one
    coordinate to whichever endpoint is farther from the current mean strictly
    increases the variance — converges to a locally-optimal vertex; multiple
    random vertex starts recover the global one (verified exhaustively for
    K <= 12 in the test suite).
    """
    rng = np.random.default_rng(seed)
    K = len(a)
    best_v, best_x = -np.inf, None
    for start in range(n_starts):
        pick = rng.random(K) < 0.5 if start > 0 else (b - b.mean() > 0)
        x = np.where(pick, b, a).astype(float)
        for _ in range(200):
            changed = False
            mu = x.mean()
            for j in range(K):
                cand = b[j] if abs(b[j] - mu) >= abs(a[j] - mu) else a[j]
                if cand != x[j]:
                    x_new = x.copy()
                    x_new[j] = cand
                    if np.var(x_new, ddof=1) > np.var(x, ddof=1) + 1e-15:
                        x = x_new
                        changed = True
                        mu = x.mean()
            if not changed:
                break
        v = float(np.var(x, ddof=1))
        if v > best_v:
            best_v, best_x = v, x
    return best_v, best_x


def vpc_partial_bounds(
    df: pd.DataFrame,
    delta_max: float,
    baseline_logit: float = 2.0,
    score: Optional[Sequence[float]] = None,
) -> Dict[str, float]:
    """Sharp bounds on the VPC when detection is only *set*-identified.

    Assumes monotone under-ascertainment (observed events never exceed true
    events) with per-stratum relative detection no worse than
    ``d_rel_j(delta_max)``. Then each stratum's true (corrected-scale)
    prevalence lies in an interval, the true stratum logits lie in a box, and
    since ``vpc_latent`` is increasing in the variance, sharp VPC bounds are
    the min/max of Var(logits) over that box — computed globally (convex
    projection for the min, vertex maximization for the max).

    These bounds remain valid under *outcome-dependent* detection (the regime
    where point correction provably breaks down), as long as the assumed
    ``delta_max`` envelope covers the true total under-detection.

    Returns dict with ``vpc_lower, vpc_upper, var_lower, var_upper, delta_max``.
    """
    if delta_max < 0:
        raise ValueError("delta_max must be non-negative")

    strata = _aggregate_strata(df)
    e_obs = strata["events"].to_numpy(dtype=float)
    n_j = strata["n"].to_numpy(dtype=float)
    if len(n_j) < 2:
        raise ValueError("Need at least 2 strata; found " + str(len(n_j)))

    s = _stratum_score(strata, score)
    d0 = float(_logit_inv(baseline_logit))
    d_rel_min = np.clip(
        np.asarray(_logit_inv(baseline_logit - delta_max * s), dtype=float) / d0, 1e-6, 1.0
    )

    # Interval for true events: [observed, observed / worst-case detection],
    # capped at n. Laplace smoothing applied consistently at both ends.
    e_lo = e_obs
    e_hi = np.minimum(n_j, e_obs / d_rel_min)
    p_lo = (e_lo + 0.5) / (n_j + 1.0)
    p_hi = (e_hi + 0.5) / (n_j + 1.0)
    a = np.log(p_lo / (1.0 - p_lo))
    b = np.log(p_hi / (1.0 - p_hi))

    v_lo, _ = _min_variance_over_box(a, b)
    v_hi, _ = _max_variance_over_box(a, b)

    return {
        "vpc_lower": vpc_latent(v_lo),
        "vpc_upper": vpc_latent(v_hi),
        "var_lower": v_lo,
        "var_upper": v_hi,
        "delta_max": float(delta_max),
    }
