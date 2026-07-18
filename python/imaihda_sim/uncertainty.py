"""Small-sample bias correction and confidence intervals for sparse strata.

The fast null-model VPC estimator (`fit_imaihda(method="fast")`) computes
between-stratum variance from Laplace-smoothed empirical logits
`(events + 0.5) / (n + 1)`. When a stratum has few individuals, that smoothing
pulls its logit toward the population mean, which *shrinks* the observed
between-stratum spread below its true value. The estimator's noise-subtraction
step does not correct for this shrinkage, so in sparse-strata regimes the fast
VPC is systematically **underestimated** relative to the true value — and gets
worse exactly where analysts most need a warning. No confidence interval or
small-sample correction for this estimator has been reported in the I-MAIHDA
literature; CRAN `MAIHDA` sidesteps the issue by using full GLMM asymptotics
instead of a fast estimator.

This module corrects it by **calibration**: for the observed stratum sizes, it
simulates the null-model estimator's own expected value under a grid of
candidate true between-stratum variances (holding stratum sizes fixed), and
inverts that curve to find the true variance whose expected naive estimate
matches what was actually observed — an indirect-inference / SIMEX-style
correction. Bootstrap replicates of the naive estimator at the calibrated
variance are pushed through the same inverse mapping to build a confidence
interval for the bias-corrected VPC.

Note on the variance formula: this module's "naive" estimator uses
**unweighted** (sample) variance of the stratum residuals, matching
``imaihda/R/diagnostics.R`` — the R package switched to this formula in
v0.2.1 specifically because precision-weighted variance downweights the
extreme strata that carry the most between-stratum signal, which introduces
its own persistent bias that does not vanish even at very large stratum
sizes (verified empirically; not a small-sample effect). Python's
``fit_imaihda()`` keeps the older precision-weighted formula as its
*default* (for benchmark continuity) but now also exposes
``fit_imaihda(df, weighting="unweighted")``, which matches this module's
``vpc_null_naive`` exactly (asserted in ``tests/test_fit_weighting.py``).
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from .fit import vpc_latent

__all__ = ["sparse_strata_vpc"]

# A stratum smaller than this is where Laplace-smoothing shrinkage is largest
# and the naive fast estimator's downward bias becomes material.
SPARSE_THRESHOLD = 20


def _logit(p):
    return np.log(p / (1.0 - p))


def _expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def _var_null_from_counts(events: np.ndarray, n: np.ndarray, p_overall: float) -> np.ndarray:
    """Null-model between-stratum variance, unweighted formula (matches
    ``imaihda/R/diagnostics.R``'s ``between_stratum_variance()``).

    Vectorized over a leading replicate axis: ``events``/``n`` have shape
    ``(n_rep, n_strata)`` and the return has shape ``(n_rep,)``.
    """
    p_s = (events + 0.5) / (n + 1.0)
    logit_s = np.log(p_s / (1.0 - p_s))
    sampling_var = 1.0 / (events + 0.5) + 1.0 / (n - events + 0.5)
    null_pred = math.log(p_overall / (1.0 - p_overall))
    resid = logit_s - null_pred
    raw_var = np.var(resid, axis=-1, ddof=1)  # unweighted, (K - 1) divisor
    expected_noise = np.mean(sampling_var, axis=-1)
    return np.maximum(raw_var - expected_noise, 0.0)


def _simulate_naive_var(
    sigma2: float, n_j: np.ndarray, p_overall: float, n_rep: int, rng: np.random.Generator
) -> np.ndarray:
    """Sampling distribution of the naive var_null estimator under a candidate
    true between-stratum variance `sigma2`, holding stratum sizes `n_j` fixed."""
    sd = math.sqrt(max(sigma2, 0.0))
    K = len(n_j)
    effects = rng.normal(0.0, sd, size=(n_rep, K)) if sd > 0 else np.zeros((n_rep, K))
    p_j = np.clip(_expit(_logit(p_overall) + effects), 1e-6, 1 - 1e-6)
    events = rng.binomial(n_j[None, :], p_j).astype(float)
    n_rep_arr = np.tile(n_j.astype(float), (n_rep, 1))
    return _var_null_from_counts(events, n_rep_arr, p_overall)


def _build_calibration_curve(
    n_j: np.ndarray, p_overall: float, sigma2_max: float, n_grid: int, n_sim: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Monte Carlo calibration curve E[naive var_null | true sigma2] over a grid."""
    grid = np.linspace(0.0, sigma2_max, n_grid)
    means = np.array([_simulate_naive_var(s, n_j, p_overall, n_sim, rng).mean() for s in grid])
    return grid, means


def sparse_strata_vpc(
    df: pd.DataFrame,
    sigma2_max: float = 2.0,
    n_grid: int = 25,
    n_sim: int = 200,
    n_boot: int = 300,
    level: float = 0.95,
    seed: Optional[int] = None,
    ci_method: str = "bootstrap",
) -> Dict[str, float]:
    """Bias-corrected null-model VPC and confidence interval for sparse strata.

    Parameters
    ----------
    df : DataFrame
        Individual-level data with columns ``stratum, y`` at minimum (as
        produced by :func:`simulate_intersectional_data`).
    sigma2_max : float, default 2.0
        Upper bound of the calibration grid on the between-stratum variance
        (logit scale). Increase if ``vpc_null_corrected`` saturates near the
        implied VPC at ``sigma2_max``.
    n_grid : int, default 25
        Number of grid points spanning ``[0, sigma2_max]`` for the
        calibration curve.
    n_sim : int, default 200
        Monte Carlo replicates per grid point used to estimate the expected
        naive estimator (and, for ``ci_method="test_inversion"``, its
        per-grid-point quantiles).
    n_boot : int, default 300
        Bootstrap replicates (at the calibrated variance) used to build the
        confidence interval when ``ci_method="bootstrap"``. Unused by
        ``"test_inversion"``.
    level : float, default 0.95
        Confidence level.
    seed : int, optional
        RNG seed for reproducibility.
    ci_method : {"bootstrap", "test_inversion"}, default "bootstrap"
        ``"bootstrap"`` (the historical default, kept for continuity of
        published results) resamples the naive statistic at the calibrated
        variance and maps replicates through the inverse curve. It
        degenerates to a zero-width interval when the naive statistic sits
        at the ``max(0, .)`` truncation floor — which happens routinely in
        extreme-K / extreme-sparsity regimes (many strata, few individuals
        each). ``"test_inversion"`` instead reports
        ``{sigma2 : q_{a/2}(T|sigma2) <= T_obs <= q_{1-a/2}(T|sigma2)}`` —
        the set of true variances under which the observed statistic is not
        extreme — which stays honestly wide at the floor (a one-sided
        ``[0, upper]`` interval rather than ``[0, 0]``).

    Returns
    -------
    dict
        ``vpc_null_naive, vpc_null_corrected, var_null_naive,
        var_null_corrected, ci_lower, ci_upper, level, at_floor, ci_method,
        n_strata, min_stratum_n, sparse, n_sim, n_boot``.
    """
    if not (0.0 < level < 1.0):
        raise ValueError("level must be in (0, 1)")
    if ci_method not in ("bootstrap", "test_inversion"):
        raise ValueError(f"ci_method must be 'bootstrap' or 'test_inversion', got {ci_method!r}")

    strata = df.groupby("stratum", observed=True).agg(n=("y", "size"), events=("y", "sum")).reset_index()
    n_j = strata["n"].to_numpy(dtype=np.int64)
    events_j = strata["events"].to_numpy(dtype=float)
    if len(n_j) < 2:
        raise ValueError("Need at least 2 strata; found " + str(len(n_j)))
    p_overall = float(np.clip(events_j.sum() / n_j.sum(), 1e-6, 1 - 1e-6))

    var_naive = float(_var_null_from_counts(events_j[None, :], n_j.astype(float)[None, :], p_overall)[0])

    rng = np.random.default_rng(seed)
    grid = np.linspace(0.0, sigma2_max, n_grid)
    sims = [_simulate_naive_var(s, n_j, p_overall, n_sim, rng) for s in grid]
    means = np.array([s.mean() for s in sims])

    at_floor = bool(var_naive <= means[0])
    if at_floor:
        # The naive estimate is already at or below the estimator's own floor
        # at sigma2 = 0; no downward correction is warranted.
        var_corrected = 0.0
    elif var_naive >= means[-1]:
        # Off the top of the calibration grid; extrapolating is unreliable.
        var_corrected = float(grid[-1])
    else:
        var_corrected = float(np.interp(var_naive, means, grid))

    alpha = 1.0 - level
    if ci_method == "bootstrap":
        boot_naive = _simulate_naive_var(var_corrected, n_j, p_overall, n_boot, rng)
        boot_corrected = np.interp(np.clip(boot_naive, means[0], means[-1]), means, grid)
        lo, hi = np.percentile(boot_corrected, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    else:
        # Monte Carlo test inversion: keep every grid sigma2 under which the
        # observed statistic is within the central (1 - alpha) band of the
        # statistic's sampling distribution.
        q_lo = np.array([np.quantile(s, alpha / 2) for s in sims])
        q_hi = np.array([np.quantile(s, 1 - alpha / 2) for s in sims])
        accept = (q_lo <= var_naive) & (var_naive <= q_hi)
        if accept.any():
            lo = float(grid[accept][0])
            hi = float(grid[accept][-1])
        else:
            # Degenerate MC corner case: fall back to the point estimate.
            lo = hi = var_corrected

    return {
        "vpc_null_naive": vpc_latent(var_naive),
        "vpc_null_corrected": vpc_latent(var_corrected),
        "var_null_naive": var_naive,
        "var_null_corrected": var_corrected,
        "ci_lower": vpc_latent(float(lo)),
        "ci_upper": vpc_latent(float(hi)),
        "level": level,
        "at_floor": at_floor,
        "ci_method": ci_method,
        "n_strata": int(len(n_j)),
        "min_stratum_n": int(n_j.min()),
        "sparse": bool(n_j.min() < SPARSE_THRESHOLD),
        "n_sim": n_sim,
        "n_boot": n_boot,
    }
