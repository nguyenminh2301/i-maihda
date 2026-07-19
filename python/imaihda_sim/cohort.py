"""Cross-cohort VPC/PCV gap decomposition — the reproducibility companion's
answer to the PhD's central methodological question.

The exposé (`nguyenminh2301/i-maihda`, WZB Doctoral Position 400) asks whether
the additive-dominant intersectional pattern of late-life multimorbidity
generalizes across high- and middle-income aging cohorts (HRS, ELSA, SHARE,
CHARLS), and warns that a raw HIC–MIC difference in the variance partition
coefficient (VPC) or proportional change in variance (PCV) must not be
mistaken for a genuine structural difference when it could instead arise from
four artefact channels: **outcome prevalence, sparse strata, selective
attrition, and SES-patterned under-detection (differential diagnosis).**

This module operationalizes exactly that. Given two cohort specifications
(same intersectional design; per-cohort values for each channel plus the true
between-stratum structure), it decomposes the observed VPC/PCV gap into a
share for each of the four artefact channels plus a **structural residual** —
the part attributable to a genuine difference in intersectional structure. The
attribution is a **Shapley decomposition** (Shorrocks 2013): the unique
attribution that is symmetric across channels, exactly additive (shares sum to
the observed gap), and order-independent — unlike a single sequential
(Oaxaca-style) path, whose components depend on the arbitrary channel ordering.

This is a *design/sensitivity* tool, not an estimator for real cohort VPCs
(the empirical PhD work uses full GLMM I-MAIHDA). Its role, as the exposé
states, is to pre-specify when VPC/PCV can and cannot support cross-cohort
comparison: it answers "how much of a gap of this size could these nuisance
differences alone produce, and is the residual structural difference
distinguishable from zero?" — self-validated below because the true structural
gap is known by construction.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from math import factorial
from typing import Dict, List, Optional, Sequence

import numpy as np

from .fit import vpc_latent

__all__ = ["CohortSpec", "cross_cohort_decomposition"]

_BETAS = np.array([0.20, 0.35, 0.30, 0.25])

# The five channels flipped between cohort A and cohort B. "structure" is the
# scientific quantity of interest; the other four are the exposé's artefacts.
CHANNELS = ["prevalence", "structure", "sparsity", "detection", "attrition"]
_ARTEFACT_CHANNELS = ["prevalence", "sparsity", "detection", "attrition"]


def _stratum_table() -> np.ndarray:
    rows = []
    for sex in (0, 1):
        for edu in (0, 1, 2):
            for wealth in (0, 1, 2):
                for rural in (0, 1):
                    rows.append((sex, edu, wealth, rural))
    return np.array(rows, dtype=float)


_TAB = _stratum_table()
_K = len(_TAB)
_SCORE = _TAB[:, 1] + _TAB[:, 2] + 0.40 * _TAB[:, 3]


@dataclass
class CohortSpec:
    """One cohort's data-generating configuration on the shared 36-stratum
    (2×3×3×2) design. The four artefact channels plus the structural
    parameters are what the decomposition flips between cohorts.

    prevalence_shift : intercept on the logit scale (event prevalence).
    interaction_sd, structure_seed : the true between-stratum structure
        (residual interaction magnitude and its RNG draw). Two cohorts share
        structure iff they share BOTH of these.
    sparse : gamma-skewed stratum allocation (True) vs. equal (False).
    detection_strength : SES-patterned under-detection (0 = none).
    attrition_strength : SES/outcome-patterned selective dropout (0 = none).
    n : individuals drawn before attrition.
    """

    prevalence_shift: float = -2.10
    interaction_sd: float = 0.30
    structure_seed: int = 999
    sparse: bool = False
    detection_strength: float = 0.0
    attrition_strength: float = 0.0
    n: int = 4000


def _true_eta(prevalence_shift: float, interaction_sd: float, structure_seed: int) -> np.ndarray:
    add = prevalence_shift + _TAB @ _BETAS
    rng = np.random.default_rng(structure_seed)
    raw = rng.normal(0.0, interaction_sd, _K)
    raw += 0.90 * ((_TAB[:, 1] == 2) & (_TAB[:, 2] == 2) & (_TAB[:, 3] == 1))
    raw -= 0.60 * ((_TAB[:, 1] == 0) & (_TAB[:, 2] == 2) & (_TAB[:, 3] == 0))
    resid = raw - raw.mean()
    return add + resid


def _simulate(spec: CohortSpec, seed: int):
    rng = np.random.default_rng(seed)
    eta = _true_eta(spec.prevalence_shift, spec.interaction_sd, spec.structure_seed)
    if spec.sparse:
        w = rng.gamma(0.35, 1.0, _K)
        w = w / w.sum()
    else:
        w = np.repeat(1.0 / _K, _K)
    sid = rng.choice(_K, size=spec.n, p=w)
    p_true = 1.0 / (1.0 + np.exp(-eta[sid]))
    y_true = rng.binomial(1, p_true)
    y = y_true
    if spec.detection_strength > 0:
        dp = 1.0 / (1.0 + np.exp(-(2.0 - spec.detection_strength * (_TAB[sid, 1] + _TAB[sid, 2])
                                   - 0.40 * spec.detection_strength * _TAB[sid, 3])))
        y = y_true * rng.binomial(1, dp)
    if spec.attrition_strength > 0:
        keep_logit = 2.0 - spec.attrition_strength * _SCORE[sid] - 0.30 * spec.attrition_strength * y_true
        keep = rng.binomial(1, 1.0 / (1.0 + np.exp(-keep_logit))).astype(bool)
        sid, y, y_true = sid[keep], y[keep], y_true[keep]
    return sid, y, y_true


def _stratum_stats(sid: np.ndarray, y: np.ndarray):
    e = np.zeros(_K)
    n = np.zeros(_K)
    np.add.at(e, sid, y)
    np.add.at(n, sid, 1)
    mask = n > 0
    return e[mask], n[mask]


def _metric(sid: np.ndarray, y: np.ndarray, metric: str) -> float:
    """Fast unweighted null-model VPC (matches R diagnostics.R); PCV adds an
    additive main-effects model at the stratum level via weighted least
    squares on the empirical logits (the design covariates are recoverable
    from the stratum index)."""
    e, n = _stratum_stats(sid, y)
    if len(n) < 2:
        return float("nan")
    p = (e + 0.5) / (n + 1.0)
    lg = np.log(p / (1.0 - p))
    sv = 1.0 / (e + 0.5) + 1.0 / (n - e + 0.5)
    p0 = float(np.clip(e.sum() / n.sum(), 1e-6, 1 - 1e-6))
    null_pred = math.log(p0 / (1.0 - p0))
    var0 = max(float(np.var(lg - null_pred, ddof=1)) - float(np.mean(sv)), 0.0)
    if metric == "vpc":
        return vpc_latent(var0)

    # PCV: residual variance after additive main effects.
    present = np.unique(sid)
    tab = _TAB[present]
    cols = [np.ones(len(tab))]
    for j in range(4):
        vals = np.unique(tab[:, j])
        for lv in vals[1:]:
            cols.append((tab[:, j] == lv).astype(float))
    X = np.column_stack(cols)
    w = 1.0 / sv
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], lg * sw, rcond=None)
    main_pred = X @ beta
    var1 = max(float(np.var(lg - main_pred, ddof=1)) - float(np.mean(sv)), 0.0)
    if var0 <= 0:
        return float("nan")
    return 100.0 * (var0 - var1) / var0


def _spec_from_bits(bits: Sequence[int], a: CohortSpec, b: CohortSpec) -> CohortSpec:
    """bit=0 -> channel takes cohort A's value, bit=1 -> cohort B's value.
    Channel order matches CHANNELS."""
    pick = {c: (b if bit else a) for c, bit in zip(CHANNELS, bits)}
    return CohortSpec(
        prevalence_shift=pick["prevalence"].prevalence_shift,
        interaction_sd=pick["structure"].interaction_sd,
        structure_seed=pick["structure"].structure_seed,
        sparse=pick["sparsity"].sparse,
        detection_strength=pick["detection"].detection_strength,
        attrition_strength=pick["attrition"].attrition_strength,
        n=b.n if bits[0] else a.n,  # sample size follows the prevalence/base cohort
    )


def _value(bits, a, b, metric, n_rep, seed0) -> float:
    spec = _spec_from_bits(bits, a, b)
    vals = []
    for r in range(n_rep):
        sid, y, _ = _simulate(spec, seed0 + r)
        vals.append(_metric(sid, y, metric=metric))
    vals = [v for v in vals if not math.isnan(v)]
    if not vals:
        # PCV is genuinely undefined when the null variance collapses to zero
        # (Study 4's PCV-collapse regime); treat as NaN so callers can flag it.
        return float("nan")
    return float(np.mean(vals))


def _shapley(a, b, metric, n_rep, seed0) -> Dict[str, float]:
    m = len(CHANNELS)
    cache = {bits: _value(bits, a, b, metric, n_rep, seed0)
             for bits in itertools.product((0, 1), repeat=m)}
    phi = {}
    for i, ch in enumerate(CHANNELS):
        total = 0.0
        for bits in itertools.product((0, 1), repeat=m):
            if bits[i] == 1:
                continue
            s = sum(bits)
            weight = factorial(s) * factorial(m - s - 1) / factorial(m)
            b1 = list(bits)
            b1[i] = 1
            total += weight * (cache[tuple(b1)] - cache[bits])
        phi[ch] = total
    phi["_vpc_a" if metric == "vpc" else "_pcv_a"] = cache[tuple([0] * m)]
    phi["_vpc_b" if metric == "vpc" else "_pcv_b"] = cache[tuple([1] * m)]
    return phi


def true_structural_gap(a: CohortSpec, b: CohortSpec, metric: str = "vpc") -> float:
    """The exact structural gap known by construction (no simulation): the
    difference in the population functional of the true stratum logits. Used
    to self-validate the decomposition."""
    eta_a = _true_eta(a.prevalence_shift, a.interaction_sd, a.structure_seed)
    eta_b = _true_eta(b.prevalence_shift, b.interaction_sd, b.structure_seed)
    if metric == "vpc":
        return vpc_latent(float(np.var(eta_b, ddof=1))) - vpc_latent(float(np.var(eta_a, ddof=1)))
    raise ValueError("true_structural_gap currently supports metric='vpc' only")


def cross_cohort_decomposition(
    spec_a: CohortSpec,
    spec_b: CohortSpec,
    metric: str = "vpc",
    n_rep: int = 12,
    n_boot: int = 40,
    level: float = 0.95,
    seed: int = 100,
) -> Dict[str, object]:
    """Decompose the observed cohort B − cohort A gap in ``metric`` ("vpc" or
    "pcv") into a Shapley share for each of the four artefact channels
    (prevalence, sparsity, detection, attrition) plus a structural residual.

    Shares sum exactly to ``value_b - value_a``. A bootstrap over independent
    simulation streams gives a confidence interval on the **structural** share
    — the inferential answer to "is the residual structural difference
    distinguishable from zero, or could nuisance differences alone explain the
    observed gap?"

    Returns a dict with per-channel shares, ``value_a``/``value_b``/``gap``,
    ``structure_share``, ``structure_ci`` (level CI), ``artefact_total``, and
    (for metric="vpc") ``true_structural_gap`` for self-validation.
    """
    if metric not in ("vpc", "pcv"):
        raise ValueError("metric must be 'vpc' or 'pcv'")
    if not (0.0 < level < 1.0):
        raise ValueError("level must be in (0, 1)")

    phi = _shapley(spec_a, spec_b, metric, n_rep, seed)
    va = phi.pop(f"_{metric}_a")
    vb = phi.pop(f"_{metric}_b")

    # Bootstrap the structural share across independent simulation streams.
    boot = []
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        s = int(rng.integers(1, 2**31 - 1))
        val = _shapley(spec_a, spec_b, metric, max(2, n_rep // 3), s)["structure"]
        if not math.isnan(val):
            boot.append(val)
    alpha = 1.0 - level
    # PCV can be undefined in nuisance-heavy subsets (var_null -> 0); if too
    # many bootstrap replicates are undefined the CI is unreliable.
    reliable = len(boot) >= max(5, int(0.5 * n_boot))
    if reliable:
        lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
        distinguishable: Optional[bool] = not (lo <= 0.0 <= hi)
    else:
        lo = hi = float("nan")
        distinguishable = None

    out: Dict[str, object] = {
        "metric": metric,
        "value_a": va,
        "value_b": vb,
        "gap": vb - va,
        "shares": {c: phi[c] for c in CHANNELS},
        "structure_share": phi["structure"],
        "structure_ci": (float(lo), float(hi)),
        "artefact_total": sum(phi[c] for c in _ARTEFACT_CHANNELS),
        "level": level,
        "structure_distinguishable_from_zero": distinguishable,
        "reliable": reliable,
    }
    if metric == "vpc":
        out["true_structural_gap"] = true_structural_gap(spec_a, spec_b, "vpc")
    return out
