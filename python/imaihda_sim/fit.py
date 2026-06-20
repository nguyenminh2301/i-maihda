from __future__ import annotations

import math
from typing import Dict
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Binomial

LOGISTIC_L1_VARIANCE = math.pi ** 2 / 3.0


def vpc_latent(stratum_variance: float) -> float:
    """Latent-response VPC for logistic MAIHDA."""
    if stratum_variance < 0:
        raise ValueError("stratum_variance must be non-negative")
    return 100.0 * stratum_variance / (stratum_variance + LOGISTIC_L1_VARIANCE)


def pcv(var_null: float, var_main: float) -> float:
    """Proportional change in stratum-level variance from null to main-effects model."""
    if var_null <= 0:
        return float("nan")
    return 100.0 * (var_null - var_main) / var_null


def _aggregate_strata(df: pd.DataFrame) -> pd.DataFrame:
    g = (
        df.groupby(["stratum", "sex", "education", "wealth", "rural"], observed=True)
        .agg(events=("y", "sum"), n=("y", "size"))
        .reset_index()
    )
    g["empirical_p"] = (g["events"] + 0.5) / (g["n"] + 1.0)
    g["empirical_logit"] = np.log(g["empirical_p"] / (1.0 - g["empirical_p"]))
    g["logit_sampling_var"] = 1.0 / (g["events"] + 0.5) + 1.0 / (g["n"] - g["events"] + 0.5)
    g["precision_weight"] = 1.0 / g["logit_sampling_var"]
    return g


def _weighted_variance(x: np.ndarray, w: np.ndarray) -> float:
    mu = np.average(x, weights=w)
    return float(np.average((x - mu) ** 2, weights=w))


def _between_stratum_variance(logit: np.ndarray, pred: np.ndarray, sampling_var: np.ndarray, weights: np.ndarray) -> float:
    """Fast diagnostic estimate of residual between-stratum logit variance.

    This is a synthetic-data diagnostic approximation, not a replacement for the
    full random-intercept GLMM used in empirical I-MAIHDA. It is intentionally
    fast enough for repeated tests and CI-style checks.
    """
    residual = logit - pred
    raw_var = _weighted_variance(residual, weights)
    expected_noise = float(np.average(sampling_var, weights=weights))
    return max(raw_var - expected_noise, 0.0)


def _stratum_diagnostics(df: pd.DataFrame, strata: pd.DataFrame) -> Dict[str, float]:
    prev = df.groupby("stratum", observed=True)["y"].mean()
    return {
        "n": int(len(df)),
        "n_strata": int(df["stratum"].nunique()),
        "min_stratum_n": int(strata["n"].min()),
        "median_stratum_n": float(strata["n"].median()),
        "max_stratum_n": int(strata["n"].max()),
        "overall_prevalence": float(df["y"].mean()),
        "min_stratum_prevalence": float(prev.min()),
        "max_stratum_prevalence": float(prev.max()),
        "true_prevalence": float(df["y_true"].mean()) if "y_true" in df else float("nan"),
        "mean_detection_probability": float(df["detection_probability"].mean())
        if "detection_probability" in df
        else float("nan"),
    }


def fit_imaihda(df: pd.DataFrame) -> Dict[str, float | str]:
    """Fit fast synthetic I-MAIHDA diagnostics.

    The function estimates the null and main-effects between-stratum variance
    from empirical stratum logits after subtracting binomial sampling variance.
    It is designed for simulation and benchmarking only. The PhD empirical work
    should use full random-intercept MAIHDA models.
    """
    strata = _aggregate_strata(df)
    logit = strata["empirical_logit"].to_numpy()
    sampling_var = strata["logit_sampling_var"].to_numpy()
    weights = strata["precision_weight"].to_numpy()

    p_overall = min(max(float(df["y"].mean()), 1e-6), 1 - 1e-6)
    null_pred = np.repeat(math.log(p_overall / (1.0 - p_overall)), len(strata))
    var0 = _between_stratum_variance(logit, null_pred, sampling_var, weights)

    main = smf.glm(
        "y ~ C(sex) + C(education) + C(wealth) + C(rural)",
        data=df,
        family=Binomial(),
    ).fit()
    main_pred = np.asarray(main.predict(strata, which="linear"))
    var1 = _between_stratum_variance(logit, main_pred, sampling_var, weights)

    diagnostics = _stratum_diagnostics(df, strata)
    out: Dict[str, float | str] = {
        **diagnostics,
        "var_null": var0,
        "vpc_null": vpc_latent(var0),
        "var_main": var1,
        "vpc_main": vpc_latent(var1),
        "pcv": pcv(var0, var1),
        "estimator": "fast_empirical_logit_diagnostic",
        "warnings": "diagnostic approximation; use full random-intercept MAIHDA for empirical analysis",
    }
    return out


def fit_scenario(name: str, scenario) -> Dict[str, float | str]:
    from .simulate import simulate_intersectional_data

    df = simulate_intersectional_data(
        n=scenario.n,
        prevalence_shift=scenario.prevalence_shift,
        interaction_sd=scenario.interaction_sd,
        detection_strength=scenario.detection_strength,
        sparse=scenario.sparse,
        seed=scenario.seed,
    )
    res = fit_imaihda(df)
    res.update(
        {
            "scenario": name,
            "description": scenario.description,
            "interaction_sd": scenario.interaction_sd,
            "detection_strength": scenario.detection_strength,
            "sparse": scenario.sparse,
        }
    )
    return res
