from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable
import numpy as np
import pandas as pd


def _logit_inv(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


@dataclass(frozen=True)
class Scenario:
    name: str
    description: str
    n: int = 6000
    prevalence_shift: float = -2.10
    interaction_sd: float = 0.0
    detection_strength: float = 0.0
    sparse: bool = False
    seed: int = 42


def _stratum_table() -> pd.DataFrame:
    rows = []
    for sex in [0, 1]:
        for education in [0, 1, 2]:          # 0 high, 2 low
            for wealth in [0, 1, 2]:         # 0 high, 2 low
                for rural in [0, 1]:         # 1 rural / less-resourced
                    rows.append((sex, education, wealth, rural))
    return pd.DataFrame(rows, columns=["sex", "education", "wealth", "rural"])


def simulate_intersectional_data(
    n: int = 6000,
    prevalence_shift: float = -2.10,
    interaction_sd: float = 0.0,
    detection_strength: float = 0.0,
    sparse: bool = False,
    seed: int = 42,
) -> pd.DataFrame:
    """Simulate individual-level binary outcomes nested in intersectional strata.

    The data are synthetic. They are designed to test how I-MAIHDA summary
    statistics react to additive social gradients, residual interaction, sparse
    strata, rare outcomes, and SES-patterned under-detection. They are not meant
    to represent any real cohort.
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

    # Additive data-generating component. Higher education/wealth values denote disadvantage.
    linear_predictor = (
        prevalence_shift
        + 0.20 * x["sex"].to_numpy()
        + 0.35 * x["education"].to_numpy()
        + 0.30 * x["wealth"].to_numpy()
        + 0.25 * x["rural"].to_numpy()
    )

    # True residual intersectional component not reducible to additive main effects.
    # This creates the benchmark where PCV should be lower after main effects are included.
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

    # Differential detection: true cases are less likely to be recorded when education/wealth
    # disadvantage is high. This can attenuate, exaggerate, or redistribute observed inequality.
    if detection_strength > 0:
        detect_p = _logit_inv(
            2.0
            - detection_strength * x["education"].to_numpy()
            - detection_strength * x["wealth"].to_numpy()
            - 0.40 * detection_strength * x["rural"].to_numpy()
        )
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


def scenario_grid(seed: int = 42) -> Dict[str, Scenario]:
    return {
        "A": Scenario(
            name="A",
            description="Additive social gradient, equal detection",
            seed=seed,
        ),
        "B": Scenario(
            name="B",
            description="Residual intersectional heterogeneity, equal detection",
            interaction_sd=0.90,
            seed=seed,
        ),
        "C": Scenario(
            name="C",
            description="Additive structure with SES-patterned under-detection",
            detection_strength=0.80,
            seed=seed,
        ),
        "D": Scenario(
            name="D",
            description="Residual interaction plus SES-patterned under-detection",
            interaction_sd=0.90,
            detection_strength=0.80,
            seed=seed,
        ),
        "E": Scenario(
            name="E",
            description="Residual interaction with rare outcome and sparse strata",
            n=3500,
            prevalence_shift=-3.00,
            interaction_sd=0.90,
            sparse=True,
            seed=seed,
        ),
    }
