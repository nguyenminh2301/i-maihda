"""Robustness/simulation study for sparse_strata_vpc() under violated
assumptions (Study 2 of the methods-note robustness study).

Runs two scenarios against an analytic ground truth (the true population
variance of the 36 stratum logits, computed directly from the generative
formula -- no simulation, no estimator bias):

  M2-1  True random effects deviate from the calibration model's i.i.d.
        Gaussian assumption: escalating structured spikes, and alternative
        tail shapes (heavy-tailed, right-skewed) at fixed dispersion.
  M2-2  Rare prevalence combined with sparse allocation (interacts with
        Laplace-smoothing behavior as p -> 0).

Prints a results table with pass/fail against the success criteria in
docs/METHODS_NOTE_ROBUSTNESS.md, writes figures/sparse_strata_robustness.png
and figures/sparse_strata_robustness_results.csv. Uses synthetic data only.
Run from the repo root:

    python scripts/validation/sparse_strata_robustness.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from imaihda_sim import simulate_intersectional_data, sparse_strata_vpc  # noqa: E402
from imaihda_sim.fit import vpc_latent  # noqa: E402
from imaihda_sim.robustness import simulate_structured_random_effects  # noqa: E402

N_REP = 30
FIG_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "sparse_strata_robustness.png")
CSV_OUT = os.path.join(os.path.dirname(__file__), "..", "..", "figures", "sparse_strata_robustness_results.csv")

NONINFERIORITY_MARGIN = 0.5  # pp
BIAS_REDUCTION_TARGET = 0.40
COVERAGE_TARGET = 0.80


def _analytic_truth_structured(prevalence_shift, severity, tail_family, interaction_sd, seed):
    """Exact population VPC of the true stratum logits for
    simulate_structured_random_effects -- mirrors its residual-generation
    logic exactly (same RNG calls) but operates on the full 36-row stratum
    table, independent of which strata happen to be sampled."""
    from imaihda_sim.simulate import _stratum_table

    strata = _stratum_table()
    k = len(strata)
    eta_add = (
        prevalence_shift
        + 0.20 * strata["sex"].to_numpy()
        + 0.35 * strata["education"].to_numpy()
        + 0.30 * strata["wealth"].to_numpy()
        + 0.25 * strata["rural"].to_numpy()
    )
    rng_re = np.random.default_rng(seed + 999)
    if tail_family == "gaussian":
        raw = rng_re.normal(0.0, interaction_sd, size=k)
    elif tail_family == "student_t3":
        raw = rng_re.standard_t(df=3, size=k) * (interaction_sd / np.sqrt(3.0))
    else:
        raw = rng_re.lognormal(mean=0.0, sigma=0.5, size=k)
        raw = (raw - raw.mean()) / raw.std() * interaction_sd

    if severity > 0:
        high_mag = 0.90 * severity
        low_mag = -0.60 * severity
        raw = raw + high_mag * (
            (strata["education"].to_numpy() == 2)
            & (strata["wealth"].to_numpy() == 2)
            & (strata["rural"].to_numpy() == 1)
        )
        raw = raw + low_mag * (
            (strata["education"].to_numpy() == 0)
            & (strata["wealth"].to_numpy() == 2)
            & (strata["rural"].to_numpy() == 0)
        )
    residual = raw - raw.mean()
    eta_true = eta_add + residual
    return vpc_latent(np.var(eta_true, ddof=1))


def _analytic_truth_plain(prevalence_shift, interaction_sd, seed):
    from imaihda_sim.simulate import _stratum_table

    strata = _stratum_table()
    k = len(strata)
    eta_add = (
        prevalence_shift
        + 0.20 * strata["sex"].to_numpy()
        + 0.35 * strata["education"].to_numpy()
        + 0.30 * strata["wealth"].to_numpy()
        + 0.25 * strata["rural"].to_numpy()
    )
    rng_re = np.random.default_rng(seed + 999)
    raw = rng_re.normal(0.0, interaction_sd, size=k)
    raw += 0.90 * ((strata["education"].to_numpy() == 2) & (strata["wealth"].to_numpy() == 2)
                   & (strata["rural"].to_numpy() == 1))
    raw -= 0.60 * ((strata["education"].to_numpy() == 0) & (strata["wealth"].to_numpy() == 2)
                   & (strata["rural"].to_numpy() == 0))
    residual = raw - raw.mean()
    eta_true = eta_add + residual
    return vpc_latent(np.var(eta_true, ddof=1))


def _bias_reduction(bias_naive, bias_corrected):
    if abs(bias_naive) < 1.0:
        return None
    return 1.0 - abs(bias_corrected) / abs(bias_naive)


def _classify(bias_naive, bias_corrected, coverage, require_full_pass):
    noninferior = abs(bias_corrected) <= abs(bias_naive) + NONINFERIORITY_MARGIN
    reduction = _bias_reduction(bias_naive, bias_corrected)
    if not noninferior:
        return "BREAKDOWN"
    if require_full_pass:
        ok_reduction = reduction is not None and reduction >= BIAS_REDUCTION_TARGET
        ok_coverage = coverage is None or coverage >= COVERAGE_TARGET
        if ok_reduction and ok_coverage:
            return "PASS"
        return "FAIL (below bias-reduction or coverage target)"
    if reduction is not None and reduction >= BIAS_REDUCTION_TARGET:
        return "attenuated-but-positive (exceeds mild-target anyway)"
    return "attenuated-but-net-neutral-or-positive"


def run_m2_1():
    rows = []
    # Spike-magnitude escalation (tail_family fixed at gaussian).
    for severity in (0, 1, 2, 3):
        for seed in range(N_REP):
            df = simulate_structured_random_effects(
                n=3500, prevalence_shift=-3.00, severity=severity, tail_family="gaussian",
                interaction_sd=0.15, sparse=True, seed=seed,
            )
            truth = _analytic_truth_structured(-3.00, severity, "gaussian", 0.15, seed)
            r = sparse_strata_vpc(df, seed=1000 + seed)
            rows.append(dict(scenario="M2-1", arm=f"severity={severity}", seed=seed,
                              vpc_true=truth, vpc_naive=r["vpc_null_naive"],
                              vpc_corrected=r["vpc_null_corrected"],
                              ci_lower=r["ci_lower"], ci_upper=r["ci_upper"]))
    # Tail-shape variation (severity fixed at the baseline=1 spike pattern).
    for tail_family in ("student_t3", "skew_lognormal"):
        for seed in range(N_REP):
            df = simulate_structured_random_effects(
                n=3500, prevalence_shift=-3.00, severity=1, tail_family=tail_family,
                interaction_sd=0.15, sparse=True, seed=seed,
            )
            truth = _analytic_truth_structured(-3.00, 1, tail_family, 0.15, seed)
            r = sparse_strata_vpc(df, seed=1000 + seed)
            rows.append(dict(scenario="M2-1", arm=f"tail={tail_family}", seed=seed,
                              vpc_true=truth, vpc_naive=r["vpc_null_naive"],
                              vpc_corrected=r["vpc_null_corrected"],
                              ci_lower=r["ci_lower"], ci_upper=r["ci_upper"]))
    return pd.DataFrame(rows)


def run_m2_2():
    rows = []
    diag = []
    for prevalence_shift in (-2.10, -3.00, -4.00, -5.00):
        min_ns, zero_event_flags = [], []
        for seed in range(N_REP):
            df = simulate_intersectional_data(
                n=3500, prevalence_shift=prevalence_shift, interaction_sd=0.15,
                sparse=True, seed=seed,
            )
            truth = _analytic_truth_plain(prevalence_shift, 0.15, seed)
            r = sparse_strata_vpc(df, seed=2000 + seed)
            rows.append(dict(scenario="M2-2", arm=f"prevalence_shift={prevalence_shift}", seed=seed,
                              vpc_true=truth, vpc_naive=r["vpc_null_naive"],
                              vpc_corrected=r["vpc_null_corrected"],
                              ci_lower=r["ci_lower"], ci_upper=r["ci_upper"]))
            min_ns.append(r["min_stratum_n"])
            strata_counts = df.groupby("stratum", observed=True)["y"].agg(["sum", "size"])
            zero_event_flags.append(bool((strata_counts["sum"] == 0).any()))
        diag.append(dict(prevalence_shift=prevalence_shift, mean_min_stratum_n=np.mean(min_ns),
                          frac_with_zero_event_stratum=np.mean(zero_event_flags)))
    return pd.DataFrame(rows), pd.DataFrame(diag)


def summarize(df, require_full_pass_arms):
    out = []
    for (scenario, arm), g in df.groupby(["scenario", "arm"], sort=False):
        bias_naive = float((g["vpc_naive"] - g["vpc_true"]).mean())
        bias_corrected = float((g["vpc_corrected"] - g["vpc_true"]).mean())
        coverage = float(((g["ci_lower"] <= g["vpc_true"]) & (g["vpc_true"] <= g["ci_upper"])).mean())
        reduction = _bias_reduction(bias_naive, bias_corrected)
        classification = _classify(bias_naive, bias_corrected, coverage, arm in require_full_pass_arms)
        out.append(dict(
            scenario=scenario, arm=arm, n_rep=len(g),
            mean_vpc_true=round(float(g["vpc_true"].mean()), 2),
            bias_naive=round(bias_naive, 2),
            bias_corrected=round(bias_corrected, 2),
            bias_reduction_pct=None if reduction is None else round(100 * reduction, 1),
            coverage=round(coverage, 2),
            classification=classification,
        ))
    return pd.DataFrame(out)


def main():
    m2_1 = run_m2_1()
    m2_2, diag = run_m2_2()

    require_full_pass = {"severity=0", "severity=1", "prevalence_shift=-2.1", "prevalence_shift=-3.0"}
    summary = pd.concat([summarize(m2_1, require_full_pass), summarize(m2_2, require_full_pass)],
                         ignore_index=True)

    pd.set_option("display.width", 140)
    print(summary.to_string(index=False))
    print("\nM2-2 sparsity/rarity diagnostics:")
    print(diag.to_string(index=False))

    os.makedirs(os.path.dirname(CSV_OUT), exist_ok=True)
    summary.to_csv(CSV_OUT, index=False)
    print(f"\nwrote {os.path.abspath(CSV_OUT)}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    sev_summary = summarize(m2_1[m2_1["arm"].str.startswith("severity=")], require_full_pass)
    sev = [int(a.split("=")[1]) for a in sev_summary["arm"]]
    ax.plot(sev, sev_summary["bias_naive"], "o-", color="gray", label="naive")
    ax.plot(sev, sev_summary["bias_corrected"], "o-", color="#440154", label="corrected")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("spike severity (0=pure Gaussian, 1=baseline Scenario E)")
    ax.set_ylabel("Bias in VPC (pp)")
    ax.set_title("M2-1: non-Gaussian random effects")
    ax.legend(fontsize=8)

    ax = axes[1]
    prev_summary = summarize(m2_2, require_full_pass)
    prevs = [float(a.split("=")[1]) for a in prev_summary["arm"]]
    ax.plot(prevs, prev_summary["bias_naive"], "o-", color="gray", label="naive")
    ax.plot(prevs, prev_summary["bias_corrected"], "o-", color="#440154", label="corrected")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("prevalence_shift (more negative = rarer outcome)")
    ax.set_ylabel("Bias in VPC (pp)")
    ax.set_title("M2-2: rare prevalence x sparsity")
    ax.legend(fontsize=8)
    ax.invert_xaxis()

    fig.suptitle("Sparse-strata bias correction under assumption violations", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.abspath(FIG_OUT), dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.abspath(FIG_OUT)}")


if __name__ == "__main__":
    main()
