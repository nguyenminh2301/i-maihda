from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime, timezone
import platform
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from imaihda_sim import scenario_grid, simulate_intersectional_data, fit_imaihda
from imaihda_sim.fit import fit_scenario
from imaihda_sim.benchmarks import evaluate_benchmarks

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
FIG = ROOT / "figures"
LOG = ROOT / "logs"
CHK = ROOT / "checkpoints"
for p in [OUT, FIG, LOG, CHK]:
    p.mkdir(exist_ok=True)


def detection_sweep(seed: int = 42) -> pd.DataFrame:
    rows = []
    for strength in [0.0, 0.3, 0.6, 0.9, 1.2]:
        df = simulate_intersectional_data(interaction_sd=0.90, detection_strength=strength, seed=seed)
        res = fit_imaihda(df)
        res["detection_strength"] = strength
        rows.append(res)
    return pd.DataFrame(rows)


def plot_results(results: pd.DataFrame, sweep: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.scatter(results["vpc_null"], results["pcv"], s=70)
    for _, row in results.iterrows():
        ax.annotate(row["scenario"], (row["vpc_null"], row["pcv"]), xytext=(5, 4), textcoords="offset points")
    ax.set_xlabel("Null-model VPC (%)")
    ax.set_ylabel("PCV after additive main effects (%)")
    ax.set_title("Synthetic I-MAIHDA scenarios: VPC/PCV behaviour")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "scenario_vpc_pcv.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(sweep["detection_strength"], sweep["vpc_null"], marker="o", label="VPC null")
    ax2 = ax.twinx()
    ax2.plot(sweep["detection_strength"], 100 * sweep["overall_prevalence"], marker="s", linestyle="--", label="Observed prevalence")
    ax.set_xlabel("SES-patterned under-detection strength")
    ax.set_ylabel("VPC null (%)")
    ax2.set_ylabel("Observed prevalence (%)")
    ax.set_title("Detection bias can mask residual intersectional heterogeneity")
    fig.tight_layout()
    fig.savefig(FIG / "detection_sweep.png", dpi=160)
    plt.close(fig)


def main() -> None:
    started = datetime.now(timezone.utc).isoformat()
    rows = []
    for name, scenario in scenario_grid().items():
        rows.append(fit_scenario(name, scenario))
    results = pd.DataFrame(rows)
    results = results[
        [
            "scenario", "description", "n", "n_strata", "min_stratum_n", "median_stratum_n",
            "max_stratum_n", "overall_prevalence", "true_prevalence", "mean_detection_probability",
            "var_null", "vpc_null", "var_main", "vpc_main", "pcv", "interaction_sd",
            "detection_strength", "sparse", "warnings"
        ]
    ]
    sweep = detection_sweep()
    benchmarks = evaluate_benchmarks(results)

    results.to_csv(OUT / "results.csv", index=False)
    sweep.to_csv(OUT / "detection_sweep.csv", index=False)
    benchmarks.to_csv(OUT / "benchmark_checks.csv", index=False)
    plot_results(results, sweep)

    metadata = {
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "all_benchmarks_passed": bool(benchmarks["passed"].all()),
    }
    (OUT / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    log = [
        "# Run log",
        f"Started: {metadata['started_utc']}",
        f"Finished: {metadata['finished_utc']}",
        f"Python: {metadata['python']}",
        "",
        "## Benchmark checks",
        benchmarks.to_markdown(index=False),
        "",
        "## Scenario results",
        results[["scenario", "overall_prevalence", "vpc_null", "vpc_main", "pcv", "min_stratum_n"]].to_markdown(index=False),
    ]
    (LOG / "run_log.md").write_text("\n".join(log), encoding="utf-8")
    (CHK / "checkpoint_03_final_run.txt").write_text(
        "Final run completed. All benchmark checks passed: " + str(bool(benchmarks["passed"].all())) + "\n",
        encoding="utf-8",
    )
    if not benchmarks["passed"].all():
        raise SystemExit("One or more benchmark checks failed; inspect outputs/benchmark_checks.csv")


if __name__ == "__main__":
    main()
