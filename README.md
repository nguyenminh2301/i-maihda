# imaihda-hic-mic-simulation v3.1

Synthetic-data workflow for stress-testing one narrow methodological issue in a proposed HIC-MIC I-MAIHDA PhD project: whether VPC and PCV can be compared across aging contexts without diagnostics for prevalence, sparse strata, and differential detection.

This repository does **not** use HRS, ELSA, SHARE, CHARLS, UK Biobank, Vietnam, or any other real data. It makes no empirical claim about any population. It is a reproducible demonstration of a methodological caveat.

## Question addressed

If a middle-income cohort shows higher VPC or lower PCV than high-income cohorts, does that necessarily mean the intersectional structure of inequality is different? This demo shows the answer is no. VPC and PCV are informative, but they can also move with outcome prevalence, sparse strata, and SES-patterned under-detection.

## What changed in v3.1

- Replaced the previous untested R-only workflow with an executable Python implementation, because the current runtime has Python but not R.
- Added scenario benchmarks with pass/fail criteria.
- Added a detection-bias sweep.
- Added stratum-size and prevalence diagnostics beside every VPC/PCV estimate.
- Added figures and run metadata.
- Added pytest tests and a final run log.

## Method

The demo simulates individuals nested in 36 intersectional strata defined by sex, education, wealth, and rural/less-resourced setting. It computes fast synthetic I-MAIHDA diagnostics using empirical stratum logits and a main-effects logistic model. This is a simulation diagnostic, not a substitute for full random-intercept MAIHDA in empirical work:

1. Null intersectional diagnostic: residual between-stratum logit variance around the overall logit
2. Main-effects diagnostic: residual between-stratum logit variance after additive main effects

For a binary outcome, VPC is calculated on the latent logistic scale:

`VPC = sigma2_stratum / (sigma2_stratum + pi^2/3) * 100`

PCV is calculated as:

`PCV = (sigma2_null - sigma2_main) / sigma2_null * 100`

## Scenarios

| Scenario | Purpose |
|---|---|
| A | Additive social gradient, equal detection |
| B | Residual intersectional heterogeneity, equal detection |
| C | Additive structure with SES-patterned under-detection |
| D | Residual interaction plus SES-patterned under-detection |
| E | Residual interaction with rare outcome and sparse strata |

## Reproduce

```bash
python -m pytest -q
python scripts/run_all.py
```

Outputs:

- `outputs/results.csv`: scenario results
- `outputs/detection_sweep.csv`: detection-bias sweep
- `outputs/benchmark_checks.csv`: pass/fail checks
- `figures/scenario_vpc_pcv.png`: VPC-PCV scenario plot
- `figures/detection_sweep.png`: detection bias plot
- `logs/run_log.md`: final run log
- `outputs/run_metadata.json`: runtime metadata

## Interpretation of a successful run

A successful run should show that:

1. In the additive-only scenario, most between-stratum variance is removed by additive main effects.
2. True residual intersectional heterogeneity raises null-model VPC and leaves more residual variance after main effects.
3. SES-patterned under-detection lowers observed prevalence and can mask or redistribute VPC/PCV.
4. Sparse strata must be diagnosed before any HIC-MIC VPC/PCV comparison.

## Scope limits

This is not a new estimator and not a replacement for established I-MAIHDA workflows. It uses a fast empirical-logit diagnostic rather than full GLMM fitting, because the purpose is repeated stress testing rather than empirical estimation. For empirical work, estimates should be checked against the modelling strategy used by the target research group and should include sensitivity analyses for convergence, sparse strata, outcome prevalence, missingness, and differential detection.

## References

- Evans CR, Williams DR, Onnela JP, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM - Population Health*. 2018.
- O'Sullivan JL, Alonso-Perez E, et al. Onset of Type 2 diabetes in adults aged 50 and older in Europe: an intersectional multilevel analysis of individual heterogeneity and discriminatory accuracy. *Diabetology & Metabolic Syndrome*. 2024. doi:10.1186/s13098-024-01533-3.
- Elff M, Heisig JP, Schaeffer M, Shikano S. Multilevel analysis with few clusters. *British Journal of Political Science*. 2021;51:412-426.
