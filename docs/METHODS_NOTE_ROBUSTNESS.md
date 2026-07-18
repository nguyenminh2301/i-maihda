# Robustness of Detection-Bias and Sparse-Strata Corrections to Violated Assumptions: A Simulation Study

*Prepared as groundwork for a methods + simulation manuscript. No real data — synthetic simulation only, as throughout this repository.*

## 1. Purpose and scope

`imaihda` provides two quantitative-bias-analysis (QBA) tools with no equivalent in CRAN `MAIHDA`: `correct_detection_bias()` (sensitivity correction for SES-patterned under-detection of the outcome) and `sparse_strata_vpc()` (bias correction and confidence intervals for the fast VPC estimator in sparse intersectional strata). Both were validated under the exact assumptions their correction procedure encodes — the analyst's detection `score=` formula matches the true mechanism; true random effects are i.i.d. Gaussian. Real applications will not satisfy these assumptions exactly.

This note reports a systematic simulation study of what happens when each method's core assumptions are violated: does the correction still help, does it merely stop helping (graceful degradation), or does it actively make estimates worse than doing nothing (breakdown)? All findings below are from code that was actually executed (`scripts/validation/detection_bias_robustness.py`, `scripts/validation/sparse_strata_robustness.py`, `scripts/validation/joint_robustness_capstone.py`); no numbers here are hypothetical.

## 2. Background

**Detection-bias correction.** The observed outcome is `y = y_true × detected`, with `logit(P(detect)) = baseline_logit − δ·score`. Given an analyst-specified `score` (default `education + wealth + 0.4·rural`) and sensitivity parameter `δ`, `correct_detection_bias()` inverts the stratum-level counts (`e_true_j = min(n_j, e_obs_j/d_j(δ))`) and refits the null and additive main-effects models on the corrected logits. `δ = 0` is the identity (no correction, matches observed VPC exactly).

**Sparse-strata correction.** The fast null-model VPC estimator smooths each stratum's empirical logit with a Laplace prior `(events + 0.5)/(n + 1)`, which shrinks the observed between-stratum spread below its true value when strata are sparse. `sparse_strata_vpc()` corrects this by *calibration*: simulating the estimator's own expected value across a grid of candidate true variances (i.i.d. Gaussian random effects, holding the observed stratum sizes fixed), then inverting that curve to find the variance consistent with the observed naive estimate.

## 3. Success criteria

For a scenario `s` with `n_rep` replicates and true value `θ_true`:

- `Bias(s) = mean(θ̂) − θ_true` (percentage points, VPC scale); `BiasReduction(s) = 1 − |Bias_corrected(s)|/|Bias_naive(s)|` (reported as N/A when `|Bias_naive(s)| < 1.0`pp — the ratio is uninformative when the naive estimator is already nearly unbiased).
- `Coverage(s)` = fraction of 95% CIs containing `θ_true` (`sparse_strata_vpc()` only; `correct_detection_bias()` returns no CI).

**Mild-to-moderate violation** (the smallest 1–2 severity levels tested per scenario):
1. **Non-inferiority**: `|Bias_corrected(s)| ≤ |Bias_naive(s)| + 0.5`pp.
2. **Minimum bias reduction**: `BiasReduction(s) ≥ 40%` (calibrated against the pre-existing, correctly-specified validation of `sparse_strata_vpc()`, which achieved ~93% reduction — 40% leaves deliberate headroom for degradation under misspecification while still requiring a materially useful correction).
3. **CI calibration** (Method 2 only): `Coverage(s) ≥ 80%` (a deliberately looser floor than nominal 95%, matching the existing test suite's Monte-Carlo-noise allowance).

**Severe violation** (the largest severity level tested per scenario):
4. Criterion 1 (non-inferiority) must **still hold** — the line between "graceful degradation" and "breakdown."
5. Bias reduction is not required to meet the 40% target; classified as **attenuated-but-net-neutral-or-positive** if `0% ≤ BiasReduction(s) < 40%`, or as **BREAKDOWN** if criterion 1 fails, with the mechanism stated explicitly.

**Non-goal.** This study does not claim nominal-coverage or near-zero-bias performance under severe or adversarial misspecification. It establishes bounds on harm and conditions for benefit — an honest scope statement, not a universal guarantee, in the spirit of quantitative bias analysis and E-value-style sensitivity frameworks (VanderWeele & Ding, 2017; see References).

## 4. Study 1: Detection-bias correction under assumption violations

Fixed scenario throughout: `n=12000`, `interaction_sd=0.9`, true `detection_strength (δ) = 0.8`, 20 replicates (seeds 0–19). Mean true VPC (on `y_true`) = **19.63%**; naive (uncorrected observed) bias = **−7.55**pp in every arm that holds the true DGP fixed.

### 4.1 M1-1 — Score functional-form misspecification

The analyst's `score=` passed to `correct_detection_bias()` differs from the true detection-driving covariates (`education + wealth + 0.4·rural`).

| Arm | Bias (naive) | Bias (corrected) | Bias reduction | Classification |
|---|---:|---:|---:|---|
| correctly specified | −7.55 | −1.31 | 82.6% | (reference) |
| `omit_rural` (drops rural term) | −7.55 | −3.64 | 51.7% | **PASS** |
| `wrong_covariate` (sex swapped for education) | −7.55 | −5.63 | 25.4% | attenuated-but-net-neutral-or-positive |
| `quadratic` (wrong functional form) | −7.55 | **+9.73** | −28.8% | **BREAKDOWN** |

Omitting a true covariate (`omit_rural`) degrades but still passes the mild-violation bar. Swapping in an irrelevant covariate (`wrong_covariate`) still helps somewhat (non-inferior) but falls short of the 40% target. **Getting the functional *form* wrong (`quadratic`) is qualitatively different**: the correction overshoots past the truth to the opposite sign and violates non-inferiority — the score's nonlinear spread over-corrects the most-disadvantaged strata disproportionately. This is the single clearest finding of Study 1: robustness to *which* covariates are included is much better than robustness to the *shape* of the score function.

### 4.2 M1-2 — Detection depends on the true outcome (identification breakdown)

`simulate_severity_dependent_detection()` adds `−severity_weight · y_true` to the detection logit: true cases are additionally under-detected *regardless of SES*, a classic verification-bias mechanism not identifiable from stratum covariates alone. Score and `δ` are correctly specified (best case for the method).

| `severity_weight` | Bias (naive) | Bias (corrected) | Bias reduction | Classification |
|---:|---:|---:|---:|---|
| 0.0 (regression check) | −7.55 | −1.31 | 82.6% | PASS |
| 0.3 | −7.22 | −3.93 | 45.5% | PASS |
| 0.6 | −7.18 | −6.01 | 16.3% | attenuated-but-net-neutral-or-positive |
| 0.8 | −7.14 | −7.17 | −0.4% | attenuated-but-net-neutral-or-positive (at the edge) |
| 1.0 | −7.04 | −7.79 | −10.7% | **BREAKDOWN** |

Degradation is smooth and monotonic, crossing from net-benefit to net-harm at **`severity_weight ≈ 0.8–0.9`**. This precisely locates the assumption's breaking point: the correction remains useful as long as outcome-dependent detection is a *minority* contributor relative to the SES-patterned mechanism it was designed to correct, and stops helping once it becomes comparably strong. This is the intended, honestly-characterized failure mode — no `score=` choice can fix it, because the violation is at the level of *what detection depends on* (an unobserved outcome), not *which covariates* the analyst picked.

### 4.3 M1-3 — `baseline_logit` (detection curvature) misspecification

Because detection is normalized relative to the `score = 0` stratum, `δ = 0` is the identity regardless of `baseline_logit`; this tests only curvature.

| `baseline_logit` | Bias (corrected) | Bias reduction |
|---:|---:|---:|
| 1.0 | +4.10 | 45.7% |
| 1.5 | +2.46 | 67.5% |
| 2.0 (true value) | −1.31 | 82.6% |
| 2.5 | −4.97 | 34.2% |
| 3.0 | −6.90 | 8.7% |

Degrades symmetrically around the true value and never approaches non-inferiority failure across the tested range — confirming this nuisance parameter is far less consequential than the score's functional form (§4.1).

### 4.4 Study 1 summary

| Criterion | Result |
|---|---|
| Non-inferiority holds under mild-to-moderate violation | Yes, in all tested arms |
| Non-inferiority holds under severe violation | **No** — `quadratic` score (M1-1) and `severity_weight=1.0` (M1-2) both breach it |
| ≥40% bias reduction under mild violation | Met for `omit_rural`, `severity_weight∈{0,0.3}`, `baseline_logit∈{1.5,2.0}`; not met for `wrong_covariate` |

## 5. Study 2: Sparse-strata calibration under assumption violations

Fixed scenario: `sparse=True`, `n=3500`, 30 replicates. Ground truth computed analytically (exact population variance of the 36 true stratum logits — no simulation needed for the truth itself).

### 5.1 M2-1 — Non-Gaussian true random effects

`simulate_structured_random_effects()` varies (a) spike-magnitude severity (0 = pure Gaussian, 1 = matches the package's original Scenario E exactly, 2–3 = escalating magnitude) and (b) base-noise shape at fixed dispersion (`student_t3` heavy-tailed, `skew_lognormal` right-skewed) — both violating the calibration model's i.i.d.-Gaussian assumption.

| Arm | True VPC | Bias (naive) | Bias (corrected) | Bias reduction | Coverage | Classification |
|---|---:|---:|---:|---:|---:|---|
| `severity=0` (pure Gaussian) | 5.60 | −3.02 | −0.77 | 74.4% | 0.97 | **PASS** |
| `severity=1` (baseline) | 9.69 | −5.11 | −1.35 | 73.6% | 0.80 | **PASS** |
| `severity=2` | 16.37 | −6.71 | +0.20 | 97.0% | 0.87 | attenuated-but-positive |
| `severity=3` | 24.51 | −6.94 | +1.29 | 81.5% | 0.93 | attenuated-but-positive |
| `tail=student_t3` | 9.61 | −5.41 | −1.77 | 67.2% | 0.83 | attenuated-but-positive |
| `tail=skew_lognormal` | 9.68 | −5.02 | −1.37 | 72.6% | 0.87 | attenuated-but-positive |

**Non-inferiority held in every arm** — the strongest robustness result in this study. Bias reduction did not degrade with escalating spike severity; if anything it *improved* at severity 2–3, and heavy-tailed/skewed noise shapes still achieved 67–73% reduction with coverage well above the 80% floor. The Gaussian calibration model is far more forgiving of shape misspecification than Study 1's detection score is of functional-form misspecification.

### 5.2 M2-2 — Rare prevalence × sparsity

`prevalence_shift ∈ {−2.10, −3.00, −4.00, −5.00}` (intercept-only baseline prevalence ≈10.9%, 4.7%, 1.8%, 0.7% respectively, before the additive social gradient) combined with `sparse=True`. All arms had ≥1 zero-event stratum in every replicate (mean smallest-stratum size ≈1.4).

| `prevalence_shift` | True VPC | Bias (naive) | Bias (corrected) | Bias reduction | Coverage | Classification |
|---:|---:|---:|---:|---:|---:|---|
| −2.10 | 9.69 | −4.67 | −0.98 | 78.9% | 0.87 | **PASS** |
| −3.00 | 9.69 | −5.11 | −1.48 | 71.0% | 0.83 | **PASS** |
| −4.00 | 9.69 | −6.70 | −3.07 | 54.1% | 0.87 | attenuated-but-positive |
| −5.00 (rarest, ~0.7%) | 9.69 | −6.58 | −4.87 | 25.9% | 0.87 | attenuated-but-net-neutral-or-positive |

Non-inferiority held at every prevalence level tested, including the most extreme (~0.7% marginal prevalence, universal zero-event strata). Bias reduction degrades smoothly from ~79% to ~26% as the outcome becomes rarer, but coverage remains stable (83–87%) throughout — the CI stays honest even as the point correction's benefit shrinks.

### 5.3 Study 2 summary

| Criterion | Result |
|---|---|
| Non-inferiority holds under mild-to-moderate violation | Yes |
| Non-inferiority holds under severe violation | **Yes, in every tested arm** — no breakdown found for `sparse_strata_vpc()` in this study |
| ≥40% bias reduction / ≥80% coverage under mild violation | Met in all "mild" arms; exceeded in most "severe" arms too |

### 5.4 Study 3 — K-generalization (previously deferred, now completed)

`simulate_k_strata()` (in `robustness.py`) generalizes the factorial design to an arbitrary number of strata (i.i.d. Gaussian residuals, no spikes — isolating the K effect). K ∈ {8 (2×2×2), 36 (2×3×3×2), 108 (2×3×3×3×2)} at fixed `n=3500`, sparse allocation, 20 replicates each, analytic ground truth:

| K | mean min stratum n | frac. at truncation floor | Bias (naive) | Bias (corrected) | Coverage (bootstrap CI) | Coverage (test-inversion CI) | Non-inferiority |
|---:|---:|---:|---:|---:|---:|---:|---|
| 8 | 19.1 | 0.65 | +0.01 | −0.03 | 0.95 | 0.90 | PASS |
| 36 | 1.6 | 0.65 | −2.76 | +0.01 | 1.00 | 0.80 | PASS |
| 108 | 1.0 | 0.90 | −6.33 | −5.31 | **0.10** | **0.95** | PASS |

Two findings. (1) **Point non-inferiority holds at every K tested** — the calibration never makes things worse, and at K=36 it removes essentially all bias (99.6% reduction). At K=108 with only ~32 individuals per stratum, most of the information is simply gone: the naive statistic sits at its `max(0,·)` truncation floor in 90% of replicates, and calibration can no longer recover much (16% reduction). (2) That floor regime exposed a **defect in the original bootstrap CI**: when the observed statistic is exactly at the floor, resampling at the (zero) calibrated variance degenerates the interval to `[0, 0]` — coverage collapsed to 0.10. The fix, added as `sparse_strata_vpc(ci_method="test_inversion")`, reports the set of true variances under which the observed statistic is not extreme (Monte Carlo test inversion); at the floor this is an honest one-sided `[0, upper]` interval, restoring coverage to 0.95 at K=108. The default `ci_method="bootstrap"` is unchanged for continuity of the results published above; the `at_floor` flag in the output tells users when the test-inversion interval is the one to report. Reproduce: `scripts/validation/k_generalization_study.py`.

### 5.5 Study 4 — PCV robustness pass (previously deferred, now completed)

The PCV — the literature's standard evidence statistic for "intersectional effects" — is a *ratio* of two variances, and behaves very differently from the null-model VPC tracked above. Population truth is analytic (`V_null = Var(η_true)`, `V_main = Var(true residual)`); n_rep=20 per arm (`scripts/validation/pcv_robustness.py`). Three findings, two of them honest negatives:

| Arm | mean true PCV | Bias naive | Bias corrected |
|---|---:|---:|---:|
| detection, correctly specified score/δ | 22.5 | **+3.7** | **+15.5** |
| detection, `omit_rural` score | 22.5 | +3.7 | +11.4 |
| detection, `quadratic` score | 22.5 | +3.7 | +26.9 |
| outcome-dependent detection, `severity_weight=1.0` | 22.5 | +15.3 | +3.1 |
| pure sparsity (n=3500, `sparse=True`), no correction exists | 22.3 | **+44.6** | — |

1. **Naive PCV is surprisingly robust to covariate-only detection bias** (+3.7pp): the ratio's numerator and denominator are distorted in partially cancelling ways.
2. **The detection correction *degrades* the PCV even under correct specification** (+15.5 vs +3.7): the correction's stratum-level WLS main-effects model is not the same estimator as the individual-level GLM the naive PCV uses, and on the ratio scale that difference is amplified. A PCV-specific correction is a genuinely open problem — `correct_detection_bias()`'s `pcv` output should be treated as diagnostic, not corrected. (Guarded by regression test so this documented limitation can't be silently forgotten.) Exception: under *strong outcome-dependent* detection the sign flips — naive PCV degrades badly (+15.3) and the corrected value is closer to truth (+3.1).
3. **Sparsity destroys the PCV entirely** (+44.6pp): with sparse strata, `var_main` truncates to zero and the PCV reads ~100% — a spurious "everything is additive, no intersectional effects" conclusion exactly in the settings where intersectional analysis is most wanted. No correction exists in this package (or, to our knowledge, anywhere); flagging this inflation is itself an actionable warning for applied work.

## 6. Joint capstone: detection bias and sparsity co-occurring (M2-3)

The realistic worst case: `sparse=True`, `interaction_sd=0.9`, `detection_strength=0.8`, `n=3500` (the package's Scenario D+E combined). Four arms on the *same* data, 20 replicates, true VPC = 20.15%:

| Arm | Bias | RMSE | Non-inferior vs. naive? |
|---|---:|---:|---|
| naive (no correction) | −10.30 | 12.28 | (reference) |
| detection-only (`correct_detection_bias`) | **−0.75** | 1.95 | Yes |
| sparse-only (`sparse_strata_vpc` on detection-biased data) | −3.74 | 9.21 | Yes |
| both, composed sequentially | **+8.79** | 11.67 | Yes (barely — within the 0.5pp margin of naive) |

**Detection correction alone removes almost all of the joint bias** (−10.30 → −0.75pp) — in this scenario the detection-bias mechanism dominates, and `correct_detection_bias()` handles it well even though the data are also sparse (a violation of neither method's core assumption, since detection correction doesn't assume dense strata). Sparse-strata correction alone, applied directly to detection-biased observed counts (a genuine assumption violation for `sparse_strata_vpc()`, whose calibration model assumes no detection process), still helps but leaves substantial residual bias.

**The composed correction (both applied sequentially) overcorrects past the truth to the opposite sign** (+8.79pp) and has worse RMSE than either individual correction. **The two corrections do not compose additively.** This is the study's most important finding for future methodological work: naively chaining a detection-bias correction into a sparse-strata calibration double-counts part of the adjustment, because the sparse-strata calibration curve is built assuming the *input* naive VPC's bias comes entirely from sampling noise — once that input has already been detection-corrected, the remaining "sparsity gap" the calibration tries to close is smaller than the calibration curve assumes, and it over-corrects. A jointly-calibrated correction that models both processes simultaneously (rather than composing two independently-calibrated corrections) is flagged as a concrete direction for future work (§7).

### 6.1 Generalization grid (previously single-scenario; now completed)

`scripts/validation/composition_grid.py` extends the capstone to a 3×2 grid — `detection_strength ∈ {0.4, 0.8, 1.2}` × `sparse ∈ {False, True}` (n=3500, n_rep=12), five arms per cell, all against the **population** estimand (this differs from the capstone table above, which used the finite-sample `y_true` refit; see `PHASE3_CAUSAL_IDENTIFICATION.md` §3.1 for why the two truths differ). Mean bias (pp):

| δ | sparse | naive | det-only | sparse-only | composed | **joint** |
|---:|---|---:|---:|---:|---:|---:|
| 0.4 | no | −13.81 | −9.32 | −7.27 | −2.13 | **−2.69** |
| 0.4 | yes | −13.75 | −7.55 | −7.16 | −0.38 | **−0.39** |
| 0.8 | no | −15.96 | −8.69 | −11.08 | −0.55 | **−2.65** |
| 0.8 | yes | −17.53 | −6.03 | −10.04 | +3.35 | **+1.49** |
| 1.2 | no | −12.52 | −7.09 | −6.39 | +7.55 | **−2.84** |
| 1.2 | yes | −14.12 | −4.47 | −9.97 | +7.27 | **+2.54** |

The single-scenario conclusion generalizes, with an important refinement: the composed estimator's bias **climbs monotonically with δ and changes sign** (−2.1 → +7.6), so its apparent accuracy at low δ is two opposite-signed errors cancelling, not reliability. The jointly-calibrated estimator is the only arm that stays stable (within ±3pp) in **every** cell, and ranks best-or-second-best by |bias| in 6/6 cells. Against the population estimand, detection-only correction is uniformly biased low (−4.5 to −9.3pp) because it inherits the sparse-shrinkage the joint calibration removes.

## 7. Limitations and deferred work

- **`baseline_logit` misspecification (§4.3)** was folded into Study 1 as a lightweight table rather than a full scenario: because detection is normalized relative to the `score=0` stratum, this parameter only rescales curvature and is empirically far less consequential than the score's covariates or functional form.
- **Number of strata (K)** — *completed since first writing*: Study 3 (§5.4) now covers K ∈ {8, 36, 108} for `sparse_strata_vpc()` via the generalized `simulate_k_strata()` generator. K-generalization for the *joint* estimator (`joint_calibrated_vpc`, which requires named covariates for its detection score) remains deferred.
- **PCV and `vpc_main`** — *completed since first writing*: Study 4 (§5.5) now covers the PCV directly. Its headline results are negative — the detection correction degrades the PCV even under correct specification, and sparsity inflates the naive PCV by ~45pp toward a spurious "purely additive" reading — so a PCV-specific *correction* (as opposed to this characterization) remains an open problem.
- **Sequential composition (§6)** — *completed since first writing*: the generalization grid (§6.1) now covers `detection_strength ∈ {0.4, 0.8, 1.2}` × both allocation regimes. The composed estimator's bias is direction-unstable in δ; the jointly-calibrated estimator (`joint_calibrated_vpc`, `PHASE3_CAUSAL_IDENTIFICATION.md` §3) is stable in every cell.
- All simulations use the package's existing 4-covariate (sex, education, wealth, rural) intersectional design; results should not be assumed to transfer to intersectional designs with different numbers or types of dimensions without re-running the corresponding scenarios.

## 8. References

1. Evans CR, Williams DR, Onnela J-P, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM - Population Health*. 2018;6:149–157. doi:10.1016/j.ssmph.2018.08.005
2. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. *Annals of Internal Medicine*. 2017;167(4):268–274. doi:10.7326/M16-2607
3. Elff M, Heisig JP, Schaeffer M, Shikano S. Multilevel analysis with few clusters: improving likelihood-based methods to provide unbiased estimates and accurate inference. *British Journal of Political Science*. 2021;51(1):412–426. doi:10.1017/S0007123419000097

---

## Appendix: reproducing this study

```bash
cd python && pip install -r requirements.txt
python ../scripts/validation/detection_bias_robustness.py     # Study 1 (~15s)
python ../scripts/validation/sparse_strata_robustness.py      # Study 2 (~16s)
python ../scripts/validation/joint_robustness_capstone.py     # Capstone (~5s)
pytest -q tests/test_detection_robustness.py tests/test_sparse_strata_robustness.py tests/test_joint_robustness.py
```

Each driver script writes its results table to a CSV alongside its figure in `figures/`, and prints the same pass/fail classification shown in this note. All numbers in §4–§6 were generated by these scripts, not transcribed by hand from other sources.

*No real data. Synthetic simulation only — this document makes no empirical claim about any population.*
