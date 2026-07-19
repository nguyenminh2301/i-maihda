# Changelog

All notable changes to this repository. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versions refer to
the R package `imaihda` (`imaihda/DESCRIPTION`); the Python implementation
(`python/imaihda_sim/`) evolves in lockstep. Release tags will be created
once the corresponding work is merged to `main`.

## [Unreleased]

### Added
- **Phase 4 — cross-cohort comparability decomposition** (`docs/PHASE4_CROSS_COHORT_DECOMPOSITION.md`,
  `python/imaihda_sim/cohort.py`): `cross_cohort_decomposition()` splits a raw
  HIC–MIC VPC/PCV gap into a Shapley contribution for each of the four artefact
  channels named in the PhD exposé (prevalence, sparse strata, selective
  attrition, differential detection) plus a structural residual, with a CI on
  the structural share and a distinguishable-from-zero verdict. Self-validated
  against a known-by-construction true gap: null-structure → 0.00pp structural
  share; masked structure → correctly unmasked. Adds the fourth artefact
  generator `simulate_selective_attrition()`. `docs/PROJECT_PROGRESS_REPORT.md`
  maps the whole methodological programme onto the PhD subprojects.
- **Robustness/simulation study** for both quantitative-bias-analysis methods
  under violated assumptions (`docs/METHODS_NOTE_ROBUSTNESS.md`;
  `python/imaihda_sim/robustness.py`; `scripts/validation/*_robustness.py`,
  `joint_robustness_capstone.py`): score misspecification, outcome-dependent
  detection (identification breakdown at `severity_weight ≈ 0.8–0.9`),
  non-Gaussian random effects, rare prevalence × sparsity, and the joint
  detection+sparsity capstone showing the two corrections do not compose
  additively. Formal non-inferiority / bias-reduction / coverage success
  criteria (§3).
- **Phase 3 — causal identification** (`docs/PHASE3_CAUSAL_IDENTIFICATION.md`;
  `python/imaihda_sim/causal.py`): identification propositions (point-ID
  under covariate-only detection; provable non-ID under outcome-dependent
  detection); `joint_calibrated_vpc()` (correlation-aware joint calibration
  targeting the population estimand; capstone bias −16.3 → +0.9pp);
  `vpc_partial_bounds()` (sharp, globally-verified partial-identification
  bounds valid where point correction fails; truth-in-bounds 45/45).
- **K-generalization study** (`simulate_k_strata()`;
  `scripts/validation/k_generalization_study.py`): K ∈ {8, 36, 108};
  non-inferiority holds at every K; exposed a truncation-floor CI defect at
  K=108 (bootstrap coverage 0.10), fixed by the new
  `sparse_strata_vpc(ci_method="test_inversion")` option (coverage 0.95)
  with an `at_floor` output flag. Default CI method unchanged.
- **PCV robustness pass** (`scripts/validation/pcv_robustness.py`): naive
  PCV is robust to covariate-only detection (+3.7pp) but the detection
  correction *degrades* PCV even under correct specification (+15.5pp) —
  treat `correct_detection_bias()`'s `pcv` as diagnostic; sparsity inflates
  naive PCV by ~45pp toward a spurious "purely additive" reading. Guarded by
  regression tests.
- **Composition grid** (`scripts/validation/composition_grid.py`): 3×2
  detection×allocation grid; sequential composition is direction-unstable in
  δ; `joint_calibrated_vpc()` stable within ±3pp in all 6 cells.
- `fit_imaihda(weighting="unweighted")` (Python): matches R's
  v0.2.1-corrected variance formula exactly; default `"precision"` unchanged
  for benchmark continuity.
- GitHub Actions CI running the Python test suite (46 tests).

## [0.4.0] — 2026-07

### Added
- `sparse_strata_vpc()` (Python `uncertainty.py` + R `uncertainty.R`):
  bias-corrected null-model VPC with confidence interval for sparse strata,
  via simulation calibration (SIMEX-style curve inversion). Validated against
  an analytic ground truth: mean bias 4.2pp → 0.3pp under Scenario-E
  sparsity, 10/12 nominal-95% CI coverage.

### Fixed / documented
- Documented that Python's `fit_imaihda()` precision-weighted variance
  formula carries a small persistent bias absent from R's v0.2.1 unweighted
  formula (closed in Unreleased via the `weighting=` option).

## [0.3.0] — 2026-07

### Added
- Detection-bias sensitivity analysis (Python `detection.py` + R
  `detection.R`): `correct_detection_bias()`, `vpc_detection_bounds()`,
  `detection_tipping_point()` (E-value analogue), `plot_detection_bounds()`.
  Self-validating recovery test: correcting at the true δ recovers the true
  VPC within ~0.5pp; swept bounds bracket the truth.

## [0.2.2] — 2026-06

### Fixed
- Synced README version references with `DESCRIPTION` (0.2.1 → 0.2.2).

## [0.2.x and earlier] — 2026-06

- Initial public history (single development burst, 2026-06-20): I-MAIHDA
  HIC-MIC simulation v3.1 → R package restructuring → `method="glmer"`
  estimator and CRAN MAIHDA cross-validation (v0.2.0) → corrected unweighted
  fast estimator within <1pp of GLMM (v0.2.1) → project compaction (v0.2.2).
  See `git log` and `docs/VERSION_3.1_QC_REPORT.md`.
