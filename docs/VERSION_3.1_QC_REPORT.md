# v3.1 QC report

## Document corrections

- Removed the over-specific phrase "beyond high-income Europe" from the main title, because the empirical design includes HRS as a high-income non-European comparator.
- Reframed the project as a structured HIC-MIC comparison, not a formal causal test of welfare-regime effects.
- Removed the REML + Student-t sentence from the exposé because binary I-MAIHDA is a logistic mixed-model setting and REML is not directly applicable to GLMMs.
- Replaced "inflates middle-income inequality" with "attenuates, exaggerates, or redistributes" because differential detection can bias observed gradients in either direction.
- Moved biological aging to a conditional mechanistic extension after a biomarker harmonization audit.
- Removed overclaiming about causal mediation unless defensible temporal ordering is available.
- Corrected the compression-of-morbidity reference to *Archives of Gerontology and Geriatrics*.

## Code corrections

- Original R workflow could not be executed in this runtime because R/Rscript is unavailable.
- Built a tested Python workflow using synthetic data only.
- Added benchmark checks, detection sweep, figures, metadata, and tests.
- The code demonstrates the narrow issue it claims to demonstrate: VPC/PCV are sensitive to prevalence, stratum sparsity, and differential detection, so raw HIC-MIC comparisons require diagnostics.

## Remaining methodological caution

The repository is a competence demo, not an empirical estimator for the PhD. It should be used to support interview discussion, not as a published method claim.
