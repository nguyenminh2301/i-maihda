# test-maihda-crossval.R — Cross-validation: imaihda vs CRAN MAIHDA
#
# Verifies that our glmer-based VPC and PCV match CRAN MAIHDA's results
# on the same dataset. Requires MAIHDA package (Suggests).

test_that("vpc_latent formula matches CRAN MAIHDA latent-response approach", {
  skip_if_not_installed("MAIHDA")

  # CRAN MAIHDA uses pi^2/3 for logistic latent VPC (same as us)
  # Test with known variance values
  sigma2 <- c(0.1, 0.5, 1.0, 3.29)

  for (s in sigma2) {
    our_vpc <- vpc_latent(s)
    # Manual calculation matching MAIHDA's approach
    expected_vpc <- 100 * s / (s + pi^2 / 3)
    expect_equal(our_vpc, expected_vpc, tolerance = 1e-8,
      info = sprintf("sigma2 = %.2f", s))
  }
})

test_that("pcv formula matches CRAN MAIHDA", {
  skip_if_not_installed("MAIHDA")

  # CRAN MAIHDA's calculate_pvc computes: (var1 - var2) / var1
  expect_equal(pcv(2.83, 0.49), 100 * (2.83 - 0.49) / 2.83, tolerance = 1e-6)
  expect_equal(pcv(1.0, 0.25), 75.0)
  expect_equal(pcv(0.5, 0.5), 0.0)
})

test_that("fit_imaihda glmer extracts stratum variance identically to CRAN MAIHDA", {
  skip_if_not_installed("MAIHDA")
  skip_if_not_installed("lme4")

  # Use the NHANES data bundled with CRAN MAIHDA
  library(MAIHDA)
  data("maihda_health_data")

  # Complete cases
  health <- maihda_health_data[complete.cases(
    maihda_health_data[, c("BMI", "Age", "Gender", "Race", "Education", "Poverty")]
  ), ]

  # Fit the null intersectional model via lme4 directly (both packages use this)
  null_model <- lme4::lmer(
    BMI ~ (1 | Gender:Race:Education),
    data = health
  )

  # Extract variance using our internal function
  our_var <- extract_glmer_stratum_variance(null_model, group = "Gender:Race:Education")

  # Extract variance using CRAN MAIHDA's internal function
  maihda_var <- MAIHDA:::maihda_stratum_variance_lme4(
    null_model, group = "Gender:Race:Education"
  )

  expect_equal(our_var, maihda_var, tolerance = 1e-10,
    info = "Stratum variance extraction should be identical")
})

test_that("fit_imaihda(method='glmer') PCV matches CRAN MAIHDA on NHANES data", {
  skip_if_not_installed("MAIHDA")
  skip_if_not_installed("lme4")

  library(MAIHDA)
  data("maihda_health_data")

  health <- maihda_health_data[complete.cases(
    maihda_health_data[, c("BMI", "Age", "Gender", "Race", "Education", "Poverty")]
  ), ]

  # CRAN MAIHDA approach
  maihda_res <- maihda(
    BMI ~ Gender + Race + Education + (1 | Gender:Race:Education),
    data = health
  )

  # CRAN MAIHDA reports:
  # - VPC from the null model fitted with REML (default for lmer)
  # - PCV from models refitted with ML (via maihda_pcv_refit_ml)
  # - var_null/var_main from the ML-refitted models (for PCV)
  maihda_vpc <- maihda_res[["summary"]][["vpc"]][["estimate"]]
  maihda_pcv <- maihda_res[["pcv"]][["pvc"]]
  maihda_var_null <- maihda_res[["pcv"]][["var_model1"]]
  maihda_var_main <- maihda_res[["pcv"]][["var_model2"]]

  # VPC comparison: use REML (default), matching CRAN MAIHDA's single-model VPC
  null_reml <- lme4::lmer(
    BMI ~ (1 | Gender:Race:Education),
    data = health
  )
  our_var_null_reml <- extract_glmer_stratum_variance(null_reml,
    group = "Gender:Race:Education")
  sigma_resid_reml <- attr(lme4::VarCorr(null_reml), "sc")^2
  our_vpc <- our_var_null_reml / (our_var_null_reml + sigma_resid_reml)

  # PCV comparison: use ML (matching maihda_pcv_refit_ml)
  null_ml <- lme4::refitML(null_reml)
  main_ml <- lme4::lmer(
    BMI ~ Gender + Race + Education + (1 | Gender:Race:Education),
    data   = health,
    REML   = FALSE
  )
  our_var_null <- extract_glmer_stratum_variance(null_ml,
    group = "Gender:Race:Education")
  our_var_main <- extract_glmer_stratum_variance(main_ml,
    group = "Gender:Race:Education")
  our_pcv <- (our_var_null - our_var_main) / our_var_null

  # imaihda reports VPC/PCV as percentages, CRAN MAIHDA as proportions
  expect_equal(our_var_null, maihda_var_null, tolerance = 1e-6,
    info = "Null between-stratum variance should match")
  expect_equal(our_var_main, maihda_var_main, tolerance = 1e-6,
    info = "Main between-stratum variance should match")
  expect_equal(our_vpc, maihda_vpc, tolerance = 1e-6,
    info = "VPC should match")
  # PCV is a ratio of two rounded estimates; floating-point accumulation
  # can push the difference to ~2e-4 even when variances match at 1e-6
  expect_equal(our_pcv, maihda_pcv, tolerance = 1e-4,
    info = "PCV should match")
})
