# test-detection-correction.R — Detection-bias sensitivity analysis
#
# Mirrors python/tests/test_detection_correction.py. The recovery property is
# the key claim: correcting observed data at the true under-detection strength
# moves the VPC back toward the value computed on the true (pre-detection)
# outcome, and the swept bounds bracket that truth.

vpc_on_true <- function(df) {
  dft <- df
  dft$y <- df$y_true
  fit_imaihda(dft)[["vpc_null"]]
}

test_that("delta = 0 reproduces the observed null-model VPC exactly", {
  df <- simulate_intersectional_data(n = 6000, interaction_sd = 0.9,
                                     detection_strength = 0.8, seed = 7)
  obs <- fit_imaihda(df)
  corr <- correct_detection_bias(df, delta = 0)
  expect_equal(corr[["vpc_null"]], obs[["vpc_null"]], tolerance = 1e-8)
  # Main-effects VPC uses a stratum-level WLS fit; matches the GLM closely.
  expect_lt(abs(corr[["vpc_main"]] - obs[["vpc_main"]]), 0.5)
})

test_that("correction at the true delta recovers the true VPC", {
  delta_true <- 0.8
  df <- simulate_intersectional_data(n = 12000, interaction_sd = 0.9,
                                     detection_strength = delta_true, seed = 7)
  vpc_true <- vpc_on_true(df)
  vpc_obs <- fit_imaihda(df)[["vpc_null"]]
  vpc_corr <- correct_detection_bias(df, delta = delta_true)[["vpc_null"]]

  expect_lt(vpc_obs, vpc_true)                              # detection masks VPC
  expect_lt(abs(vpc_corr - vpc_true), abs(vpc_obs - vpc_true))
  expect_lt(abs(vpc_corr - vpc_true), 2.0)                  # within ~2 pp
})

test_that("sensitivity bounds bracket the true VPC", {
  df <- simulate_intersectional_data(n = 12000, interaction_sd = 0.9,
                                     detection_strength = 0.8, seed = 7)
  vpc_true <- vpc_on_true(df)
  bounds <- vpc_detection_bounds(df, delta_max = 1.2, n_grid = 25)

  expect_equal(nrow(bounds), 25)
  expect_equal(bounds$delta[1], 0)
  expect_equal(bounds$vpc_null[1], fit_imaihda(df)[["vpc_null"]], tolerance = 1e-8)
  expect_true(min(bounds$vpc_null) <= vpc_true &&
                vpc_true <= max(bounds$vpc_null))
  expect_true(all(diff(bounds$implied_true_prevalence) >= -1e-9))
})

test_that("tipping point finds the crossing and returns NA when unreachable", {
  df <- simulate_intersectional_data(n = 12000, interaction_sd = 0.9,
                                     detection_strength = 0.8, seed = 7)
  vpc_obs <- fit_imaihda(df)[["vpc_null"]]
  target <- vpc_obs + 2.0
  tp <- detection_tipping_point(df, threshold = target, quantity = "vpc_null",
                                delta_max = 3.0)
  expect_false(is.na(tp))
  expect_gt(tp, 0)
  expect_lte(tp, 3.0)
  expect_lt(abs(correct_detection_bias(df, tp)[["vpc_null"]] - target), 1e-2)

  expect_true(is.na(detection_tipping_point(df, threshold = vpc_obs - 5.0,
                                            quantity = "vpc_null", delta_max = 1.0)))
})

test_that("score length is validated", {
  df <- simulate_intersectional_data(n = 2000, seed = 1)
  expect_error(correct_detection_bias(df, delta = 0.5, score = c(1, 2)),
               "score has length")
})
