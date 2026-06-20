# test-fit.R — Unit tests for fit_imaihda() and simulation

test_that("fit_imaihda validates input", {
  expect_error(fit_imaihda(data.frame(x = 1)), "missing required columns")
  expect_error(
    fit_imaihda(data.frame(stratum = factor(1), sex = 0, education = 0,
                           wealth = 0, rural = 0, y = 0)),
    "at least 2 strata"
  )
  expect_error(
    fit_imaihda(data.frame(stratum = factor(c(1, 2)), sex = c(0, 1),
                           education = c(0, 0), wealth = c(0, 0),
                           rural = c(0, 0), y = c(0.5, 1))),
    "binary"
  )
})

test_that("fit_imaihda returns expected structure", {
  df <- simulate_intersectional_data(n = 1000, seed = 42)
  res <- fit_imaihda(df)

  expect_type(res, "list")
  expect_true("var_null" %in% names(res))
  expect_true("vpc_null" %in% names(res))
  expect_true("var_main" %in% names(res))
  expect_true("vpc_main" %in% names(res))
  expect_true("pcv" %in% names(res))
  expect_true("n_strata" %in% names(res))
  expect_true("overall_prevalence" %in% names(res))
  expect_equal(res[["estimator"]], "fast_empirical_logit_diagnostic")
})

test_that("additive scenario is additive-dominant", {
  df <- simulate_intersectional_data(n = 1800, seed = 123)
  res <- fit_imaihda(df)
  expect_true(res$pcv > 70)
  expect_true(res$vpc_main < res$vpc_null)
})

test_that("interaction scenario increases VPC", {
  df_add <- simulate_intersectional_data(n = 1500, seed = 42)
  df_int <- simulate_intersectional_data(n = 1500, interaction_sd = 0.9,
                                         seed = 42)
  res_add <- fit_imaihda(df_add)
  res_int <- fit_imaihda(df_int)
  expect_true(res_int$vpc_null > res_add$vpc_null)
})

test_that("detection bias reduces observed prevalence", {
  d0 <- simulate_intersectional_data(n = 1200, seed = 123,
                                     detection_strength = 0)
  d1 <- simulate_intersectional_data(n = 1200, seed = 123,
                                     detection_strength = 1.0)
  expect_true(mean(d1$y) < mean(d0$y))
  expect_equal(mean(d1$y_true), mean(d0$y_true))
})

test_that("simulation produces expected structure", {
  df <- simulate_intersectional_data(n = 1000, seed = 123)

  expect_equal(nrow(df), 1000)
  expect_equal(length(unique(df$stratum)), 36)
  expect_true(all(df$y %in% c(0, 1)))
  expect_true(mean(df$y) > 0 && mean(df$y) < 1)
  expected_cols <- c("sex", "education", "wealth", "rural", "stratum", "y")
  expect_true(all(expected_cols %in% names(df)))
  
  # With detection and interaction: expect extra columns
  df2 <- simulate_intersectional_data(n = 500, detection_strength = 0.8,
                                       interaction_sd = 0.9, seed = 42)
  expect_true("y_true" %in% names(df2))
  expect_true("detection_probability" %in% names(df2))
  expect_true("true_stratum_residual" %in% names(df2))
})

test_that("sparse simulation creates uneven stratum sizes", {
  df <- simulate_intersectional_data(n = 400, sparse = TRUE, seed = 1)
  sizes <- table(df$stratum)

  # Stratum sizes should vary substantially (gamma allocation)
  ratio <- as.numeric(max(sizes)) / as.numeric(min(sizes))
  expect_true(ratio > 10,
    info = sprintf("Expected max/min > 10, got %.1f", ratio))

  # Mean stratum size should be close to n/36 but variance should be high
  cv <- sd(as.numeric(sizes)) / mean(as.numeric(sizes))
  expect_true(cv > 0.5,
    info = sprintf("Expected CV > 0.5, got %.2f", cv))
})

test_that("all benchmark scenarios pass", {
  grid <- scenario_grid()
  rows <- lapply(names(grid), function(nm) fit_scenario(nm, grid[[nm]]))
  results <- do.call(rbind, lapply(rows, as.data.frame))
  benchmarks <- evaluate_benchmarks(results)
  expect_true(all(as.logical(benchmarks$passed)))
})
