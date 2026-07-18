# test-sparse-strata-ci.R — Sparse-strata bias correction and CI
#
# Mirrors python/tests/test_sparse_strata_ci.py. Ground truth is computed
# analytically from the generative formula in simulate.R (no simulation, no
# estimator bias), then compared against the naive vs. bias-corrected VPC
# under genuinely sparse stratum allocation (sparse = TRUE, matching the
# package's own Scenario E).

interaction_sd <- 0.15
interaction_seed <- 999L

analytic_truth <- function(interaction_sd = 0.15, seed = 999L) {
  strata <- stratum_table()
  eta_add <- -2.10 +
    0.20 * strata$sex + 0.35 * strata$education +
    0.30 * strata$wealth + 0.25 * strata$rural

  set.seed(seed + 999L)
  raw <- stats::rnorm(nrow(strata), mean = 0, sd = interaction_sd)
  idx_high <- which(strata$education == 2 & strata$wealth == 2 & strata$rural == 1)
  raw[idx_high] <- raw[idx_high] + 0.90
  idx_low <- which(strata$education == 0 & strata$wealth == 2 & strata$rural == 0)
  raw[idx_low] <- raw[idx_low] - 0.60

  eta_true <- eta_add + (raw - mean(raw))
  vpc_latent(stats::var(eta_true))
}

truth <- analytic_truth(interaction_sd, interaction_seed)

test_that("naive VPC underestimates under genuine sparsity", {
  naive <- numeric(6)
  for (seed in seq_len(6) - 1) {
    df <- simulate_intersectional_data(n = 3500, interaction_sd = interaction_sd,
                                       sparse = TRUE, seed = seed)
    naive[seed + 1] <- sparse_strata_vpc(df, seed = 100 + seed)$vpc_null_naive
  }
  expect_lt(mean(naive), truth - 2.0)
})

test_that("bias correction reduces average gap to truth", {
  naive <- numeric(6)
  corrected <- numeric(6)
  for (seed in seq_len(6) - 1) {
    df <- simulate_intersectional_data(n = 3500, interaction_sd = interaction_sd,
                                       sparse = TRUE, seed = seed)
    r <- sparse_strata_vpc(df, seed = 100 + seed)
    naive[seed + 1] <- r$vpc_null_naive
    corrected[seed + 1] <- r$vpc_null_corrected
  }
  expect_lt(abs(mean(corrected) - truth), abs(mean(naive) - truth))
})

test_that("CI is well-formed and sparse flag is set", {
  df <- simulate_intersectional_data(n = 3500, interaction_sd = interaction_sd,
                                     sparse = TRUE, seed = 1)
  r <- sparse_strata_vpc(df, seed = 7)
  expect_lte(r$ci_lower, r$ci_upper)
  expect_gte(r$ci_lower, 0)
  expect_lte(r$ci_upper, 100)
  expect_true(r$sparse)
  expect_lt(r$min_stratum_n, 20)
})

test_that("dense well-powered data needs little correction", {
  df <- simulate_intersectional_data(n = 60000, interaction_sd = interaction_sd, seed = 3)
  r <- sparse_strata_vpc(df, seed = 42)
  expect_false(r$sparse)
  expect_lt(abs(r$vpc_null_corrected - r$vpc_null_naive), 1.5)
})

test_that("results are reproducible with the same seed", {
  df <- simulate_intersectional_data(n = 1500, interaction_sd = interaction_sd,
                                     sparse = TRUE, seed = 5)
  r1 <- sparse_strata_vpc(df, seed = 11)
  r2 <- sparse_strata_vpc(df, seed = 11)
  expect_equal(r1, r2)
})

test_that("level is validated", {
  df <- simulate_intersectional_data(n = 2000, seed = 1)
  expect_error(sparse_strata_vpc(df, level = 1.5), "level must be in")
  expect_error(sparse_strata_vpc(df, level = 0), "level must be in")
})
