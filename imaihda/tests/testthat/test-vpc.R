# test-vpc.R — Unit tests for vpc_latent() and pcv()

test_that("vpc_latent returns correct values", {
  expect_equal(vpc_latent(0), 0, tolerance = 1e-8)
  expect_true(vpc_latent(0.5) > 0 && vpc_latent(0.5) < 100)
  # When stratum_variance equals pi^2/3, VPC should be 50%
  expect_equal(vpc_latent(pi^2 / 3), 50, tolerance = 1e-6)
})

test_that("vpc_latent rejects invalid inputs", {
  expect_error(vpc_latent(-0.1), "must be non-negative")
  expect_error(vpc_latent("a"), "must be a single numeric value")
  expect_error(vpc_latent(c(0.1, 0.2)), "must be a single numeric value")
})

test_that("pcv returns correct values", {
  expect_equal(pcv(1.0, 0.25), 75.0)
  expect_equal(pcv(2.0, 0.5), 75.0)
  expect_equal(pcv(0.5, 0.0), 100.0)
  expect_equal(pcv(0.5, 0.5), 0.0)
})

test_that("pcv returns NaN when null variance is zero or negative", {
  expect_true(is.nan(pcv(0, 0)))
  expect_true(is.nan(pcv(-0.1, 0)))
})

test_that("pcv rejects invalid inputs", {
  expect_error(pcv("a", 1), "must be a single numeric value")
  expect_error(pcv(1, c(1, 2)), "must be a single numeric value")
})
