# test_simulation.R — Unit tests for I-MAIHDA R simulation
#
# Replicates tests/test_simulation.py
# Run with: Rscript -e "testthat::test_file('R/test_simulation.R')"
#   (from the imaihda_hic_mic_v31/ directory)

library(testthat)

# Source R modules (path relative to project root, find it)
find_root <- function() {
  # When run via Rscript, the --file= argument gives the script path
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(dirname(script_path)), winslash = "/"))
  }
  # Interactive: try relative paths
  if (file.exists("R/simulate.R")) return(normalizePath("."))
  if (file.exists("../R/simulate.R")) return(normalizePath(".."))
  # Fallback: try relative to this file's location
  tryCatch({
    test_dir <- dirname(sys.frame(1)$ofile)
    if (file.exists(file.path(test_dir, "../R/simulate.R")))
      return(normalizePath(file.path(test_dir, "..")))
  }, error = function(e) NULL)
  stop("Cannot find R/ source files")
}

ROOT <- find_root()
source(file.path(ROOT, "R", "simulate.R"))
source(file.path(ROOT, "R", "fit.R"))

test_that("vpc_pcv_formulas", {
  expect_equal(vpc_latent(0), 0, tolerance = 1e-8)
  expect_true(vpc_latent(0.5) > 0 && vpc_latent(0.5) < 100)
  expect_equal(pcv(1.0, 0.25), 75.0)
})

test_that("simulation_has_expected_structure", {
  df <- simulate_intersectional_data(n = 1000, seed = 123)
  expect_equal(length(unique(df$stratum)), 36)
  expect_true(all(df$y %in% c(0, 1)))
  expect_true(mean(df$y) > 0 && mean(df$y) < 1)
})

test_that("additive_scenario_fits_and_is_additive_dominant", {
  df <- simulate_intersectional_data(n = 1800, seed = 123)
  res <- fit_imaihda(df)
  expect_true(res$pcv > 70)
  expect_true(res$vpc_main < res$vpc_null)
})

test_that("detection_bias_reduces_observed_prevalence", {
  d0 <- simulate_intersectional_data(n = 1200, seed = 123, detection_strength = 0)
  d1 <- simulate_intersectional_data(n = 1200, seed = 123, detection_strength = 1.0)
  expect_true(mean(d1$y) < mean(d0$y))
  expect_equal(mean(d1$y_true), mean(d0$y_true))
})
