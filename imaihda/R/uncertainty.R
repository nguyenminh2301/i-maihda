# uncertainty.R — Small-sample bias correction and confidence intervals for
# sparse intersectional strata.
#
# See python/imaihda_sim/uncertainty.py for the full derivation. This module
# mirrors that Python implementation exactly. Its "naive" estimator uses the
# same unweighted variance formula as between_stratum_variance()
# (diagnostics.R), which is also what fit_imaihda()'s R implementation uses,
# so vpc_null_naive here matches fit_imaihda(df)$vpc_null on the same data.

# A stratum smaller than this is where Laplace-smoothing shrinkage is largest
# and the naive fast estimator's downward bias becomes material.
.sparse_threshold <- 20

.var_null_from_counts <- function(events, n, p_overall) {
  p_s <- (events + 0.5) / (n + 1.0)
  logit_s <- log(p_s / (1.0 - p_s))
  sampling_var <- 1.0 / (events + 0.5) + 1.0 / (n - events + 0.5)
  null_pred <- rep(log(p_overall / (1.0 - p_overall)), length(events))
  between_stratum_variance(logit_s, null_pred, sampling_var)
}

.simulate_naive_var <- function(sigma2, n_j, p_overall, n_rep) {
  sd <- sqrt(max(sigma2, 0))
  k <- length(n_j)
  out <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    effects <- if (sd > 0) stats::rnorm(k, mean = 0, sd = sd) else rep(0, k)
    p_j <- pmin(pmax(logit_inv(logit(p_overall) + effects), 1e-6), 1 - 1e-6)
    events <- stats::rbinom(k, size = n_j, prob = p_j)
    out[i] <- .var_null_from_counts(events, n_j, p_overall)
  }
  out
}

.build_calibration_curve <- function(n_j, p_overall, sigma2_max, n_grid, n_sim) {
  grid <- seq(0, sigma2_max, length.out = n_grid)
  means <- vapply(
    grid,
    function(s) mean(.simulate_naive_var(s, n_j, p_overall, n_sim)),
    numeric(1)
  )
  list(grid = grid, means = means)
}

#' Bias-corrected null-model VPC and confidence interval for sparse strata
#'
#' The fast null-model VPC estimator computes between-stratum variance from
#' Laplace-smoothed empirical logits \code{(events + 0.5) / (n + 1)}. When a
#' stratum has few individuals, that smoothing pulls its logit toward the
#' population mean, which shrinks the observed between-stratum spread below
#' its true value. In sparse-strata regimes the fast VPC is therefore
#' systematically underestimated — worst exactly where analysts most need a
#' warning. No confidence interval or small-sample correction for this
#' estimator has been reported in the I-MAIHDA literature.
#'
#' This function corrects it by calibration: for the observed stratum sizes,
#' it simulates the null-model estimator's own expected value under a grid of
#' candidate true between-stratum variances (holding stratum sizes fixed), and
#' inverts that curve to find the true variance whose expected naive estimate
#' matches what was actually observed (an indirect-inference / SIMEX-style
#' correction). Bootstrap replicates of the naive estimator at the calibrated
#' variance are pushed through the same inverse mapping to build a confidence
#' interval for the bias-corrected VPC.
#'
#' @param df A data.frame with at least columns \code{stratum} and \code{y}
#'   (binary outcome 0/1), e.g. from \code{\link{simulate_intersectional_data}}.
#' @param sigma2_max Numeric. Upper bound of the calibration grid on the
#'   between-stratum variance (logit scale). Default 2.0.
#' @param n_grid Integer. Number of grid points spanning
#'   \code{[0, sigma2_max]} for the calibration curve. Default 25.
#' @param n_sim Integer. Monte Carlo replicates per grid point used to
#'   estimate the expected naive estimator. Default 200.
#' @param n_boot Integer. Bootstrap replicates (at the calibrated variance)
#'   used to build the confidence interval. Default 300.
#' @param level Numeric in (0, 1). Confidence level. Default 0.95.
#' @param seed Integer, optional. RNG seed for reproducibility.
#'
#' @return A named list: vpc_null_naive, vpc_null_corrected, var_null_naive,
#'   var_null_corrected, ci_lower, ci_upper, level, n_strata, min_stratum_n,
#'   sparse, n_sim, n_boot.
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 3500, interaction_sd = 0.15,
#'                                    sparse = TRUE, seed = 1)
#' sparse_strata_vpc(df, seed = 7)
#' }
#'
#' @export
sparse_strata_vpc <- function(df, sigma2_max = 2.0, n_grid = 25L, n_sim = 200L,
                              n_boot = 300L, level = 0.95, seed = NULL) {
  if (!(level > 0 && level < 1)) stop("level must be in (0, 1)")
  if (!is.null(seed)) set.seed(seed)

  agg <- stats::aggregate(
    df$y,
    by = list(stratum = df$stratum),
    FUN = function(x) c(events = sum(x), n = length(x))
  )
  n_j <- as.numeric(agg$x[, "n"])
  events_j <- as.numeric(agg$x[, "events"])
  if (length(n_j) < 2) stop("Need at least 2 strata; found ", length(n_j))

  p_overall <- min(max(sum(events_j) / sum(n_j), 1e-6), 1 - 1e-6)
  var_naive <- .var_null_from_counts(events_j, n_j, p_overall)

  curve <- .build_calibration_curve(n_j, p_overall, sigma2_max, n_grid, n_sim)
  grid <- curve$grid
  means <- curve$means

  if (var_naive <= means[1]) {
    var_corrected <- 0
  } else if (var_naive >= means[length(means)]) {
    var_corrected <- grid[length(grid)]
  } else {
    var_corrected <- stats::approx(means, grid, xout = var_naive)$y
  }

  boot_naive <- .simulate_naive_var(var_corrected, n_j, p_overall, n_boot)
  boot_naive_clipped <- pmin(pmax(boot_naive, means[1]), means[length(means)])
  boot_corrected <- stats::approx(means, grid, xout = boot_naive_clipped)$y
  alpha <- 1 - level
  ci <- stats::quantile(boot_corrected, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

  list(
    vpc_null_naive     = vpc_latent(var_naive),
    vpc_null_corrected = vpc_latent(var_corrected),
    var_null_naive     = var_naive,
    var_null_corrected = var_corrected,
    ci_lower           = vpc_latent(ci[1]),
    ci_upper           = vpc_latent(ci[2]),
    level              = level,
    n_strata           = length(n_j),
    min_stratum_n      = min(n_j),
    sparse             = min(n_j) < .sparse_threshold,
    n_sim              = n_sim,
    n_boot             = n_boot
  )
}
