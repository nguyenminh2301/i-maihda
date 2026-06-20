# diagnostics.R — Method-of-moments estimator for between-stratum variance
#
# Internal functions implementing the fast empirical-logit diagnostic
# for MAIHDA. Not exported.

#' Weighted variance
#'
#' Computes the weighted variance of a vector:
#' \eqn{\sum w_i (x_i - \bar{x}_w)^2 / \sum w_i}
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of non-negative weights (same length as x).
#'
#' @return Weighted variance (scalar).
#'
#' @keywords internal
weighted_variance <- function(x, w) {
  mu <- weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

#' Between-stratum variance estimator (method of moments)
#'
#' Estimates the true between-stratum variance on the logit scale by
#' subtracting the expected binomial sampling noise from the observed
#' weighted variance of stratum-level residuals. This is a fast
#' synthetic-data diagnostic, not a replacement for full GLMM estimation.
#'
#' The estimator:
#' \deqn{\hat{\sigma}^2_{stratum} = \max(0,\; \mathrm{Var}_w(\mathrm{residual})
#'       - \mathbb{E}_w[\mathrm{Var}(\mathrm{logit}_j)])}
#'
#' where \eqn{\mathrm{Var}_w} is the precision-weighted variance and
#' \eqn{\mathbb{E}_w} is the precision-weighted mean of stratum-level
#' delta-method sampling variances.
#'
#' @param logit_vec Numeric vector. Empirical logits for each stratum.
#' @param pred_vec Numeric vector. Model-predicted logits for each stratum.
#' @param sampling_var Numeric vector. Delta-method sampling variance
#'   of each empirical logit.
#' @param weights Numeric vector. Precision weights (1 / sampling_var).
#'
#' @return Non-negative numeric scalar: estimated between-stratum variance.
#'
#' @keywords internal
between_stratum_variance <- function(logit_vec, pred_vec, sampling_var, weights) {
  residual     <- logit_vec - pred_vec
  raw_var      <- weighted_variance(residual, weights)
  expected_noise <- weighted.mean(sampling_var, weights)
  max(raw_var - expected_noise, 0)
}
