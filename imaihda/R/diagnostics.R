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
#' variance of stratum-level residuals. Uses unweighted (sample) variance
#' rather than precision-weighted variance, because precision weights
#' downweight extreme strata that carry the most between-stratum signal.
#'
#' The estimator:
#' \deqn{\hat{\sigma}^2_{stratum} = \max(0,\; \mathrm{Var}(\mathrm{residual})
#'       - \overline{\mathrm{Var}(\mathrm{logit}_j)})}
#'
#' where \eqn{\mathrm{Var}} is the standard (unweighted) sample variance and
#' \eqn{\overline{\mathrm{Var}(\mathrm{logit}_j)}} is the arithmetic mean of
#' stratum-level delta-method sampling variances.
#'
#' @param logit_vec Numeric vector. Empirical logits for each stratum.
#' @param pred_vec Numeric vector. Model-predicted logits for each stratum.
#' @param sampling_var Numeric vector. Delta-method sampling variance
#'   of each empirical logit.
#' @param weights Numeric vector. Kept for backward compatibility; not used.
#'
#' @return Non-negative numeric scalar: estimated between-stratum variance.
#'
#' @keywords internal
between_stratum_variance <- function(logit_vec, pred_vec, sampling_var, weights = NULL) {
  residual     <- logit_vec - pred_vec
  raw_var      <- stats::var(residual)            # unweighted, (n-1) divisor
  expected_noise <- mean(sampling_var)             # arithmetic mean
  max(raw_var - expected_noise, 0)
}
