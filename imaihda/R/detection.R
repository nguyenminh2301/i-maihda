# detection.R — Detection-bias sensitivity analysis for intersectional MAIHDA
#
# Quantitative bias analysis (QBA) for the between-stratum VPC and PCV under
# SES-patterned under-detection of the outcome. Unlike CRAN MAIHDA — which
# assumes the outcome is measured equally well across strata — these tools ask
# how far the observed VPC/PCV could be from the truth if disadvantaged strata
# under-record true cases. No equivalent exists in CRAN MAIHDA.

# Per-stratum disadvantage score, matching the simulator's detection model
# (simulate.R): logit(P(detect)) = baseline - delta * (edu + wealth + 0.4*rural).
.rural_weight <- 0.40

.stratum_score <- function(strata, score = NULL) {
  if (!is.null(score)) {
    if (length(score) != nrow(strata)) {
      stop("score has length ", length(score),
           " but there are ", nrow(strata), " strata")
    }
    return(as.numeric(score))
  }
  as.numeric(strata$education) +
    as.numeric(strata$wealth) +
    .rural_weight * as.numeric(strata$rural)
}

#' Correct VPC/PCV for hypothesized SES-patterned under-detection
#'
#' Quantitative bias analysis for intersectional MAIHDA. The observed outcome is
#' \code{y = y_true * detected}; detection removes only true positives, so within
#' a stratum the observed prevalence deflates by the stratum detection
#' probability: \eqn{E[p_{obs}] = p_{true} \times d(\delta)}, hence
#' \eqn{p_{true} = p_{obs} / d(\delta)}. \code{delta} is a sensitivity parameter;
#' \code{delta = 0} applies no correction and returns the observed VPC/PCV.
#'
#' Detection is taken \emph{relative} to the most-advantaged stratum
#' (\code{score = 0}), because only the SES-patterned \emph{differential}
#' under-detection is identifiable from observed data; a uniform ascertainment
#' level that hits every stratum equally is not.
#'
#' Operates entirely at the stratum level, so it is fast and needs no
#' individual-level refit. The main-effects model is a weighted least-squares fit
#' of the corrected stratum logits (the stratum-level analogue of the additive
#' GLM in \code{\link{fit_imaihda}}).
#'
#' @param df A data.frame with columns stratum, sex, education, wealth, rural, y.
#' @param delta Numeric >= 0. Strength of SES-patterned under-detection.
#' @param baseline_logit Numeric. Controls how steeply detection declines with
#'   disadvantage (default 2.0, matching the simulator). \code{delta = 0} is the
#'   identity regardless of this value.
#' @param score Optional numeric vector (one value per stratum) overriding the
#'   default disadvantage score \code{education + wealth + 0.4 * rural}.
#'
#' @return A named list: delta, baseline_logit, var_null, vpc_null, var_main,
#'   vpc_main, pcv, observed_prevalence, implied_true_prevalence,
#'   min_relative_detection.
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 8000, interaction_sd = 0.9,
#'                                    detection_strength = 0.8, seed = 7)
#' correct_detection_bias(df, delta = 0.0)   # observed
#' correct_detection_bias(df, delta = 0.8)   # corrected
#' }
#'
#' @export
correct_detection_bias <- function(df, delta, baseline_logit = 2.0, score = NULL) {
  if (delta < 0) stop("delta must be non-negative")

  strata <- aggregate_strata(df)
  e_obs <- as.numeric(strata$events)
  n <- as.numeric(strata$n)

  s <- .stratum_score(strata, score)
  d0 <- logit_inv(baseline_logit)
  d <- logit_inv(baseline_logit - delta * s) / d0
  d <- pmin(pmax(d, 1e-6), 1.0)

  # Invert the detection process: reconstruct implied true events, capped at n.
  e_true <- pmin(n, e_obs / d)

  # Recompute Laplace-smoothed logits and delta-method sampling variance from
  # the corrected counts (same formulas as aggregate_strata()).
  p_corr <- (e_true + 0.5) / (n + 1.0)
  logit_corr <- log(p_corr / (1.0 - p_corr))
  sampling_var <- 1.0 / (e_true + 0.5) + 1.0 / (n - e_true + 0.5)
  weights <- 1.0 / sampling_var

  p_overall <- min(max(sum(e_true) / sum(n), 1e-6), 1 - 1e-6)
  null_pred <- rep(log(p_overall / (1 - p_overall)), nrow(strata))
  var0 <- between_stratum_variance(logit_corr, null_pred, sampling_var)

  # Weighted additive fit of corrected stratum logits (main-effects model).
  fit_df <- data.frame(
    logit_corr = logit_corr,
    sex        = factor(strata$sex),
    education  = factor(strata$education),
    wealth     = factor(strata$wealth),
    rural      = factor(strata$rural)
  )
  main <- stats::lm(logit_corr ~ sex + education + wealth + rural,
                    data = fit_df, weights = weights)
  main_pred <- as.numeric(stats::fitted(main))
  var1 <- between_stratum_variance(logit_corr, main_pred, sampling_var)

  list(
    delta                   = as.numeric(delta),
    baseline_logit          = as.numeric(baseline_logit),
    var_null                = var0,
    vpc_null                = vpc_latent(var0),
    var_main                = var1,
    vpc_main                = vpc_latent(var1),
    pcv                     = pcv(var0, var1),
    observed_prevalence     = sum(e_obs) / sum(n),
    implied_true_prevalence = p_overall,
    min_relative_detection  = min(d)
  )
}

#' Sensitivity bounds on VPC/PCV over a range of under-detection strengths
#'
#' Sweeps \code{delta} from 0 to \code{delta_max} and tabulates the corrected
#' VPC/PCV. The range of \code{vpc_null} across rows is the sensitivity bound on
#' the true VPC.
#'
#' @param df A data.frame as for \code{\link{correct_detection_bias}}.
#' @param delta_max Numeric. Largest under-detection strength to consider.
#' @param n_grid Integer >= 2. Number of grid points (default 21).
#' @param baseline_logit,score Passed to \code{\link{correct_detection_bias}}.
#'
#' @return A data.frame with one row per delta and columns delta, var_null,
#'   vpc_null, var_main, vpc_main, pcv, observed_prevalence,
#'   implied_true_prevalence, min_relative_detection.
#'
#' @export
vpc_detection_bounds <- function(df, delta_max = 1.0, n_grid = 21L,
                                 baseline_logit = 2.0, score = NULL) {
  if (delta_max < 0) stop("delta_max must be non-negative")
  if (n_grid < 2) stop("n_grid must be at least 2")

  deltas <- seq(0, delta_max, length.out = n_grid)
  rows <- lapply(deltas, function(d) {
    r <- correct_detection_bias(df, d, baseline_logit = baseline_logit, score = score)
    r$baseline_logit <- NULL
    as.data.frame(r)
  })
  do.call(rbind, rows)
}

#' Detection-bias tipping point (E-value analogue for MAIHDA)
#'
#' Finds the smallest under-detection strength at which the corrected quantity
#' (default null-model VPC) first reaches \code{threshold}. Returns \code{NA} if
#' the threshold is not crossed within \code{[0, delta_max]}.
#'
#' @param df A data.frame as for \code{\link{correct_detection_bias}}.
#' @param threshold Numeric. Target value of \code{quantity}.
#' @param quantity Character. Which returned quantity to track (default
#'   \code{"vpc_null"}).
#' @param delta_max Numeric. Upper bound of the search (default 3.0).
#' @param baseline_logit,score Passed to \code{\link{correct_detection_bias}}.
#'
#' @return Numeric scalar: the tipping-point delta, or \code{NA_real_}.
#'
#' @export
detection_tipping_point <- function(df, threshold, quantity = "vpc_null",
                                    delta_max = 3.0, baseline_logit = 2.0,
                                    score = NULL) {
  f <- function(d) {
    correct_detection_bias(df, d, baseline_logit = baseline_logit,
                           score = score)[[quantity]] - threshold
  }
  f0 <- f(0)
  if (isTRUE(f0 == 0)) return(0)
  f_max <- f(delta_max)
  if (is.na(f0) || is.na(f_max) || sign(f0) == sign(f_max)) return(NA_real_)
  stats::uniroot(f, lower = 0, upper = delta_max)$root
}
