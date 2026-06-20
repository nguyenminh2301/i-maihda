# aggregate.R — Stratum-level aggregation and diagnostic summary
#
# Internal functions for preparing intersectional data for MAIHDA
# variance estimation. Not exported.

#' Aggregate individual-level data to stratum-level empirical logits
#'
#' Groups individuals by intersectional stratum and computes empirical
#' logits with Laplace smoothing (events + 0.5) / (n + 1), along with
#' delta-method sampling variances and precision weights.
#'
#' @param df A data.frame with columns: stratum (factor), sex, education,
#'   wealth, rural (all factor or integer), and y (binary outcome 0/1).
#'
#' @return A data.frame with stratum-level statistics: events, n,
#'   empirical_p, empirical_logit, logit_sampling_var, precision_weight.
#'
#' @keywords internal
aggregate_strata <- function(df) {
  # Group by stratum and covariates
  g <- aggregate(
    df$y,
    by = list(
      stratum   = df$stratum,
      sex       = df$sex,
      education = df$education,
      wealth    = df$wealth,
      rural     = df$rural
    ),
    FUN = function(x) c(events = sum(x), n = length(x))
  )

  # Flatten matrix column from aggregate()
  g$events <- g$x[, "events"]
  g$n      <- g$x[, "n"]
  g$x      <- NULL

  # Laplace-smoothed empirical probabilities and logits
  g$empirical_p     <- (g$events + 0.5) / (g$n + 1)
  g$empirical_logit <- log(g$empirical_p / (1 - g$empirical_p))

  # Delta-method sampling variance of the empirical logit
  g$logit_sampling_var <- 1 / (g$events + 0.5) + 1 / (g$n - g$events + 0.5)

  # Precision (inverse-variance) weights
  g$precision_weight <- 1 / g$logit_sampling_var

  g
}

#' Compute stratum-level diagnostic summaries
#'
#' Reports sample sizes, prevalence, and detection statistics for
#' quality control of intersectional MAIHDA models.
#'
#' @param df A data.frame from simulate_intersectional_data() or
#'   similarly structured individual-level data.
#' @param strata A data.frame from aggregate_strata().
#'
#' @return A named list of diagnostic values: n, n_strata, min/median/max
#'   stratum sizes, overall/min/max prevalence, true prevalence (if
#'   y_true column present), and mean detection probability (if
#'   detection_probability column present).
#'
#' @keywords internal
stratum_diagnostics <- function(df, strata) {
  prev <- tapply(df$y, df$stratum, mean)

  list(
    n                        = nrow(df),
    n_strata                 = length(unique(df$stratum)),
    min_stratum_n            = as.integer(min(strata$n)),
    median_stratum_n         = as.numeric(median(strata$n)),
    max_stratum_n            = as.integer(max(strata$n)),
    overall_prevalence       = mean(df$y),
    min_stratum_prevalence   = as.numeric(min(prev)),
    max_stratum_prevalence   = as.numeric(max(prev)),
    true_prevalence          = if ("y_true" %in% names(df)) mean(df$y_true) else NaN,
    mean_detection_probability = if ("detection_probability" %in% names(df))
      mean(df$detection_probability) else NaN
  )
}
