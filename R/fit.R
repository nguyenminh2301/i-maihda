# fit.R — Fast synthetic I-MAIHDA diagnostics (VPC, PCV) for benchmarking
#
# Replicates the logic of imaihda_sim/fit.py in R.
# Uses empirical-logit method-of-moments estimators, NOT full GLMM.
# For empirical PhD work, use lme4::glmer with random intercepts.

LOGISTIC_L1_VARIANCE <- pi^2 / 3

# ── Core formulas ─────────────────────────────────────────────────────────────

#' Latent-response VPC for logistic MAIHDA
vpc_latent <- function(stratum_variance) {
  if (stratum_variance < 0) stop("stratum_variance must be non-negative")
  100 * stratum_variance / (stratum_variance + LOGISTIC_L1_VARIANCE)
}

#' Proportional change in stratum-level variance (null → main-effects)
pcv <- function(var_null, var_main) {
  if (var_null <= 0) return(NaN)
  100 * (var_null - var_main) / var_null
}

# ── Aggregation and variance estimation ───────────────────────────────────────

#' Aggregate individual-level data to stratum-level empirical logits
#' Uses Laplace smoothing (events + 0.5) / (n + 1)
aggregate_strata <- function(df) {
  g <- aggregate(
    df$y,
    by = list(
      stratum    = df$stratum,
      sex        = df$sex,
      education  = df$education,
      wealth     = df$wealth,
      rural      = df$rural
    ),
    FUN = function(x) c(events = sum(x), n = length(x))
  )
  # Flatten matrix column
  g$events <- g$x[, "events"]
  g$n      <- g$x[, "n"]
  g$x      <- NULL

  g$empirical_p     <- (g$events + 0.5) / (g$n + 1)
  g$empirical_logit <- log(g$empirical_p / (1 - g$empirical_p))
  g$logit_sampling_var <- 1 / (g$events + 0.5) + 1 / (g$n - g$events + 0.5)
  g$precision_weight   <- 1 / g$logit_sampling_var
  g
}

#' Weighted variance: sum(w * (x - weighted.mean(x, w))^2) / sum(w)
weighted_variance <- function(x, w) {
  mu <- weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

#' Fast diagnostic estimate of residual between-stratum logit variance
#'
#' Subtracts expected binomial sampling noise from the weighted variance
#' of stratum-level residuals. Method-of-moments estimator.
between_stratum_variance <- function(logit_vec, pred_vec, sampling_var, weights) {
  residual  <- logit_vec - pred_vec
  raw_var   <- weighted_variance(residual, weights)
  expected_noise <- weighted.mean(sampling_var, weights)
  max(raw_var - expected_noise, 0)
}

#' Stratum-level diagnostics
stratum_diagnostics <- function(df, strata) {
  prev <- tapply(df$y, df$stratum, mean)
  out <- list(
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
  out
}

# ── Main fitting function ─────────────────────────────────────────────────────

#' Fit fast synthetic I-MAIHDA diagnostics
#'
#' Estimates null and main-effects between-stratum variance using
#' empirical stratum logits after subtracting binomial sampling variance.
#' Designed for simulation and benchmarking only.
fit_imaihda <- function(df) {
  strata   <- aggregate_strata(df)
  logit_v  <- strata$empirical_logit
  sampling_var <- strata$logit_sampling_var
  weights  <- strata$precision_weight

  # Null model: overall logit
  p_overall <- mean(df$y)
  p_overall <- min(max(p_overall, 1e-6), 1 - 1e-6)
  null_pred <- rep(log(p_overall / (1 - p_overall)), nrow(strata))
  var0 <- between_stratum_variance(logit_v, null_pred, sampling_var, weights)

  # Main-effects model: additive main effects only
  main <- glm(
    y ~ factor(sex) + factor(education) + factor(wealth) + factor(rural),
    data    = df,
    family  = binomial()
  )
  main_pred <- predict(main, newdata = strata, type = "link")
  var1 <- between_stratum_variance(logit_v, main_pred, sampling_var, weights)

  diag <- stratum_diagnostics(df, strata)
  c(
    diag,
    list(
      var_null  = var0,
      vpc_null  = vpc_latent(var0),
      var_main  = var1,
      vpc_main  = vpc_latent(var1),
      pcv       = pcv(var0, var1),
      estimator = "fast_empirical_logit_diagnostic",
      warnings  = "diagnostic approximation; use full random-intercept MAIHDA for empirical analysis"
    )
  )
}

#' Run a single scenario end-to-end
fit_scenario <- function(name, scenario) {
  # Use [[ for exact match (avoid R $ partial matching: $n matches $name)
  args <- list(
    n                  = if (is.null(scenario[["n"]])) 6000L else scenario[["n"]],
    prevalence_shift   = if (is.null(scenario[["prevalence_shift"]])) -2.10 else scenario[["prevalence_shift"]],
    interaction_sd     = if (is.null(scenario[["interaction_sd"]])) 0.0 else scenario[["interaction_sd"]],
    detection_strength = if (is.null(scenario[["detection_strength"]])) 0.0 else scenario[["detection_strength"]],
    sparse             = if (is.null(scenario[["sparse"]])) FALSE else scenario[["sparse"]],
    seed               = if (is.null(scenario[["seed"]])) 42L else scenario[["seed"]]
  )
  df <- do.call(simulate_intersectional_data, args)
  res <- fit_imaihda(df)
  res$scenario          <- name
  res$description       <- scenario$description
  res$interaction_sd    <- if (is.null(scenario$interaction_sd)) 0.0 else scenario$interaction_sd
  res$detection_strength <- if (is.null(scenario$detection_strength)) 0.0 else scenario$detection_strength
  res$sparse            <- if (is.null(scenario$sparse)) FALSE else scenario$sparse
  res
}
