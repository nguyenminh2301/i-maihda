# fit.R — Main MAIHDA fitting function with dual estimators
#
#' Fit I-MAIHDA diagnostics (fast or full GLMM)
#'
#' Estimates null and main-effects between-stratum variance and computes
#' VPC and PCV. Two estimation methods are available:
#'
#' \describe{
#'   \item{\code{"fast"} (default)}{Empirical-logit method-of-moments
#'         estimator. Subtracts expected binomial sampling noise from
#'         the weighted variance of stratum-level residuals. Fast enough
#'         for repeated stress-testing. Approximate.}
#'   \item{\code{"glmer"}}{Full random-intercept logistic GLMM via
#'         \code{lme4::glmer}. Extracts between-stratum variance from
#'         the estimated random-effects variance component. Gold-standard
#'         for empirical research. Equivalent to CRAN \code{MAIHDA}.}
#' }
#'
#' @param df A data.frame with columns:
#'   \describe{
#'     \item{stratum}{Factor. Intersectional stratum membership.}
#'     \item{sex}{Factor or integer (0/1).}
#'     \item{education}{Factor or integer (0/1/2).}
#'     \item{wealth}{Factor or integer (0/1/2).}
#'     \item{rural}{Factor or integer (0/1).}
#'     \item{y}{Integer (0/1). Binary outcome.}
#'   }
#' @param method Character. \code{"fast"} for method-of-moments
#'   diagnostic (default), or \code{"glmer"} for full GLMM via lme4.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{n, n_strata}{Sample size and number of strata.}
#'     \item{min_stratum_n, median_stratum_n, max_stratum_n}{Stratum size distribution.}
#'     \item{overall_prevalence}{Mean of y.}
#'     \item{var_null, vpc_null}{Null-model between-stratum variance and VPC.}
#'     \item{var_main, vpc_main}{Main-effects residual variance and VPC.}
#'     \item{pcv}{Proportional Change in Variance (\%).}
#'     \item{estimator}{Which method was used.}
#'     \item{warnings}{Methodological caveats.}
#'   }
#'
#' @references
#' Evans CR, Williams DR, Onnela JP, Subramanian SV (2018).
#' "A multilevel approach to modeling health inequalities at the
#' intersection of multiple social identities."
#' \emph{SSM - Population Health}, 6, 149--157.
#'
#' Bulut H (2026). \emph{MAIHDA: Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy}. R package version 0.1.11.
#' \url{https://CRAN.R-project.org/package=MAIHDA}
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 2000, seed = 123)
#'
#' # Fast diagnostic (default, suitable for simulation stress-tests)
#' res_fast <- fit_imaihda(df)
#'
#' # Full GLMM estimator (for empirical research)
#' res_glmer <- fit_imaihda(df, method = "glmer")
#' }
#'
#' @export
fit_imaihda <- function(df, method = c("fast", "glmer")) {
  method <- match.arg(method)

  # Validate input
  required_cols <- c("stratum", "sex", "education", "wealth", "rural", "y")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (length(unique(df$stratum)) < 2) {
    stop("Need at least 2 strata; found ", length(unique(df$stratum)))
  }
  if (!all(df$y %in% c(0, 1))) {
    stop("y must be binary (0/1)")
  }

  if (method == "fast") {
    res <- fit_imaihda_fast(df)
  } else {
    res <- fit_imaihda_glmer(df)
  }

  res
}

# ── Fast method-of-moments estimator ──────────────────────────────────────────

fit_imaihda_fast <- function(df) {
  strata <- aggregate_strata(df)
  logit_v <- strata$empirical_logit
  sampling_var <- strata$logit_sampling_var
  weights <- strata$precision_weight

  # Null model: overall logit
  p_overall <- mean(df$y)
  p_overall <- min(max(p_overall, 1e-6), 1 - 1e-6)
  null_pred <- rep(log(p_overall / (1 - p_overall)), nrow(strata))
  var0 <- between_stratum_variance(logit_v, null_pred, sampling_var, weights)

  # Main-effects additive GLM
  main <- stats::glm(
    y ~ factor(sex) + factor(education) + factor(wealth) + factor(rural),
    data    = df,
    family  = stats::binomial()
  )
  main_pred <- stats::predict(main, newdata = strata, type = "link")
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
      warnings  = paste(
        "diagnostic approximation;",
        "use method='glmer' or full random-intercept MAIHDA for empirical analysis"
      )
    )
  )
}

# ── Full GLMM estimator (lme4::glmer) ─────────────────────────────────────────

#' Extract between-stratum random-intercept variance from a glmer fit
#'
#' Mirrors \code{MAIHDA:::maihda_stratum_variance_lme4()}.
#'
#' @param model A fitted \code{\link[lme4]{glmer}} model.
#' @param group Character. Name of the random-effect grouping factor.
#' @return Numeric scalar: intercept variance.
#' @keywords internal
extract_glmer_stratum_variance <- function(model, group = "stratum") {
  vc <- lme4::VarCorr(model)
  if (!group %in% names(vc)) {
    stop("No '", group, "' random-effect variance found in the model.")
  }
  group_vc <- as.matrix(vc[[group]])
  effect_names <- rownames(group_vc)
  intercept_name <- intersect(c("(Intercept)", "Intercept"), effect_names)
  if (length(intercept_name) == 0) {
    stop("The '", group,
         "' random effect must include an intercept for MAIHDA variance calculations.")
  }
  as.numeric(group_vc[intercept_name[1], intercept_name[1]])
}

fit_imaihda_glmer <- function(df) {
  # Ensure stratum is a factor
  df$stratum <- as.factor(df$stratum)

  # Null model: intersectional strata only (no fixed effects)
  null_model <- lme4::glmer(
    y ~ (1 | stratum),
    data    = df,
    family  = stats::binomial()
  )

  # Main-effects model: additive social determinants + random intercept
  main_model <- lme4::glmer(
    y ~ factor(sex) + factor(education) + factor(wealth) + factor(rural) +
      (1 | stratum),
    data    = df,
    family  = stats::binomial()
  )

  var0 <- extract_glmer_stratum_variance(null_model)
  var1 <- extract_glmer_stratum_variance(main_model)

  # Diagnostics from aggregated strata
  strata <- aggregate_strata(df)
  diag <- stratum_diagnostics(df, strata)

  c(
    diag,
    list(
      var_null  = var0,
      vpc_null  = vpc_latent(var0),
      var_main  = var1,
      vpc_main  = vpc_latent(var1),
      pcv       = pcv(var0, var1),
      estimator = "lme4_glmer_full_GLMM",
      warnings  = paste(
        "Full GLMM via lme4::glmer.",
        "Equivalent to CRAN MAIHDA package.",
        "Verify convergence with summary(model)."
      )
    )
  )
}
