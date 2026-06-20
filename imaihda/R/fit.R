# fit.R — Main MAIHDA fitting function
#
#' Fit fast synthetic I-MAIHDA diagnostics
#'
#' Estimates null and main-effects between-stratum variance using empirical
#' stratum logits and a method-of-moments estimator that subtracts binomial
#' sampling noise. This is a fast diagnostic designed for simulation and
#' benchmarking, NOT a replacement for full random-intercept MAIHDA
#' (e.g., \code{lme4::glmer}) in empirical research.
#'
#' The function:
#' \enumerate{
#'   \item Aggregates individual-level data to intersectional strata
#'         with Laplace-smoothed empirical logits.
#'   \item Estimates null-model between-stratum variance (strata only).
#'   \item Fits an additive main-effects logistic GLM and estimates
#'         residual between-stratum variance after main effects.
#'   \item Computes VPC (null and main) and PCV.
#'   \item Reports stratum diagnostics (size, prevalence, detection).
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
#'   Optional columns for diagnostic reporting: \code{y_true},
#'   \code{detection_probability}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{n, n_strata}{Sample size and number of strata.}
#'     \item{min_stratum_n, median_stratum_n, max_stratum_n}{Stratum size distribution.}
#'     \item{overall_prevalence}{Mean of y.}
#'     \item{true_prevalence}{Mean of y_true (if available).}
#'     \item{mean_detection_probability}{Mean detection probability (if available).}
#'     \item{var_null}{Null-model between-stratum variance.}
#'     \item{vpc_null}{Null-model VPC (\%)}.}
#'     \item{var_main}{Main-effects residual between-stratum variance.}
#'     \item{vpc_main}{Main-effects VPC (\%)}.}
#'     \item{pcv}{Proportional Change in Variance (\%)}.}
#'     \item{estimator}{Always \code{"fast_empirical_logit_diagnostic"}.}
#'     \item{warnings}{Methodological caveat about diagnostic approximation.}
#'   }
#'
#' @references
#' Evans CR, Williams DR, Onnela JP, Subramanian SV (2018).
#' "A multilevel approach to modeling health inequalities at the
#' intersection of multiple social identities."
#' \emph{SSM - Population Health}, 6, 149--157.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' df <- simulate_intersectional_data(n = 2000, seed = 123)
#' res <- fit_imaihda(df)
#' print(res$vpc_null)
#' print(res$pcv)
#' }
#'
#' @export
fit_imaihda <- function(df) {
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

  # Aggregate to strata
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

  # Assemble results
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
        "use full random-intercept MAIHDA for empirical analysis"
      )
    )
  )
}
