# discriminatory.R — Discriminatory accuracy metrics for binary MAIHDA
#
#' Discriminatory Accuracy for Binary MAIHDA
#'
#' Computes Area Under the ROC Curve (AUC) and Median Odds Ratio (MOR)
#' for a binary-outcome MAIHDA model. Supports two methods: fast
#' (approximate from predicted probabilities) and glmer (from fitted
#' lme4 model, matches CRAN MAIHDA).
#'
#' @param model A fitted model object. For \code{method="fast"}, a
#'   result list from \code{fit_imaihda()}. For \code{method="glmer"},
#'   a fitted \code{lme4::glmer} object.
#' @param data Optional data.frame with observed outcomes (if not
#'   extractable from model).
#' @param method Character. \code{"fast"} for approximate AUC from
#'   stratum-level predicted probabilities, or \code{"glmer"} for
#'   full individual-level predictions from the GLMM.
#'
#' @return A list with components:
#'   \describe{
#'     \item{AUC}{Area Under the ROC Curve (0.5 = no discrimination, 1 = perfect).}
#'     \item{MOR}{Median Odds Ratio. Compares two randomly chosen individuals
#'           from different strata. MOR = 1 indicates no between-stratum variation.}
#'   }
#'
#' @details
#' \strong{AUC} quantifies how well the model discriminates between
#' individuals with and without the outcome. For intersectional MAIHDA,
#' AUC reflects the strata's ability to separate cases from non-cases.
#'
#' \strong{MOR} is defined as
#' \eqn{\exp(\sqrt{2\sigma^2_{\text{stratum}}} \cdot \Phi^{-1}(0.75))},
#' where \eqn{\Phi^{-1}(0.75) \approx 0.6745}. It is the median of the
#' distribution of odds ratios comparing two individuals from
#' randomly selected strata.
#'
#' @references
#' Merlo J et al. (2006). "A brief conceptual tutorial of multilevel
#' analysis in social epidemiology." \emph{J Epidemiol Community Health},
#' 60(4), 290–297.
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 2000, interaction_sd = 0.9, seed = 42)
#' res <- fit_imaihda(df)
#' da <- discriminatory_accuracy(res, method = "fast")
#' print(da$AUC)
#' print(da$MOR)
#'
#' # Using glmer method (matches CRAN MAIHDA)
#' fit <- lme4::glmer(y ~ (1 | stratum), data = df, family = binomial())
#' da_glmer <- discriminatory_accuracy(fit, method = "glmer", data = df)
#' }
#'
#' @export
discriminatory_accuracy <- function(model, data = NULL, method = c("fast", "glmer")) {
  method <- match.arg(method)

  if (method == "fast") {
    return(discriminatory_accuracy_fast(model, data))
  } else {
    return(discriminatory_accuracy_glmer(model, data))
  }
}

# ── Fast method: approximate from aggregated stratum predictions ─────────────

discriminatory_accuracy_fast <- function(res, data = NULL) {
  # res is a fit_imaihda result list
  # We approximate AUC from stratum-level observed and predicted rates
  if (is.null(data) && !is.null(res[["n_strata"]])) {
    # Can't compute individual-level AUC without individual data
    # Use stratum-level approximation with precision weights
    p_overall <- res[["overall_prevalence"]]

    # Approximate AUC as 0.5 + 0.5 * (VPC-based discriminatory power)
    # This is a heuristic approximation; glmer method is preferred
    vpc <- res[["vpc_null"]] / 100
    auc_approx <- 0.5 + 0.4 * sqrt(vpc)  # Heuristic calibration
    auc_approx <- min(max(auc_approx, 0.5), 1.0)

    # MOR from stratum variance
    sigma2 <- res[["var_null"]]
    mor <- exp(sqrt(2 * sigma2) * stats::qnorm(0.75))

    return(list(
      AUC = auc_approx,
      MOR = mor,
      method = "fast_approximation",
      warning = paste(
        "Fast AUC is a heuristic approximation based on VPC.",
        "Use method='glmer' for exact AUC."
      )
    ))
  }

  # If individual data is available, compute ROC
  if (!is.null(data) && "y" %in% names(data)) {
    # Predict from aggregated strata
    y_true <- data$y
    # Simple approximation: use overall mean as predictor
    # (fast method doesn't fit individual-level model)
    y_pred <- rep(mean(y_true), length(y_true))
    auc <- 0.5  # no individual discrimination at fast level
    mor <- exp(sqrt(2 * res[["var_null"]]) * stats::qnorm(0.75))
    return(list(AUC = auc, MOR = mor, method = "fast_approximation",
                warning = "Fast method cannot compute individual-level AUC. Use method='glmer'."))
  }

  list(AUC = NA_real_, MOR = NA_real_, method = "fast_approximation",
       warning = "No individual data available for AUC computation.")
}

# ── GLMM method: from lme4::glmer (matches CRAN MAIHDA) ──────────────────────

#' @importFrom stats qnorm predict
discriminatory_accuracy_glmer <- function(model, data = NULL) {
  if (!inherits(model, "glmerMod")) {
    stop("method='glmer' requires a fitted lme4::glmer object")
  }

  # Extract stratum variance
  re_names <- names(lme4::VarCorr(model))
  stratum_group <- re_names[re_names != "stratum"]
  if (length(stratum_group) == 0) stratum_group <- re_names[1]
  sigma2 <- extract_glmer_stratum_variance(model, group = stratum_group)

  # MOR
  mor <- exp(sqrt(2 * sigma2) * stats::qnorm(0.75))

  # AUC: requires individual predictions
  if (!is.null(data) && "y" %in% names(data)) {
    y_true <- data$y
    y_pred <- stats::predict(model, type = "response")
    # Compute AUC using trapezoidal rule (no pROC dependency)
    auc <- .compute_auc(y_true, y_pred)
  } else {
    # Try extracting from model frame
    mf <- model@frame
    y_col <- names(mf)[1]
    if (y_col %in% names(mf)) {
      y_true <- mf[[y_col]]
      y_pred <- stats::predict(model, type = "response")
      auc <- .compute_auc(y_true, y_pred)
    } else {
      auc <- NA_real_
    }
  }

  list(AUC = auc, MOR = mor, method = "lme4_glmer")
}

# Simple AUC via trapezoidal rule (avoids pROC dependency)
.compute_auc <- function(y_true, y_pred) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  n <- length(y_true)
  ord <- order(y_pred, decreasing = TRUE)
  y_sorted <- y_true[ord]
  tpr <- cumsum(y_sorted) / sum(y_sorted)
  fpr <- cumsum(1 - y_sorted) / sum(1 - y_sorted)
  # Trapezoidal integration
  auc <- sum(diff(fpr) * (tpr[-1] + tpr[-n]) / 2)
  auc
}
