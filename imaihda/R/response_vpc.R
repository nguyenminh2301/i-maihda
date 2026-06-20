# response_vpc.R — Response-scale VPC for binary MAIHDA
#
#' Response-Scale (Probability) VPC for Binary MAIHDA
#'
#' Computes the VPC on the probability (response) scale for a binary
#' outcome. Unlike the latent-scale VPC (\code{vpc_latent()}) which uses
#' the logistic variance \eqn{\pi^2/3}, the response-scale VPC is
#' computed by simulation and represents the proportion of variance in
#' observed probabilities that lies between strata.
#'
#' @param model For \code{method="glmer"}: a fitted \code{lme4::glmer} object.
#'   For \code{method="fast"}: a result list from \code{fit_imaihda()}.
#' @param method Character. \code{"fast"} for delta-method approximation,
#'   or \code{"glmer"} for simulation-based (matches CRAN MAIHDA).
#' @param n_sim Integer. Number of simulation draws (for glmer method only).
#' @param seed Integer. RNG seed for reproducibility.
#'
#' @return Numeric scalar: VPC on the probability scale (0--1).
#'
#' @details
#' The latent-scale VPC (logistic) is the standard in I-MAIHDA because
#' it is parameterization-invariant. The response-scale VPC is sometimes
#' preferred for substantive interpretation, but it depends on the
#' overall prevalence and the fixed-effect coefficients. This function
#' provides both the fast delta-method approximation and the full
#' simulation-based estimator (matching CRAN MAIHDA's
#' \code{maihda_vpc_response()}).
#'
#' @references
#' Goldstein H, Browne W, Rasbash J (2002). "Partitioning variation in
#' multilevel models." \emph{Understanding Statistics}, 1(4), 223–231.
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 3000, interaction_sd = 0.9, seed = 42)
#' res <- fit_imaihda(df)
#'
#' # Fast approximation
#' response_vpc(res, method = "fast")
#'
#' # Simulation-based (matches CRAN MAIHDA)
#' fit <- lme4::glmer(y ~ (1 | stratum), data = df, family = binomial())
#' response_vpc(fit, method = "glmer", n_sim = 5000, seed = 42)
#' }
#'
#' @export
response_vpc <- function(model, method = c("fast", "glmer"),
                         n_sim = 5000L, seed = NULL) {
  method <- match.arg(method)

  if (method == "fast") {
    return(response_vpc_fast(model))
  } else {
    return(response_vpc_glmer(model, n_sim, seed))
  }
}

# ── Fast: delta-method approximation ─────────────────────────────────────────

response_vpc_fast <- function(res) {
  # res is a fit_imaihda() result list
  sigma2 <- res[["var_null"]]
  p <- res[["overall_prevalence"]]

  # Delta-method: Var(p) ≈ p^2 * (1-p)^2 * sigma^2_stratum
  # This is the first-order approximation of the inverse-logit transform
  between_var <- (p * (1 - p))^2 * sigma2

  # Within-stratum (binomial) variance at the probability scale
  within_var <- p * (1 - p) / res[["n"]] * res[["n_strata"]]
  # Actually: average individual-level binomial variance
  within_var <- p * (1 - p)

  vpc_resp <- between_var / (between_var + within_var)
  vpc_resp
}

# ── GLMM: simulation-based (matches CRAN MAIHDA) ─────────────────────────────

response_vpc_glmer <- function(model, n_sim = 5000L, seed = NULL) {
  if (!inherits(model, "glmerMod")) {
    stop("method='glmer' requires a fitted lme4::glmer object")
  }

  if (!is.null(seed)) set.seed(seed)

  # Extract fixed effects and random effects variance
  fe <- lme4::fixef(model)
  re_names <- names(lme4::VarCorr(model))
  stratum_group <- re_names[1]
  sigma2 <- extract_glmer_stratum_variance(model, group = stratum_group)
  sigma <- sqrt(sigma2)

  # Draw random effects from N(0, sigma^2)
  re_draws <- stats::rnorm(n_sim, mean = 0, sd = sigma)

  # Compute predicted probabilities for each draw
  # Use the intercept-only case (null model)
  intercept <- fe[["(Intercept)"]]
  probs <- plogis(intercept + re_draws)

  # Between-stratum variance of probabilities
  between_var <- stats::var(probs)

  # Within-stratum binomial variance (average across draws)
  within_var <- mean(probs * (1 - probs))

  vpc_resp <- between_var / (between_var + within_var)
  vpc_resp
}
