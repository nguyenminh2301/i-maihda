# simulate.R — Synthetic intersectional data generation for I-MAIHDA testing
#
#' Build intersectional stratum table
#'
#' Creates the full factorial grid of 36 intersectional strata
#' defined by sex (2) × education (3) × wealth (3) × rural (2).
#'
#' @return A data.frame with 36 rows and columns: sex, education, wealth, rural.
#'
#' @keywords internal
stratum_table <- function() {
  expand.grid(
    sex       = c(0L, 1L),
    education = c(0L, 1L, 2L),   # 0 = high, 2 = low
    wealth    = c(0L, 1L, 2L),   # 0 = high, 2 = low
    rural     = c(0L, 1L)        # 1 = rural / less-resourced
  )
}

#' Simulate intersectional binary outcome data
#'
#' Generates synthetic individual-level data nested in intersectional
#' strata for stress-testing I-MAIHDA diagnostics. The data are
#' synthetic and make no empirical claim about any population.
#'
#' The data-generating process:
#' \enumerate{
#'   \item Assign individuals to 36 intersectional strata
#'         (sex × education × wealth × rural).
#'   \item Compute additive linear predictor with configurable
#'         prevalence shift and social gradients.
#'   \item Optionally add residual intersectional heterogeneity
#'         (structured interaction effects not reducible to main effects).
#'   \item Optionally apply SES-patterned under-detection
#'         (true cases less likely recorded in disadvantaged strata).
#'   \item Optionally use sparse stratum allocation
#'         (gamma-distributed stratum sizes).
#' }
#'
#' @param n Integer. Number of individuals (default 6000).
#' @param prevalence_shift Numeric. Intercept on logit scale (default -2.10).
#' @param interaction_sd Numeric. SD of residual intersectional random
#'   effects. 0 = purely additive (default).
#' @param detection_strength Numeric. Strength of SES-patterned
#'   under-detection. 0 = equal detection (default).
#' @param sparse Logical. If TRUE, stratum sizes follow a skewed gamma
#'   distribution instead of equal allocation.
#' @param seed Integer. RNG seed for reproducibility (default 42).
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{sex, education, wealth, rural}{Social dimensions (factors).}
#'     \item{stratum}{Factor. Intersectional stratum ID.}
#'     \item{y}{Integer (0/1). Observed binary outcome.}
#'     \item{y_true}{Integer (0/1). True outcome before detection bias.}
#'     \item{detection_probability}{Numeric. Individual detection probability.}
#'     \item{true_stratum_residual}{Numeric. True residual interaction effect
#'           for the individual's stratum.}
#'   }
#'
#' @references
#' Evans CR, Williams DR, Onnela JP, Subramanian SV (2018).
#' "A multilevel approach to modeling health inequalities at the
#' intersection of multiple social identities."
#' \emph{SSM - Population Health}, 6, 149--157.
#'
#' @examples
#' # Additive scenario
#' df <- simulate_intersectional_data(n = 1000, seed = 123)
#' table(df$stratum)
#'
#' # With residual intersectional heterogeneity
#' df_b <- simulate_intersectional_data(n = 1000, interaction_sd = 0.9, seed = 123)
#'
#' # With detection bias
#' df_c <- simulate_intersectional_data(n = 1000, detection_strength = 0.8, seed = 123)
#' mean(df_c$y)       # observed prevalence
#' mean(df_c$y_true)  # true prevalence (higher)
#'
#' @export
simulate_intersectional_data <- function(
    n                  = 6000L,
    prevalence_shift   = -2.10,
    interaction_sd     = 0.0,
    detection_strength = 0.0,
    sparse             = FALSE,
    seed               = 42L
) {
  set.seed(seed)
  strata <- stratum_table()
  k <- nrow(strata)

  # Stratum allocation weights
  if (sparse) {
    weights <- stats::rgamma(k, shape = 0.35, scale = 1.0)
    weights <- weights / sum(weights)
  } else {
    weights <- rep(1 / k, k)
  }

  stratum_id <- sample(seq_len(k), size = n, replace = TRUE, prob = weights)

  # Build output efficiently — only social dimension columns + stratum + y
  out <- data.frame(
    sex       = strata$sex[stratum_id],
    education = strata$education[stratum_id],
    wealth    = strata$wealth[stratum_id],
    rural     = strata$rural[stratum_id],
    stratum   = as.integer(stratum_id),
    stringsAsFactors = FALSE
  )

  # Additive linear predictor
  eta <- prevalence_shift +
    0.20 * out$sex +
    0.35 * out$education +
    0.30 * out$wealth +
    0.25 * out$rural

  # Residual intersectional heterogeneity
  if (interaction_sd > 0) {
    set.seed(seed + 999L)
    raw <- stats::rnorm(k, mean = 0, sd = interaction_sd)

    idx_high <- which(
      strata$education == 2 & strata$wealth == 2 & strata$rural == 1
    )
    raw[idx_high] <- raw[idx_high] + 0.90

    idx_low <- which(
      strata$education == 0 & strata$wealth == 2 & strata$rural == 0
    )
    raw[idx_low] <- raw[idx_low] - 0.60

    stratum_residual <- raw - mean(raw)
    out$true_stratum_residual <- stratum_residual[stratum_id]
    eta <- eta + out$true_stratum_residual
  }

  # True outcomes
  y_true <- stats::rbinom(n, size = 1L, prob = logit_inv(eta))

  # Differential detection
  if (detection_strength > 0) {
    detect_logit <- 2.0 -
      detection_strength * out$education -
      detection_strength * out$wealth -
      0.40 * detection_strength * out$rural
    detect_p <- logit_inv(detect_logit)
    detected <- stats::rbinom(n, size = 1L, prob = detect_p)
    out$y <- as.integer(y_true * detected)
    out$y_true <- as.integer(y_true)
    out$detection_probability <- detect_p
  } else {
    out$y <- as.integer(y_true)
    out$y_true <- out$y  # same as observed when no detection bias
  }

  out
}
