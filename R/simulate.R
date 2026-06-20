# simulate.R — Synthetic intersectional data generation for I-MAIHDA stress tests
#
# Replicates the logic of imaihda_sim/simulate.py in R.
# Uses Mersenne Twister RNG (set.seed), so numerical output differs from
# Python's PCG64, but statistical patterns and benchmark pass/fail are equivalent.
#
# Reference: Evans CR et al. (2018) SSM-PH; O'Sullivan JL et al. (2024) DMS

# ── Core helpers ──────────────────────────────────────────────────────────────

logit_inv <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(p) {
  log(p / (1 - p))
}

# ── Stratum table: 2(sex) × 3(education) × 3(wealth) × 2(rural) = 36 ────────

stratum_table <- function() {
  expand.grid(
    sex       = c(0L, 1L),
    education = c(0L, 1L, 2L),   # 0 = high, 2 = low
    wealth    = c(0L, 1L, 2L),   # 0 = high, 2 = low
    rural     = c(0L, 1L)        # 1 = rural / less-resourced
  )
}

# ── Main simulation ───────────────────────────────────────────────────────────

#' Simulate individual-level binary outcomes nested in intersectional strata
#'
#' @param n Number of individuals (default 6000)
#' @param prevalence_shift Intercept on logit scale (default -2.10)
#' @param interaction_sd SD of residual intersectional random effects (default 0)
#' @param detection_strength Strength of SES-patterned under-detection (default 0)
#' @param sparse If TRUE, stratum sizes follow a skewed gamma distribution
#' @param seed RNG seed
#' @return data.frame with columns: sex, education, wealth, rural, stratum,
#'         y, y_true, detection_probability, true_stratum_residual
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
    weights <- rgamma(k, shape = 0.35, scale = 1.0)
    weights <- weights / sum(weights)
  } else {
    weights <- rep(1 / k, k)
  }

  stratum_id <- sample(seq_len(k), size = n, replace = TRUE, prob = weights)
  x <- strata[stratum_id, , drop = FALSE]
  rownames(x) <- NULL

  # Additive data-generating component
  linear_predictor <- prevalence_shift +
    0.20 * x$sex +
    0.35 * x$education +
    0.30 * x$wealth +
    0.25 * x$rural

  # True residual intersectional component
  true_stratum_residual <- rep(0, k)
  if (interaction_sd > 0) {
    set.seed(seed + 999L)
    raw <- rnorm(k, mean = 0, sd = interaction_sd)
    # Structured heterogeneity: education==2 & wealth==2 & rural==1
    idx_high <- which(strata$education == 2 & strata$wealth == 2 & strata$rural == 1)
    raw[idx_high] <- raw[idx_high] + 0.90
    # Structured heterogeneity: education==0 & wealth==2 & rural==0
    idx_low <- which(strata$education == 0 & strata$wealth == 2 & strata$rural == 0)
    raw[idx_low] <- raw[idx_low] - 0.60
    true_stratum_residual <- raw - mean(raw)
    linear_predictor <- linear_predictor + true_stratum_residual[stratum_id]
  }

  p_true <- logit_inv(linear_predictor)
  y_true <- rbinom(n, size = 1L, prob = p_true)

  # Differential detection: true cases less likely recorded in disadvantaged strata
  if (detection_strength > 0) {
    detect_logit <- 2.0 -
      detection_strength * x$education -
      detection_strength * x$wealth -
      0.40 * detection_strength * x$rural
    detect_p <- logit_inv(detect_logit)
    detected <- rbinom(n, size = 1L, prob = detect_p)
    y_observed <- y_true * detected
  } else {
    detect_p <- rep(1, n)
    y_observed <- y_true
  }

  # Assemble output
  out <- x
  out$stratum                 <- factor(stratum_id, levels = seq_len(k))
  out$y                       <- as.integer(y_observed)
  out$y_true                  <- as.integer(y_true)
  out$detection_probability   <- detect_p
  out$true_stratum_residual   <- true_stratum_residual[stratum_id]
  out
}

# ── Scenario definitions ──────────────────────────────────────────────────────

#' Return a named list of scenario parameter sets
scenario_grid <- function(seed = 42L) {
  list(
    A = list(
      name        = "A",
      description = "Additive social gradient, equal detection",
      seed        = seed
    ),
    B = list(
      name            = "B",
      description     = "Residual intersectional heterogeneity, equal detection",
      interaction_sd  = 0.90,
      seed            = seed
    ),
    C = list(
      name               = "C",
      description        = "Additive structure with SES-patterned under-detection",
      detection_strength = 0.80,
      seed               = seed
    ),
    D = list(
      name               = "D",
      description        = "Residual interaction plus SES-patterned under-detection",
      interaction_sd     = 0.90,
      detection_strength = 0.80,
      seed               = seed
    ),
    E = list(
      name             = "E",
      description      = "Residual interaction with rare outcome and sparse strata",
      n                = 3500L,
      prevalence_shift = -3.00,
      interaction_sd   = 0.90,
      sparse           = TRUE,
      seed             = seed
    )
  )
}
