# scenarios.R — Predefined test scenarios and benchmark validation
#
#' Predefined MAIHDA stress-test scenarios
#'
#' Returns a list of five parameter sets (A–E) designed to test how
#' VPC and PCV respond to additive social gradients, residual
#' intersectional heterogeneity, SES-patterned under-detection,
#' and sparse strata.
#'
#' @param seed Integer. RNG seed for reproducibility (default 42).
#'
#' @return A named list of lists, each with parameters for
#'   \code{\link{simulate_intersectional_data}}.
#'
#' @details
#' \describe{
#'   \item{A}{Additive social gradient, equal detection.}
#'   \item{B}{Residual intersectional heterogeneity, equal detection.}
#'   \item{C}{Additive structure with SES-patterned under-detection.}
#'   \item{D}{Residual interaction plus SES-patterned under-detection.}
#'   \item{E}{Residual interaction with rare outcome and sparse strata.}
#' }
#'
#' @export
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

#' Run a single scenario end-to-end
#'
#' Simulates data for one scenario and fits MAIHDA diagnostics.
#'
#' @param name Character. Scenario label (e.g., "A").
#' @param scenario List. Scenario parameters from \code{\link{scenario_grid}}.
#'
#' @return A named list of results (same format as \code{\link{fit_imaihda}})
#'   with additional scenario metadata.
#'
#' @keywords internal
fit_scenario <- function(name, scenario) {
  args <- list(
    n                  = if (is.null(scenario[["n"]])) 6000L else scenario[["n"]],
    prevalence_shift   = if (is.null(scenario[["prevalence_shift"]])) -2.10
                         else scenario[["prevalence_shift"]],
    interaction_sd     = if (is.null(scenario[["interaction_sd"]])) 0.0
                         else scenario[["interaction_sd"]],
    detection_strength = if (is.null(scenario[["detection_strength"]])) 0.0
                         else scenario[["detection_strength"]],
    sparse             = if (is.null(scenario[["sparse"]])) FALSE
                         else scenario[["sparse"]],
    seed               = if (is.null(scenario[["seed"]])) 42L
                         else scenario[["seed"]]
  )
  df  <- do.call(simulate_intersectional_data, args)
  res <- fit_imaihda(df)
  res$scenario          <- name
  res$description       <- scenario$description
  res$interaction_sd    <- if (is.null(scenario[["interaction_sd"]])) 0.0
                           else scenario[["interaction_sd"]]
  res$detection_strength <- if (is.null(scenario[["detection_strength"]])) 0.0
                            else scenario[["detection_strength"]]
  res$sparse            <- if (is.null(scenario[["sparse"]])) FALSE
                           else scenario[["sparse"]]
  res
}

#' Evaluate methodological benchmark checks
#'
#' Runs six pass/fail checks against scenario results to verify
#' that the synthetic demo answers its intended question: VPC and
#' PCV are sensitive to prevalence, sparse strata, and differential
#' detection.
#'
#' @param results A data.frame with columns scenario, vpc_null, vpc_main,
#'   pcv, overall_prevalence, min_stratum_n (as returned by running
#'   all scenarios through \code{\link{fit_scenario}}).
#'
#' @return A data.frame with columns: check, passed, value, criterion.
#'
#' @export
evaluate_benchmarks <- function(results) {
  rownames(results) <- results$scenario

  checks <- list()

  add <- function(name, passed, value, criterion) {
    checks[[length(checks) + 1]] <<- data.frame(
      check     = name,
      passed    = passed,
      value     = value,
      criterion = criterion,
      stringsAsFactors = FALSE
    )
  }

  add(
    "A_additive_is_additive_dominant",
    !is.na(results["A", "pcv"]) && results["A", "pcv"] >= 80 &&
      results["A", "vpc_main"] < 1.0,
    sprintf("PCV=%.1f; VPC_main=%.2f",
            results["A", "pcv"], results["A", "vpc_main"]),
    "PCV >= 80 and VPC_main < 1"
  )

  add(
    "B_interaction_increases_vpc",
    results["B", "vpc_null"] > results["A", "vpc_null"] + 5,
    sprintf("A=%.2f; B=%.2f",
            results["A", "vpc_null"], results["B", "vpc_null"]),
    "B VPC_null > A VPC_null + 5 percentage points"
  )

  add(
    "B_interaction_leaves_residual_variance",
    !is.na(results["B", "pcv"]) && results["B", "pcv"] < 70,
    sprintf("PCV=%.1f", results["B", "pcv"]),
    "PCV < 70"
  )

  add(
    "C_detection_reduces_observed_prevalence",
    results["C", "overall_prevalence"] < results["A", "overall_prevalence"],
    sprintf("A=%.1f%%; C=%.1f%%",
            100 * results["A", "overall_prevalence"],
            100 * results["C", "overall_prevalence"]),
    "C observed prevalence < A observed prevalence"
  )

  add(
    "D_detection_can_mask_interaction_vpc",
    results["D", "vpc_null"] < results["B", "vpc_null"],
    sprintf("B=%.2f; D=%.2f",
            results["B", "vpc_null"], results["D", "vpc_null"]),
    "D VPC_null < B VPC_null despite same residual-interaction SD"
  )

  add(
    "E_sparse_strata_are_flagged",
    results["E", "min_stratum_n"] < results["B", "min_stratum_n"],
    sprintf("B min_n=%.0f; E min_n=%.0f",
            results["B", "min_stratum_n"], results["E", "min_stratum_n"]),
    "E minimum stratum size < B minimum stratum size"
  )

  do.call(rbind, checks)
}
