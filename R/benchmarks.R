# benchmarks.R — Scenario-level pass/fail checks for the synthetic demo
#
# Replicates imaihda_sim/benchmarks.py

evaluate_benchmarks <- function(results) {
  # results: data.frame with at least columns scenario, vpc_null, vpc_main, pcv,
  #          overall_prevalence, min_stratum_n
  r <- results
  rownames(r) <- r$scenario

  checks <- list()

  add <- function(name, passed, value, criterion) {
    checks[[length(checks) + 1]] <<- list(
      check     = name,
      passed    = passed,
      value     = value,
      criterion = criterion
    )
  }

  add(
    "A_additive_is_additive_dominant",
    !is.na(r["A", "pcv"]) && r["A", "pcv"] >= 80 && r["A", "vpc_main"] < 1.0,
    sprintf("PCV=%.1f; VPC_main=%.2f", r["A", "pcv"], r["A", "vpc_main"]),
    "PCV >= 80 and VPC_main < 1"
  )

  add(
    "B_interaction_increases_vpc",
    r["B", "vpc_null"] > r["A", "vpc_null"] + 5,
    sprintf("A=%.2f; B=%.2f", r["A", "vpc_null"], r["B", "vpc_null"]),
    "B VPC_null > A VPC_null + 5 percentage points"
  )

  add(
    "B_interaction_leaves_residual_variance",
    !is.na(r["B", "pcv"]) && r["B", "pcv"] < 70,
    sprintf("PCV=%.1f", r["B", "pcv"]),
    "PCV < 70"
  )

  add(
    "C_detection_reduces_observed_prevalence",
    r["C", "overall_prevalence"] < r["A", "overall_prevalence"],
    sprintf("A=%.1f%%; C=%.1f%%", 100 * r["A", "overall_prevalence"], 100 * r["C", "overall_prevalence"]),
    "C observed prevalence < A observed prevalence"
  )

  add(
    "D_detection_can_mask_interaction_vpc",
    r["D", "vpc_null"] < r["B", "vpc_null"],
    sprintf("B=%.2f; D=%.2f", r["B", "vpc_null"], r["D", "vpc_null"]),
    "D VPC_null < B VPC_null despite same residual-interaction SD"
  )

  add(
    "E_sparse_strata_are_flagged",
    r["E", "min_stratum_n"] < r["B", "min_stratum_n"],
    sprintf("B min_n=%.0f; E min_n=%.0f", r["B", "min_stratum_n"], r["E", "min_stratum_n"]),
    "E minimum stratum size < B minimum stratum size"
  )

  do.call(rbind, lapply(checks, as.data.frame, stringsAsFactors = FALSE))
}
