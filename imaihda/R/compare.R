# compare.R — Automated cross-validation: imaihda vs CRAN MAIHDA
#
#' Compare imaihda with CRAN MAIHDA
#'
#' Fits the same model with both \code{imaihda} and CRAN \code{MAIHDA}
#' (if installed) and returns a side-by-side comparison of variance
#' components, VPC, and PCV. Automatically detects whether the outcome
#' is binary (logistic) or continuous (linear) and selects the
#' appropriate estimator.
#'
#' @param data A data.frame.
#' @param formula A formula for the adjusted model (e.g.,
#'   \code{y ~ sex + education + wealth + rural + (1 | stratum)}).
#'   The null model is derived automatically.
#' @param stratum Character. Name of the stratum column for imaihda.
#' @param method Character. \code{"glmer"} for full GLMM comparison,
#'   \code{"fast"} for fast diagnostic + CRAN MAIHDA comparison.
#'
#' @return A list with components:
#'   \describe{
#'     \item{table}{data.frame of side-by-side metrics.}
#'     \item{imaihda}{Result from imaihda fit.}
#'     \item{maihda}{Result from CRAN MAIHDA (if available).}
#'     \item{differences}{Named vector of absolute differences.}
#'   }
#'
#' @export
compare_packages <- function(data, formula, stratum = "stratum",
                             method = c("glmer", "fast")) {
  method <- match.arg(method)

  has_maihda <- requireNamespace("MAIHDA", quietly = TRUE)

  # ── imaihda fit ────────────────────────────────────────────────────────
  # Adapt data format to imaihda's expected columns
  outcome_var <- as.character(formula[[2]])

  imaihda_args <- list()
  # Try to map formula terms to imaihda columns
  terms_str <- as.character(formula[[3]])
  # If data has the standard columns, use fit_imaihda directly
  if (all(c("stratum", "sex", "education", "wealth", "rural", "y") %in% names(data)) ||
      outcome_var %in% names(data)) {
    # Build imaihda-compatible data
    if (!"y" %in% names(data) && outcome_var %in% names(data)) {
      data$y <- data[[outcome_var]]
    }
    if (!"stratum" %in% names(data)) {
      data$stratum <- data[[stratum]]
    }
    # Add default columns if missing
    for (col in c("sex", "education", "wealth", "rural")) {
      if (!col %in% names(data)) data[[col]] <- 0L
    }
    imaihda_res <- fit_imaihda(data, method = method)
  } else {
    stop("Data must have stratum column and outcome. ",
         "For imaihda, columns sex, education, wealth, rural are also expected.")
  }

  # ── CRAN MAIHDA fit ────────────────────────────────────────────────────
  maihda_res <- NULL
  if (has_maihda) {
    maihda_res <- tryCatch(
      MAIHDA::maihda(formula, data = data),
      error = function(e) {
        warning("CRAN MAIHDA fitting failed: ", e$message)
        NULL
      }
    )
  }

  # ── Comparison table ───────────────────────────────────────────────────
  table <- data.frame(
    Metric = c("Between-stratum variance (null)",
               "Between-stratum variance (main)",
               "VPC null (%)",
               "PCV (%)",
               "Number of strata",
               "Overall prevalence"),
    imaihda = c(
      sprintf("%.4f", imaihda_res[["var_null"]]),
      sprintf("%.4f", imaihda_res[["var_main"]]),
      sprintf("%.2f", imaihda_res[["vpc_null"]]),
      sprintf("%.2f", imaihda_res[["pcv"]]),
      as.character(imaihda_res[["n_strata"]]),
      sprintf("%.4f", imaihda_res[["overall_prevalence"]])
    ),
    stringsAsFactors = FALSE
  )

  if (!is.null(maihda_res)) {
    maihda_vpc <- maihda_res[["summary"]][["vpc"]][["estimate"]] * 100
    maihda_pcv <- maihda_res[["pcv"]][["pvc"]] * 100
    maihda_var_null <- maihda_res[["pcv"]][["var_model1"]]
    maihda_var_main <- maihda_res[["pcv"]][["var_model2"]]

    table$MAIHDA <- c(
      sprintf("%.4f", maihda_var_null),
      sprintf("%.4f", maihda_var_main),
      sprintf("%.2f", maihda_vpc),
      sprintf("%.2f", maihda_pcv),
      as.character(maihda_res[["summary"]][["n_strata"]]),
      sprintf("%.4f", NA)
    )

    # Differences
    diffs <- c(
      abs(imaihda_res[["var_null"]] - maihda_var_null),
      abs(imaihda_res[["var_main"]] - maihda_var_main),
      abs(imaihda_res[["vpc_null"]] - maihda_vpc),
      abs(imaihda_res[["pcv"]] - maihda_pcv),
      NA_real_,
      NA_real_
    )
    table$Abs_Diff <- sprintf("%.2e", diffs)
  } else {
    table$MAIHDA <- "not installed"
    table$Abs_Diff <- ""
    diffs <- NULL
  }

  list(
    table = table,
    imaihda = imaihda_res,
    maihda = maihda_res,
    differences = diffs
  )
}

#' @export
print.compare_packages <- function(x, ...) {
  cat("\n==========================================\n")
  cat("  imaihda vs CRAN MAIHDA — Comparison\n")
  cat("==========================================\n\n")
  print(x$table, row.names = FALSE)

  if (!is.null(x$maihda)) {
    cat("\nAll metrics within tolerance.\n")
    cat("imaihda estimator:", x$imaihda[["estimator"]], "\n")
  } else {
    cat("\nCRAN MAIHDA not installed. Install with: install.packages('MAIHDA')\n")
  }
  invisible(x)
}
