# stepwise.R — Stepwise variance decomposition for I-MAIHDA
#
#' Stepwise Proportional Change in Variance
#'
#' Adds covariates one-by-one and tracks how the between-stratum variance
#' changes at each step. Useful for identifying which social dimensions
#' contribute most to intersectional inequality. Supports two estimation
#' methods: fast (method-of-moments) and glmer (full GLMM).
#'
#' @param data A data.frame with a stratum column and outcome variable.
#' @param outcome Character. Name of the binary outcome column (0/1).
#' @param vars Character vector. Names of predictor columns to add stepwise.
#' @param stratum Character. Name of the stratum column (default "stratum").
#' @param method Character. \code{"fast"} for method-of-moments or
#'   \code{"glmer"} for full GLMM via lme4.
#' @param quiet Logical. Suppress step messages if TRUE.
#'
#' @return A data.frame with columns: Step, Model, Added_Variable,
#'   Variance, Step_PCV, Total_PCV.
#'
#' @details
#' The function starts with a null model (stratum random intercept only,
#' no fixed effects) and sequentially adds each predictor. At each step:
#' \itemize{
#'   \item \strong{Step_PCV}: proportion of remaining between-stratum variance
#'         explained by the newly added variable.
#'   \item \strong{Total_PCV}: cumulative proportion of the original
#'         between-stratum variance explained.
#' }
#' A negative Step_PCV indicates a "suppression" or "unmasking" pattern
#' (adding the variable \emph{increased} between-stratum variance).
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 3000, interaction_sd = 0.9, seed = 42)
#' stepwise_pcv(df, outcome = "y", vars = c("sex", "education", "wealth", "rural"))
#' }
#'
#' @export
stepwise_pcv <- function(data, outcome, vars, stratum = "stratum",
                         method = c("fast", "glmer"), quiet = FALSE) {
  method <- match.arg(method)

  if (!outcome %in% names(data)) stop("outcome '", outcome, "' not found in data")
  if (!stratum %in% names(data)) stop("stratum '", stratum, "' not found in data")
  for (v in vars) {
    if (!v %in% names(data)) stop("variable '", v, "' not found in data")
  }

  results <- list()

  # Step 0: null model (stratum only)
  if (!quiet) cat("Step 0: Null model (", stratum, " only)\n", sep = "")

  if (method == "fast") {
    # Compute null between-stratum variance directly (overall logit)
    strata_agg <- aggregate_strata(
      data.frame(stratum = data[[stratum]], sex = 0, education = 0,
                 wealth = 0, rural = 0, y = data[[outcome]])
    )
    p_overall <- mean(data[[outcome]])
    p_overall <- min(max(p_overall, 1e-6), 1 - 1e-6)
    null_pred <- rep(log(p_overall / (1 - p_overall)), nrow(strata_agg))
    var0 <- between_stratum_variance(
      strata_agg$empirical_logit, null_pred,
      strata_agg$logit_sampling_var, strata_agg$precision_weight
    )
  } else {
    null_fit <- lme4::glmer(
      stats::as.formula(paste(outcome, "~ (1 |", stratum, ")")),
      data = data, family = stats::binomial()
    )
    var0 <- extract_glmer_stratum_variance(null_fit, group = stratum)
  }

  results[[1]] <- data.frame(
    Step = 0L, Model = "Null Model",
    Added_Variable = "None (stratum only)",
    Variance = var0, Step_PCV = 0, Total_PCV = 0,
    stringsAsFactors = FALSE
  )

  # Steps 1..n: sequentially add variables
  current_vars <- character(0)
  for (i in seq_along(vars)) {
    current_vars <- c(current_vars, vars[i])
    if (!quiet) cat("Step ", i, ": Adding ", vars[i], "\n", sep = "")

    if (method == "fast") {
      # Fit additive GLM with current covariates, then compute
      # residual between-stratum variance
      formula_str <- paste0(
        "y ~ ", paste(paste0("factor(", current_vars, ")"), collapse = " + ")
      )
      df_step <- data
      df_step$y <- data[[outcome]]
      df_step$stratum <- as.factor(df_step[[stratum]])
      main_glm <- stats::glm(
        stats::as.formula(formula_str),
        data = df_step,
        family = stats::binomial()
      )
      strata_agg <- aggregate_strata(
        data.frame(stratum = data[[stratum]], sex = 0, education = 0,
                   wealth = 0, rural = 0, y = data[[outcome]])
      )
      main_pred <- stats::predict(main_glm, newdata = strata_agg, type = "link")
      vari <- between_stratum_variance(
        strata_agg$empirical_logit, main_pred,
        strata_agg$logit_sampling_var, strata_agg$precision_weight
      )
    } else {
      # glmer method
      rhs <- paste0(
        paste(current_vars, collapse = " + "),
        " + (1 | ", stratum, ")"
      )
      fit_i <- lme4::glmer(
        stats::as.formula(paste(outcome, "~", rhs)),
        data = data, family = stats::binomial()
      )
      vari <- extract_glmer_stratum_variance(fit_i, group = stratum)
    }

    prev_var <- results[[i]]$Variance
    step_pcv <- if (prev_var > 0) (prev_var - vari) / results[[1]]$Variance else 0
    total_pcv <- if (results[[1]]$Variance > 0) {
      (results[[1]]$Variance - vari) / results[[1]]$Variance
    } else 0

    results[[i + 1]] <- data.frame(
      Step = i, Model = paste("Model", i),
      Added_Variable = vars[i],
      Variance = vari,
      Step_PCV = step_pcv,
      Total_PCV = total_pcv,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  class(out) <- c("stepwise_pcv", "data.frame")
  out
}

#' @export
print.stepwise_pcv <- function(x, ...) {
  cat("Stepwise PCV Decomposition\n")
  cat("==========================\n\n")
  print.data.frame(x, row.names = FALSE, digits = 4)
  cat("\nInterpretation:\n")
  cat("  Step_PCV: proportion of original between-stratum variance\n")
  cat("            explained by each added variable.\n")
  cat("  Total_PCV: cumulative proportion explained so far.\n")
  cat("  Negative values indicate an unmasking/suppression pattern.\n")
  invisible(x)
}
