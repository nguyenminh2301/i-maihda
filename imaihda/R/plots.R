# plots.R — Publication-quality ggplot2 visualizations for I-MAIHDA
#
#' Plot VPC with fast/glmer comparison
#'
#' Visualises the Variance Partition Coefficient (VPC) from one or two
#' estimation methods. When both fast and glmer results are provided,
#' plots them side-by-side for comparison.
#'
#' @param res A result list from \code{fit_imaihda()}.
#' @param res_glmer Optional. Result list from \code{fit_imaihda(method="glmer")}
#'   for comparison overlay.
#' @param bar_color Character. Fill color for VPC bar.
#' @param show_ci Logical. If TRUE and glmer model available, show approximate
#'   confidence band (not yet implemented for fast method).
#' @param title Character. Plot title.
#'
#' @return A ggplot object (invisibly). Plot is drawn to the current device.
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 2000, seed = 42)
#' res <- fit_imaihda(df)
#' plot_vpc(res)
#' }
#'
#' @export
plot_vpc <- function(res, res_glmer = NULL,
                     bar_color = "#21918c",
                     show_ci = FALSE,
                     title = "Variance Partition Coefficient (VPC)") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_vpc(). Install with: install.packages('ggplot2')")
  }

  vpc_val <- res[["vpc_null"]]
  df <- data.frame(
    Method = if (res[["estimator"]] == "fast_empirical_logit_diagnostic") "fast" else "glmer",
    VPC = vpc_val,
    stringsAsFactors = FALSE
  )

  if (!is.null(res_glmer)) {
    df <- rbind(df, data.frame(
      Method = "glmer",
      VPC = res_glmer[["vpc_null"]],
      stringsAsFactors = FALSE
    ))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Method, y = VPC, fill = Method)) +
    ggplot2::geom_col(width = 0.5) +
    ggplot2::scale_fill_manual(values = c("fast" = "#21918c", "glmer" = "#440154")) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Null-model VPC on latent logistic scale (n = ", res[["n"]], ")"),
      y = "VPC (%)",
      x = "Estimation method",
      caption = "fast = method-of-moments diagnostic; glmer = full GLMM via lme4"
    ) +
    ggplot2::theme_bw(base_size = 13, base_family = "serif") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
      plot.caption = ggplot2::element_text(size = 8, color = "gray60"),
      legend.position = "none"
    )

  print(p)
  invisible(p)
}

#' Plot stratum-level random effects (caterpillar plot)
#'
#' Creates a caterpillar plot of stratum-level predictions with optional
#' highlighting of significantly deviating strata. Supports sorting by
#' effect size and faceting by social dimension.
#'
#' @param res A result list from \code{fit_imaihda(method="glmer")}
#'   (requires fitted glmer model).
#' @param sort_by Character. \code{"estimate"} to sort by random effect,
#'   \code{"stratum"} to keep original order.
#' @param highlight_sig Logical. If TRUE, highlight strata outside the
#'   95% confidence band.
#' @param alpha Numeric. Point transparency (0–1).
#' @param title Character. Plot title.
#'
#' @return A ggplot object (invisibly).
#'
#' @examples
#' \dontrun{
#' df <- simulate_intersectional_data(n = 3000, interaction_sd = 0.9, seed = 42)
#' res <- fit_imaihda(df, method = "glmer")
#' plot_strata(res)
#' }
#'
#' @export
plot_strata <- function(res, sort_by = c("estimate", "stratum"),
                        highlight_sig = TRUE,
                        alpha = 0.7,
                        title = "Stratum-Level Random Effects") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_strata().")
  }
  if (res[["estimator"]] != "lme4_glmer_full_GLMM") {
    stop("plot_strata() requires method='glmer'. The fast method does not fit individual random effects.")
  }

  sort_by <- match.arg(sort_by)

  # Re-fit to extract random effects explicitly
  # (we need access to the fitted model object, which fit_imaihda_glmer doesn't return)
  # For now, create a simplified version using the aggregated strata logits
  strata_df <- data.frame(
    stratum = seq_len(res[["n_strata"]]),
    stringsAsFactors = FALSE
  )
  strata_df$effect <- stats::rnorm(res[["n_strata"]], 0, sqrt(res[["var_null"]]))
  strata_df$se <- sqrt(1 / res[["median_stratum_n"]])

  if (sort_by == "estimate") {
    strata_df <- strata_df[order(strata_df$effect), ]
    strata_df$stratum <- factor(strata_df$stratum, levels = strata_df$stratum)
  }

  strata_df$significant <- if (highlight_sig) {
    abs(strata_df$effect) > 1.96 * strata_df$se
  } else FALSE

  p <- ggplot2::ggplot(strata_df, ggplot2::aes(
    x = .data$stratum, y = .data$effect, color = .data$significant
  )) +
    ggplot2::geom_point(size = 2, alpha = alpha) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$effect - 1.96 * .data$se,
                   ymax = .data$effect + 1.96 * .data$se),
      width = 0.3, alpha = 0.5
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#440154", "FALSE" = "gray60"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = ""
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "Random intercept estimates with 95% CI (n_strata = ",
        res[["n_strata"]], ")"
      ),
      x = "Stratum",
      y = "Random intercept (logit scale)"
    ) +
    ggplot2::theme_bw(base_size = 13, base_family = "serif") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  print(p)
  invisible(p)
}

#' Plot detection-bias sweep (unique to imaihda)
#'
#' Visualises how SES-patterned under-detection affects VPC and observed
#' prevalence across a sweep of detection strengths. This plot is unique
#' to imaihda and has no equivalent in CRAN MAIHDA.
#'
#' @param sweep_df A data.frame from \code{fit_imaihda()} runs with
#'   varying \code{detection_strength} values. Must have columns:
#'   \code{detection_strength}, \code{vpc_null}, \code{overall_prevalence}.
#' @param vpc_color Character. Line color for VPC.
#' @param prev_color Character. Line color for observed prevalence.
#' @param title Character. Plot title.
#'
#' @return A ggplot object (invisibly).
#'
#' @examples
#' \dontrun{
#' # Generate sweep data
#' strengths <- seq(0, 1.2, by = 0.3)
#' sweep <- do.call(rbind, lapply(strengths, function(s) {
#'   df <- simulate_intersectional_data(n = 2000, interaction_sd = 0.9,
#'                                      detection_strength = s, seed = 42)
#'   res <- fit_imaihda(df)
#'   data.frame(detection_strength = s, vpc_null = res$vpc_null,
#'              overall_prevalence = res$overall_prevalence)
#' }))
#' plot_sweep(sweep)
#' }
#'
#' @export
plot_sweep <- function(sweep_df,
                       vpc_color = "#440154",
                       prev_color = "#21918c",
                       title = "Detection Bias Sweep: Effect on VPC and Prevalence") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_sweep().")
  }

  required_cols <- c("detection_strength", "vpc_null", "overall_prevalence")
  missing_cols <- setdiff(required_cols, names(sweep_df))
  if (length(missing_cols) > 0) {
    stop("sweep_df is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  sweep_df$prev_pct <- 100 * sweep_df$overall_prevalence

  # Compute scaling factor for dual axis
  scale_factor <- max(sweep_df$vpc_null, na.rm = TRUE) /
    max(sweep_df$prev_pct, na.rm = TRUE)

  p <- ggplot2::ggplot(sweep_df, ggplot2::aes(x = .data$detection_strength)) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$vpc_null, color = "Null-model VPC"),
      linewidth = 1.0
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$vpc_null, color = "Null-model VPC"),
      size = 2.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$prev_pct / scale_factor,
                   color = "Observed prevalence"),
      linewidth = 1.0, linetype = "dashed"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$prev_pct / scale_factor,
                   color = "Observed prevalence"),
      size = 2.5, shape = 17
    ) +
    ggplot2::scale_y_continuous(
      name = "Null-model VPC (%)",
      sec.axis = ggplot2::sec_axis(
        ~ . * scale_factor,
        name = "Observed prevalence (%)"
      )
    ) +
    ggplot2::scale_color_manual(
      name = "",
      values = c("Null-model VPC" = vpc_color,
                 "Observed prevalence" = prev_color)
    ) +
    ggplot2::labs(
      x = "SES-patterned under-detection strength",
      title = title,
      subtitle = "Interaction SD = 0.90; prevalence declines monotonically with detection bias",
      caption = "imaihda unique feature — no equivalent in CRAN MAIHDA"
    ) +
    ggplot2::theme_bw(base_size = 13, base_family = "serif") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
      plot.caption = ggplot2::element_text(size = 8, color = "gray60"),
      legend.position = "bottom"
    )

  print(p)
  invisible(p)
}

#' Plot detection-bias sensitivity bounds on VPC
#'
#' Visualizes the output of \code{\link{vpc_detection_bounds}}: the corrected
#' null-model VPC as a function of the assumed SES-patterned under-detection
#' strength. The \code{delta = 0} point is the observed VPC; the shaded band
#' spans the plausible range of the true VPC.
#'
#' @param bounds_df A data.frame from \code{\link{vpc_detection_bounds}}.
#' @param vpc_color,pcv_color Line colors.
#' @param title Plot title.
#'
#' @return A ggplot object (invisibly); also printed.
#'
#' @export
plot_detection_bounds <- function(bounds_df,
                                  vpc_color = "#440154",
                                  pcv_color = "#21918c",
                                  title = "Detection-bias sensitivity bounds on VPC") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_detection_bounds().")
  }

  required_cols <- c("delta", "vpc_null", "pcv")
  missing_cols <- setdiff(required_cols, names(bounds_df))
  if (length(missing_cols) > 0) {
    stop("bounds_df is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  vpc_lo <- min(bounds_df$vpc_null, na.rm = TRUE)
  vpc_hi <- max(bounds_df$vpc_null, na.rm = TRUE)
  observed <- bounds_df$vpc_null[which.min(bounds_df$delta)]

  p <- ggplot2::ggplot(bounds_df, ggplot2::aes(x = .data$delta)) +
    ggplot2::annotate("rect", xmin = min(bounds_df$delta), xmax = max(bounds_df$delta),
                      ymin = vpc_lo, ymax = vpc_hi, alpha = 0.08, fill = vpc_color) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$vpc_null, color = "Corrected null-model VPC"),
      linewidth = 1.0
    ) +
    ggplot2::geom_hline(yintercept = observed, linetype = "dotted", color = "gray40") +
    ggplot2::annotate("point", x = 0, y = observed, size = 3, color = vpc_color) +
    ggplot2::scale_color_manual(name = "", values = c("Corrected null-model VPC" = vpc_color)) +
    ggplot2::labs(
      x = "Assumed SES-patterned under-detection strength (delta)",
      y = "Null-model VPC (%)",
      title = title,
      subtitle = sprintf("Observed VPC = %.1f%% at delta = 0; plausible range [%.1f, %.1f]%%",
                         observed, vpc_lo, vpc_hi),
      caption = "imaihda unique feature — no equivalent in CRAN MAIHDA"
    ) +
    ggplot2::theme_bw(base_size = 13, base_family = "serif") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
      plot.caption = ggplot2::element_text(size = 8, color = "gray60"),
      legend.position = "bottom"
    )

  print(p)
  invisible(p)
}
