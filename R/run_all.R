# run_all.R — Full I-MAIHDA synthetic simulation pipeline (R version)
#
# Replicates scripts/run_all.py in R.
# Output: R_outputs/*.csv, R_figures/*.png (publication quality), R_logs/run_log.md
#
# Usage: Rscript R/run_all.R
#   (run from the imaihda_hic_mic_v31/ directory)

# ── Setup ─────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(jsonlite)
  library(scales)
})

# Resolve paths relative to project root (parent of R/)
# Works both interactively (sys.frame) and via Rscript (commandArgs)
find_root <- function() {
  # When run via Rscript, the --file= argument gives the script path
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(dirname(script_path)), winslash = "/"))
  }
  # Interactive: try sys.frame
  tryCatch({
    p <- normalizePath(dirname(dirname(sys.frame(1)$ofile)), winslash = "/")
    if (!is.na(p) && nchar(p) > 0) return(p)
  }, error = function(e) NULL)
  # Fallback: assume working directory is the project root
  normalizePath(getwd(), winslash = "/")
}
ROOT <- find_root()
OUT  <- file.path(ROOT, "R_outputs")
FIG  <- file.path(ROOT, "R_figures")
LOG  <- file.path(ROOT, "R_logs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG, showWarnings = FALSE, recursive = TRUE)

# Source R modules
source(file.path(ROOT, "R", "simulate.R"))
source(file.path(ROOT, "R", "fit.R"))
source(file.path(ROOT, "R", "benchmarks.R"))

# ── Detection sweep ───────────────────────────────────────────────────────────
detection_sweep <- function(seed = 42L) {
  strengths <- c(0.0, 0.3, 0.6, 0.9, 1.2)
  rows <- lapply(strengths, function(s) {
    df <- simulate_intersectional_data(interaction_sd = 0.90, detection_strength = s, seed = seed)
    res <- fit_imaihda(df)
    res$detection_strength <- s
    res
  })
  # Bind rows, filling missing columns with NA
  all_names <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    for (nm in setdiff(all_names, names(r))) r[[nm]] <- NA
    as.data.frame(r, stringsAsFactors = FALSE)
  }))
}

# ── Figures (publication-quality ggplot2) ─────────────────────────────────────
plot_results <- function(results, sweep_df) {
  # Figure 1: VPC-PCV scenario scatter
  p1 <- ggplot(results, aes(x = vpc_null, y = pcv)) +
    geom_point(size = 3.5, color = "#21918c") +
    geom_text(aes(label = scenario), hjust = -0.3, vjust = 0.3, size = 4,
              family = "serif") +
    labs(
      x        = "Null-model VPC (%)",
      y        = "PCV after additive main effects (%)",
      title    = "Synthetic I-MAIHDA scenarios: VPC/PCV behaviour",
      subtitle = "R reproduction — Mersenne Twister RNG"
    ) +
    theme_bw(base_size = 13, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 9, color = "gray40")
    )
  ggsave(file.path(FIG, "scenario_vpc_pcv.png"), p1, width = 8, height = 5,
         dpi = 300)

  # Figure 2: Detection bias sweep
  sweep_df$obs_prev_pct <- 100 * sweep_df$overall_prevalence

  p2 <- ggplot(sweep_df, aes(x = detection_strength)) +
    geom_line(aes(y = vpc_null, color = "Null-model VPC"), linewidth = 1.0) +
    geom_point(aes(y = vpc_null, color = "Null-model VPC"), size = 2.5) +
    geom_line(aes(y = obs_prev_pct / 5, color = "Observed prevalence (%)"),
              linewidth = 1.0, linetype = "dashed") +
    geom_point(aes(y = obs_prev_pct / 5, color = "Observed prevalence (%)"),
               size = 2.5, shape = 17) +
    scale_y_continuous(
      name   = "Null-model VPC (%)",
      sec.axis = sec_axis(~ . * 5, name = "Observed prevalence (%)")
    ) +
    scale_color_manual(
      name   = "",
      values = c("Null-model VPC" = "#440154", "Observed prevalence (%)" = "#21918c")
    ) +
    labs(
      x        = "SES-patterned under-detection strength",
      title    = "Detection bias can mask residual intersectional heterogeneity",
      subtitle = "Interaction SD = 0.90; prevalence declines with detection strength"
    ) +
    theme_bw(base_size = 13, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 9, color = "gray40"),
      legend.position  = "bottom"
    )
  ggsave(file.path(FIG, "detection_sweep.png"), p2, width = 8, height = 5,
         dpi = 300)

  invisible(list(p1 = p1, p2 = p2))
}

# ── Main ──────────────────────────────────────────────────────────────────────
main <- function() {
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  # Run all scenarios
  grid <- scenario_grid()
  rows <- lapply(names(grid), function(nm) fit_scenario(nm, grid[[nm]]))

  # Collect results
  results <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  col_order <- c(
    "scenario", "description", "n", "n_strata", "min_stratum_n",
    "median_stratum_n", "max_stratum_n", "overall_prevalence",
    "true_prevalence", "mean_detection_probability",
    "var_null", "vpc_null", "var_main", "vpc_main", "pcv",
    "interaction_sd", "detection_strength", "sparse", "warnings"
  )
  results <- results[, intersect(col_order, names(results)), drop = FALSE]

  # Detection sweep
  sweep <- detection_sweep()

  # Benchmarks
  benchmarks <- evaluate_benchmarks(results)

  # Write CSVs
  write.csv(results,    file.path(OUT, "results.csv"),          row.names = FALSE)
  write.csv(sweep,      file.path(OUT, "detection_sweep.csv"),  row.names = FALSE)
  write.csv(benchmarks, file.path(OUT, "benchmark_checks.csv"), row.names = FALSE)

  # Figures
  plot_results(results, sweep)

  # Metadata
  finished <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  metadata <- list(
    started_utc           = started,
    finished_utc          = finished,
    r_version             = as.character(R.version.string),
    platform              = R.version$platform,
    all_benchmarks_passed = all(as.logical(benchmarks$passed))
  )
  writeLines(toJSON(metadata, pretty = TRUE, auto_unbox = TRUE),
             file.path(OUT, "run_metadata.json"))

  # Log
  log_lines <- c(
    "# R Run Log",
    paste("Started: ", started),
    paste("Finished:", finished),
    paste("R:       ", R.version.string),
    "",
    "## Benchmark checks",
    capture.output(print(benchmarks, row.names = FALSE)),
    "",
    "## Scenario results",
    capture.output(print(results[, c("scenario", "overall_prevalence", "vpc_null",
                                      "vpc_main", "pcv", "min_stratum_n")],
                          row.names = FALSE)),
    "",
    "## Note on cross-language comparison",
    "R uses Mersenne Twister (set.seed); Python uses PCG64 (numpy.random.default_rng).",
    "Exact numerical values differ, but statistical patterns and benchmark pass/fail",
    "should be equivalent. If benchmark pass/fail disagrees, investigate seed sensitivity."
  )
  writeLines(log_lines, file.path(LOG, "run_log.md"))

  # Checkpoint
  cat(sprintf("Final run completed. All benchmarks passed: %s\n",
              all(as.logical(benchmarks$passed))),
      file = file.path(ROOT, "R_checkpoint_final_run.txt"))

  # Report
  cat("\n═══════════════════════════════════════════\n")
  cat(" R I-MAIHDA simulation — run complete\n")
  cat("═══════════════════════════════════════════\n\n")
  cat("Benchmark checks:\n")
  print(benchmarks, row.names = FALSE)
  cat("\nScenario results:\n")
  print(results[, c("scenario", "overall_prevalence", "vpc_null", "vpc_main", "pcv", "min_stratum_n")],
        row.names = FALSE)
  cat(sprintf("\nAll benchmarks passed: %s\n", all(as.logical(benchmarks$passed))))
  cat(sprintf("Output: %s\n", OUT))
  cat(sprintf("Figures: %s\n", FIG))
  cat(sprintf("Log: %s\n", LOG))

  if (!all(as.logical(benchmarks$passed))) {
    stop("One or more benchmark checks failed; inspect R_outputs/benchmark_checks.csv")
  }

  invisible(list(results = results, sweep = sweep, benchmarks = benchmarks))
}

main()
