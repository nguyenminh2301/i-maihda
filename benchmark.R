# benchmark.R — Comprehensive benchmarking: imaihda vs CRAN MAIHDA
#
# Measures computation time, RAM usage, and result accuracy across
# sample sizes from 10K to 2M. Compares fast vs glmer methods.

suppressPackageStartupMessages({
  library(imaihda)
  library(ggplot2)
  library(viridis)
})

has_maihda <- requireNamespace("MAIHDA", quietly = TRUE)

# ── Config ───────────────────────────────────────────────────────────────────
SIZES_FAST  <- c(10000L, 50000L, 100000L, 500000L, 1000000L, 2000000L)
SIZES_GLMER <- c(10000L, 50000L, 100000L)  # glmer too slow beyond 100K
N_RUNS <- 3L  # repeat runs for stable timing
SEED <- 42L

OUT_DIR <- "benchmark_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Timer helper ─────────────────────────────────────────────────────────────
time_run <- function(expr, gc_first = TRUE) {
  if (gc_first) gc(reset = TRUE)
  t <- system.time(expr, gcFirst = FALSE)
  # Rough memory estimate: sum of Ncells + Vcells (Mb)
  mem_info <- gc(reset = FALSE)
  mem_used <- sum(mem_info[, 2], na.rm = TRUE)  # column 2 = "used (Mb)"
  list(
    elapsed = as.numeric(t["elapsed"]),
    user    = as.numeric(t["user.self"]),
    sys     = as.numeric(t["sys.self"]),
    mem_mb  = mem_used
  )
}

# ── Single-run benchmark ─────────────────────────────────────────────────────
run_bench <- function(n, seed, method, interaction_sd = 0.9) {
  df <- simulate_intersectional_data(
    n = n, interaction_sd = interaction_sd, seed = seed
  )

  if (method == "fast") {
    t <- time_run(res <- fit_imaihda(df, method = "fast"))
    res$time <- t
    res
  } else if (method == "glmer") {
    t <- time_run(res <- fit_imaihda(df, method = "glmer"))
    res$time <- t
    res
  } else if (method == "maihda" && has_maihda) {
    t <- time_run({
      res_m <- MAIHDA::maihda(
        y ~ sex + education + wealth + rural + (1 | stratum),
        data = df, family = "binomial"
      )
    })
    list(
      var_null  = res_m[["pcv"]][["var_model1"]],
      var_main  = res_m[["pcv"]][["var_model2"]],
      vpc_null  = res_m[["summary"]][["vpc"]][["estimate"]] * 100,
      pcv       = res_m[["pcv"]][["pvc"]] * 100,
      n_strata  = res_m[["summary"]][["n_strata"]],
      overall_prevalence = mean(df$y),
      time = t
    )
  }
}

# ── Run benchmarks ───────────────────────────────────────────────────────────
cat("===== I-MAIHDA BENCHMARK: imaihda vs CRAN MAIHDA =====\n")
cat("Sample sizes (fast):", paste(SIZES_FAST, collapse = ", "), "\n")
cat("Sample sizes (glmer):", paste(SIZES_GLMER, collapse = ", "), "\n")
cat("Repeats:", N_RUNS, "\n\n")

results <- list()

for (n in SIZES_FAST) {
  cat(sprintf("\n--- n = %s ---\n", format(n, big.mark = ",")))

  # Fast method
  fast_times <- numeric(N_RUNS)
  fast_res <- NULL
  for (r in seq_len(N_RUNS)) {
    cat(sprintf("  fast run %d/%d...\n", r, N_RUNS))
    res <- run_bench(n, SEED + r, "fast")
    fast_times[r] <- res$time$elapsed
    if (r == 1) fast_res <- res
  }
  fast_avg <- mean(fast_times)

  results[[length(results) + 1]] <- data.frame(
    n = n, method = "imaihda-fast",
    time_sec = fast_avg,
    time_sd = if (N_RUNS > 1) sd(fast_times) else NA,
    vpc_null = fast_res$vpc_null,
    pcv = fast_res$pcv,
    var_null = fast_res$var_null,
    prevalence = fast_res$overall_prevalence,
    stringsAsFactors = FALSE
  )

  # Glmer method (only for smaller sizes)
  if (n %in% SIZES_GLMER) {
    glmer_times <- numeric(N_RUNS)
    glmer_res <- NULL
    for (r in seq_len(N_RUNS)) {
      cat(sprintf("  glmer run %d/%d...\n", r, N_RUNS))
      res <- tryCatch(
        run_bench(n, seed + r, "glmer"),
        error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
      )
      if (!is.null(res)) {
        glmer_times[r] <- res$time$elapsed
        if (r == 1) glmer_res <- res
      }
    }
    if (!is.null(glmer_res)) {
      glmer_avg <- mean(glmer_times[glmer_times > 0])
      results[[length(results) + 1]] <- data.frame(
        n = n, method = "imaihda-glmer",
        time_sec = glmer_avg,
        time_sd = if (sum(glmer_times > 0) > 1) sd(glmer_times[glmer_times > 0]) else NA,
        vpc_null = glmer_res$vpc_null,
        pcv = glmer_res$pcv,
        var_null = glmer_res$var_null,
        prevalence = glmer_res$overall_prevalence,
        stringsAsFactors = FALSE
      )
    }
  }

  # CRAN MAIHDA (only for smallest sizes — it's similar to glmer)
  if (n %in% SIZES_GLMER && has_maihda) {
    maihda_times <- numeric(N_RUNS)
    maihda_res <- NULL
    for (r in seq_len(N_RUNS)) {
      cat(sprintf("  MAIHDA run %d/%d...\n", r, N_RUNS))
      res <- tryCatch(
        run_bench(n, seed + r, "maihda"),
        error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
      )
      if (!is.null(res)) {
        maihda_times[r] <- res$time$elapsed
        if (r == 1) maihda_res <- res
      }
    }
    if (!is.null(maihda_res)) {
      maihda_avg <- mean(maihda_times[maihda_times > 0])
      results[[length(results) + 1]] <- data.frame(
        n = n, method = "CRAN-MAIHDA",
        time_sec = maihda_avg,
        time_sd = if (sum(maihda_times > 0) > 1) sd(maihda_times[maihda_times > 0]) else NA,
        vpc_null = maihda_res$vpc_null,
        pcv = maihda_res$pcv,
        var_null = maihda_res$var_null,
        prevalence = maihda_res$overall_prevalence,
        stringsAsFactors = FALSE
      )
    }
  }
}

# ── Save results ─────────────────────────────────────────────────────────────
bench_df <- do.call(rbind, results)
write.csv(bench_df, file.path(OUT_DIR, "benchmark_results.csv"), row.names = FALSE)

cat("\n===== RESULTS =====\n")
print(bench_df, row.names = FALSE)

# ── Plot: Time scaling ───────────────────────────────────────────────────────
if (nrow(bench_df) > 0) {
  bench_df$n_label <- paste0(format(bench_df$n / 1000, big.mark = ","), "K")

  p_time <- ggplot(bench_df, aes(x = n, y = time_sec, color = method, shape = method)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    scale_x_log10(
      breaks = unique(bench_df$n),
      labels = scales::comma
    ) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(
      values = c("imaihda-fast" = "#21918c", "imaihda-glmer" = "#440154",
                 "CRAN-MAIHDA" = "#fde725")
    ) +
    labs(
      title = "Computation Time by Method and Sample Size",
      subtitle = "Fast method-of-moments vs full GLMM (log-log scale)",
      x = "Sample size (n)",
      y = "Elapsed time (seconds)",
      color = "Method", shape = "Method"
    ) +
    theme_bw(base_size = 13, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "bottom"
    )
  ggsave(file.path(OUT_DIR, "benchmark_time.png"), p_time,
         width = 8, height = 5, dpi = 300)

  # ── Plot: VPC comparison ───────────────────────────────────────────────────
  p_vpc <- ggplot(
    bench_df[!is.na(bench_df$vpc_null), ],
    aes(x = factor(n), y = vpc_null, fill = method)
  ) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(
      values = c("imaihda-fast" = "#21918c", "imaihda-glmer" = "#440154",
                 "CRAN-MAIHDA" = "#fde725")
    ) +
    labs(
      title = "VPC Estimate by Method and Sample Size",
      subtitle = "Fast method converges to glmer estimate as n increases",
      x = "Sample size",
      y = "VPC null (%)",
      fill = "Method"
    ) +
    theme_bw(base_size = 13, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "bottom"
    )
  ggsave(file.path(OUT_DIR, "benchmark_vpc.png"), p_vpc,
         width = 8, height = 5, dpi = 300)

  cat(sprintf("\nPlots saved to %s/\n", OUT_DIR))
}

cat("\nBenchmark complete.\n")
