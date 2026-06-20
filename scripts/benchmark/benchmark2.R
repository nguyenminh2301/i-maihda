# benchmark2.R — Clean benchmarking: imaihda fast vs glmer vs CRAN MAIHDA
library(imaihda)
library(ggplot2)
has_maihda <- requireNamespace("MAIHDA", quietly = TRUE)

OUT_DIR <- "benchmark_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Benchmark one configuration ──────────────────────────────────────────────
bench_one <- function(n, run_seed, use_method) {
  set.seed(run_seed)
  df <- simulate_intersectional_data(
    n = n, interaction_sd = 0.9, seed = run_seed
  )
  
  gc(reset = TRUE)
  t <- system.time({
    if (use_method == "fast") {
      res <- fit_imaihda(df, method = "fast")
    } else if (use_method == "glmer") {
      res <- fit_imaihda(df, method = "glmer")
    } else if (use_method == "maihda") {
      res_m <- MAIHDA::maihda(
        y ~ sex + education + wealth + rural + (1 | stratum),
        data = df, family = "binomial"
      )
      res <- list(
        var_null = res_m[["pcv"]][["var_model1"]],
        var_main = res_m[["pcv"]][["var_model2"]],
        vpc_null = res_m[["summary"]][["vpc"]][["estimate"]] * 100,
        pcv      = res_m[["pcv"]][["pvc"]] * 100
      )
    }
  })
  mem <- gc(reset = FALSE)
  mem_mb <- sum(mem[, 2], na.rm = TRUE)
  
  c(
    res[c("var_null", "vpc_null", "pcv")],
    list(time = as.numeric(t["elapsed"]), mem_mb = mem_mb, n = n)
  )
}

# ── Run ──────────────────────────────────────────────────────────────────────
SEED_BASE <- 42L
sizes <- c(10000L, 50000L, 100000L, 500000L, 1000000L, 2000000L)
sizes_glmer <- c(10000L, 50000L, 100000L)

cat("===== BENCHMARK START =====\n")
rows <- list()

for (n in sizes) {
  cat(sprintf("\n--- n = %s ---\n", format(n, big.mark = ",")))
  
  # Fast
  for (r in 1:3) {
    cat(sprintf("  fast run %d...", r))
    b <- bench_one(n, SEED_BASE + r, "fast")
    rows[[length(rows) + 1]] <- data.frame(
      n = n, method = "imaihda-fast", run = r,
      time_sec = b$time, mem_mb = b$mem_mb,
      vpc_null = b$vpc_null, pcv = b$pcv, var_null = b$var_null,
      stringsAsFactors = FALSE
    )
    cat(sprintf(" %.2fs\n", b$time))
  }
  
  # Glmer (smaller sizes only)
  if (n %in% sizes_glmer) {
    for (r in 1:2) {
      cat(sprintf("  glmer run %d...", r))
      b <- tryCatch(
        bench_one(n, SEED_BASE + r, "glmer"),
        error = function(e) { cat(sprintf(" ERROR: %s\n", e$message)); NULL }
      )
      if (!is.null(b)) {
        rows[[length(rows) + 1]] <- data.frame(
          n = n, method = "imaihda-glmer", run = r,
          time_sec = b$time, mem_mb = b$mem_mb,
          vpc_null = b$vpc_null, pcv = b$pcv, var_null = b$var_null,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # CRAN MAIHDA
  if (n %in% sizes_glmer && has_maihda) {
    for (r in 1:2) {
      cat(sprintf("  MAIHDA run %d...", r))
      b <- tryCatch(
        bench_one(n, SEED_BASE + r, "maihda"),
        error = function(e) { cat(sprintf(" ERROR: %s\n", e$message)); NULL }
      )
      if (!is.null(b)) {
        rows[[length(rows) + 1]] <- data.frame(
          n = n, method = "CRAN-MAIHDA", run = r,
          time_sec = b$time, mem_mb = b$mem_mb,
          vpc_null = b$vpc_null, pcv = b$pcv, var_null = b$var_null,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

bench_df <- do.call(rbind, rows)
write.csv(bench_df, file.path(OUT_DIR, "benchmark_results.csv"), row.names = FALSE)

# Aggregate
agg <- aggregate(
  cbind(time_sec, vpc_null, pcv) ~ n + method, data = bench_df, FUN = mean
)
agg_sd <- aggregate(
  time_sec ~ n + method, data = bench_df, FUN = sd
)
names(agg_sd)[3] <- "time_sd"
agg <- merge(agg, agg_sd, by = c("n", "method"), all.x = TRUE)

cat("\n===== AGGREGATED RESULTS =====\n")
print(agg[order(agg$n, agg$method), ], row.names = FALSE)

# ── Plots ────────────────────────────────────────────────────────────────────
if (nrow(bench_df) > 0) {
  # Time scaling
  p <- ggplot(agg, aes(x = n, y = time_sec, color = method)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c(
      "imaihda-fast" = "#21918c", "imaihda-glmer" = "#440154",
      "CRAN-MAIHDA" = "#fde725"
    )) +
    labs(title = "Computation Time by Method and Sample Size",
         subtitle = "Log-log scale. Fast method scales linearly; GLMM super-linearly.",
         x = "Sample size", y = "Time (seconds)", color = "Method") +
    theme_bw(13, "serif") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(OUT_DIR, "benchmark_time.png"), p, width = 8, height = 5, dpi = 300)

  # VPC convergence
  p2 <- ggplot(agg, aes(x = factor(n), y = vpc_null, fill = method)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(
      "imaihda-fast" = "#21918c", "imaihda-glmer" = "#440154",
      "CRAN-MAIHDA" = "#fde725"
    )) +
    labs(title = "VPC Convergence: Fast vs GLMM by Sample Size",
         subtitle = "Fast method converges to GLMM estimate as n grows",
         x = "Sample size", y = "VPC null (%)", fill = "Method") +
    theme_bw(13, "serif") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(OUT_DIR, "benchmark_vpc.png"), p2, width = 8, height = 5, dpi = 300)
  
  cat(sprintf("\nPlots: %s/benchmark_*.png\n", OUT_DIR))
}

cat("\nDone.\n")
