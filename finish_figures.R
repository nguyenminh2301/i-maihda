# finish_figures.R — Generate remaining figures 9-10
devtools::load_all("imaihda")
library(ggplot2)
library(viridis)
OUT <- "figures"
theme_pub <- theme_bw(base_size = 13, base_family = "serif") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        plot.caption = element_text(size = 8, color = "gray60", hjust = 0),
        legend.position = "bottom")

SEED <- 42L

cat("Generating sweep data...\n")
strengths <- seq(0, 1.2, by = 0.2)
sweep_list <- lapply(strengths, function(ds) {
  df <- simulate_intersectional_data(
    n = 6000, interaction_sd = 0.90, detection_strength = ds, seed = SEED
  )
  res <- fit_imaihda(df, method = "fast")
  data.frame(
    detection_strength = ds,
    vpc_null = res$vpc_null,
    overall_prevalence = res$overall_prevalence,
    stringsAsFactors = FALSE
  )
})
sweep_df <- do.call(rbind, sweep_list)
cat(sprintf("Sweep data: %d rows, cols: %s\n", nrow(sweep_df), paste(names(sweep_df), collapse=", ")))

cat("Figure 9: plot_sweep() example...\n")
png(file.path(OUT, "plot_sweep_example.png"), width = 9, height = 5.5, units = "in", res = 300)
plot_sweep(sweep_df, title = "plot_sweep() — Detection Bias Diagnostic (Unique to imaihda)")
dev.off()

cat("Figure 10: stepwise_pcv() bar chart...\n")
df_ex <- simulate_intersectional_data(n = 5000, interaction_sd = 0.9, seed = 123)
sw <- stepwise_pcv(df_ex, outcome = "y", vars = c("sex", "education", "wealth", "rural"),
                   method = "fast", quiet = TRUE)

sw_plot <- sw
sw_plot$Step_PCV_pct <- sw_plot$Step_PCV * 100
sw_plot$Total_PCV_pct <- sw_plot$Total_PCV * 100

p10 <- ggplot(sw_plot[-1, ], aes(x = factor(Step), y = Step_PCV_pct, fill = Added_Variable)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%+.1f%%", Step_PCV_pct)),
            vjust = ifelse(sw_plot$Step_PCV_pct[-1] >= 0, -0.3, 1.3),
            family = "serif", size = 3.5) +
  scale_fill_viridis_d(option = "D", name = "Variable") +
  labs(
    title = "Stepwise PCV Decomposition",
    subtitle = sprintf("Stepwise contribution of each social dimension. Total PCV = %.1f%%.",
                       sw_plot$Total_PCV_pct[nrow(sw_plot)]),
    x = "Step", y = "Step PCV (%)",
    caption = "imaihda v0.2.1 | stepwise_pcv(method='fast')"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_pub + theme(legend.position = "right")
ggsave(file.path(OUT, "stepwise_pcv_example.png"), p10, width = 9, height = 5.5, dpi = 300)

cat("Done.\n")
