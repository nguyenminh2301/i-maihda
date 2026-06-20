# regenerate_stepwise_fig.R
devtools::load_all("imaihda")
library(ggplot2)
library(viridis)

theme_pub <- theme_bw(base_size = 13, base_family = "serif") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        plot.caption = element_text(size = 8, color = "gray60", hjust = 0),
        legend.position = "bottom")

set.seed(123)
df <- simulate_intersectional_data(n = 5000, interaction_sd = 0.9, seed = 123)
sw <- stepwise_pcv(df, outcome = "y", vars = c("sex", "education", "wealth", "rural"),
                   method = "fast", quiet = TRUE)

# Build stepwise plot data — exclude step 0 (null model)
sw_plot <- sw[-1, ]
sw_plot$Step_PCV_pct <- sw_plot$Step_PCV * 100
sw_plot$Total_PCV_pct <- sw_plot$Total_PCV * 100

# Create stepwise bar chart with total PCV line
p <- ggplot(sw_plot, aes(x = factor(Step), y = Step_PCV_pct, fill = Added_Variable)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%+.1f%%", Step_PCV_pct)),
            vjust = ifelse(sw_plot$Step_PCV >= 0, -0.3, 1.3),
            family = "serif", size = 4) +
  scale_fill_viridis_d(option = "D", name = "Added Variable") +
  labs(
    title = "Stepwise PCV Decomposition",
    subtitle = sprintf("Cumulative PCV = %.1f%%. Each bar = proportion of original between-stratum variance explained by that variable.",
                       sw_plot$Total_PCV_pct[nrow(sw_plot)]),
    x = "Step",
    y = "Step PCV (%)",
    caption = "imaihda v0.2.1 | stepwise_pcv(method='fast') | Negative = suppression/unmasking"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_pub + theme(legend.position = "right")

ggsave("figures/stepwise_pcv_example.png", p, width = 9, height = 5.5, dpi = 300)
cat("Saved figures/stepwise_pcv_example.png\n")
