# Quick VPC check
library(imaihda)
library(lme4)

set.seed(42)
df <- simulate_intersectional_data(n = 5000, interaction_sd = 0.9, seed = 42)
cat("Overall prevalence:", mean(df$y), "\n\n")

# Fast
r_fast <- fit_imaihda(df, method = "fast")
cat("FAST vpc_null:", r_fast$vpc_null, "%  var_null:", r_fast$var_null, "\n")
cat("FAST vpc_main:", r_fast$vpc, "%  var_main:", r_fast$var_main, "\n\n")

# Glmer
r_glm <- fit_imaihda(df, method = "glmer")
cat("GLMER vpc_null:", r_glm$vpc_null, "%  var_null:", r_glm$var_null, "\n")
cat("GLMER vpc_main:", r_glm$vpc, "%  var_main:", r_glm$var_main, "\n\n")

# Manual lme4 check
null_fit <- glmer(y ~ (1 | stratum), data = df, family = binomial())
null_var <- as.numeric(VarCorr(null_fit)$stratum)
cat("Manual null var:", null_var,
    "  VPC:", round(100 * null_var / (null_var + pi^2/3), 2), "%\n")

main_fit <- glmer(y ~ factor(sex) + factor(education) + factor(wealth) +
                  factor(rural) + (1 | stratum),
                  data = df, family = binomial())
main_var <- as.numeric(VarCorr(main_fit)$stratum)
cat("Manual main var:", main_var,
    "  VPC:", round(100 * main_var / (main_var + pi^2/3), 2), "%\n")

# Check fast method intermediate
strata_agg <- aggregate_strata(df)
cat("\n--- Fast method details ---\n")
cat("N strata:", length(unique(strata_agg$empirical_logit)), "\n")
cat("Between-stratum var (unadjusted):",
    between_stratum_variance(
      strata_agg$empirical_logit,
      rep(mean(strata_agg$empirical_logit), nrow(strata_agg)),
      strata_agg$logit_sampling_var,
      strata_agg$precision_weight
    ), "\n")
