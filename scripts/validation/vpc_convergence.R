# VPC convergence analysis
library(imaihda)
library(lme4)

do_check <- function(n, seed) {
  df <- simulate_intersectional_data(n = n, interaction_sd = 0.9, seed = seed)
  
  # Fast: null model
  strata_agg <- aggregate_strata(df)
  fast_var_null <- between_stratum_variance(
    strata_agg$empirical_logit,
    rep(mean(strata_agg$empirical_logit), nrow(strata_agg)),
    strata_agg$logit_sampling_var,
    strata_agg$precision_weight
  )
  fast_vpc_null <- 100 * fast_var_null / (fast_var_null + LOGISTIC_L1_VARIANCE)
  
  # Fast: main model (residual after controlling for additive effects)
  pred_main_glm <- glm(y ~ factor(sex) + factor(education) + factor(wealth) + factor(rural),
                       data = df, family = binomial())
  strata_agg_main <- aggregate_strata(df)
  pred_main_link <- predict(pred_main_glm, newdata = strata_agg_main, type = "link")
  fast_var_main <- between_stratum_variance(
    strata_agg_main$empirical_logit, pred_main_link,
    strata_agg_main$logit_sampling_var, strata_agg_main$precision_weight
  )
  fast_vpc_main <- 100 * fast_var_main / (fast_var_main + LOGISTIC_L1_VARIANCE)
  
  # Glmer
  null_fit <- glmer(y ~ (1 | stratum), data = df, family = binomial())
  glmer_var_null <- as.numeric(VarCorr(null_fit)$stratum)
  glmer_vpc_null <- 100 * glmer_var_null / (glmer_var_null + pi^2/3)
  
  main_fit <- glmer(y ~ factor(sex) + factor(education) + factor(wealth) + factor(rural) + (1|stratum),
                    data = df, family = binomial())
  glmer_var_main <- as.numeric(VarCorr(main_fit)$stratum)
  glmer_vpc_main <- 100 * glmer_var_main / (glmer_var_main + pi^2/3)
  
  c(n = n, 
    fast_var_null = unname(fast_var_null), fast_vpc_null = unname(fast_vpc_null),
    fast_var_main = unname(fast_var_main), fast_vpc_main = unname(fast_vpc_main),
    glmer_var_null = unname(glmer_var_null), glmer_vpc_null = unname(glmer_vpc_null),
    glmer_var_main = unname(glmer_var_main), glmer_vpc_main = unname(glmer_vpc_main),
    prevalence = unname(mean(df$y)))
}

seeds <- 42:44
sizes <- c(5000L, 10000L, 20000L, 50000L, 100000L, 200000L)

cat("=== VPC Convergence Analysis ===\n\n")
cat(sprintf("%-10s %6s %7s %7s %7s %7s\n",
            "n", "prev", "fVPC0", "gVPC0", "fVPC1", "gVPC1"))
cat(strrep("-", 50), "\n")

for (n in sizes) {
  m <- matrix(NA, nrow = length(seeds), ncol = 10)
  colnames(m) <- c("n","prev","fast_var_null","fast_vpc_null","fast_var_main","fast_vpc_main",
                    "glmer_var_null","glmer_vpc_null","glmer_var_main","glmer_vpc_main")
  for (i in seq_along(seeds)) {
    res <- tryCatch(do_check(n, seeds[i]), error = function(e) NULL)
    if (!is.null(res) && length(res) == 10) {
      m[i, ] <- res
    }
  }
  cm <- colMeans(m, na.rm = TRUE)
  cat(sprintf("%-10s %6.4f %7.2f %7.2f %7.2f %7.2f\n",
              format(n, big.mark = ","),
              cm["prev"],
              cm["fast_vpc_null"], cm["glmer_vpc_null"],
              cm["fast_vpc_main"], cm["glmer_vpc_main"]))
}
