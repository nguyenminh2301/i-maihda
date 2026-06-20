library(imaihda); library(MAIHDA)
set.seed(42)
df <- simulate_intersectional_data(n = 500, seed = 42)
mdf <- MAIHDA::make_strata(df, c("sex", "education", "wealth", "rural"))
cat("Trying mdf$data...\n")
r1 <- tryCatch(
  MAIHDA::maihda(y ~ sex + education + wealth + rural + (1 | stratum),
                 data = mdf[["data"]], family = "binomial"),
  error = function(e) e[["message"]]
)
cat("Result:", if (is.character(r1)) r1 else "SUCCESS", "\n")

cat("\nTrying mdf directly...\n")
r2 <- tryCatch(
  MAIHDA::maihda(y ~ sex + education + wealth + rural + (1 | stratum),
                 data = mdf, family = "binomial"),
  error = function(e) e[["message"]]
)
cat("Result:", if (is.character(r2)) r2 else "SUCCESS", "\n")
