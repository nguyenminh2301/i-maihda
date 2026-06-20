# utils.R — Internal helper functions (not exported)

# Logistic function
logit_inv <- function(x) {
  1 / (1 + exp(-x))
}

# Logit (log-odds) transform
logit <- function(p) {
  log(p / (1 - p))
}

# Logistic level-1 variance (variance of standard logistic distribution)
LOGISTIC_L1_VARIANCE <- pi^2 / 3
