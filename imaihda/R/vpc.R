# vpc.R — Core MAIHDA formulas: VPC and PCV
#
# These are the two fundamental metrics used in Intersectional MAIHDA
# to quantify the proportion of variance attributable to intersectional
# strata and the reduction achieved by additive main effects.

#' Variance Partition Coefficient on the Latent Logistic Scale
#'
#' Computes the VPC for a binary outcome in a multilevel logistic model
#' using the latent-response approach. The level-1 (individual) variance
#' is fixed at \eqn{\pi^2/3}, the variance of the standard logistic
#' distribution. This is the standard approach in the I-MAIHDA literature
#' (Evans et al. 2018; O'Sullivan et al. 2024).
#'
#' VPC represents the proportion of total variance in the latent outcome
#' that lies between intersectional strata. Higher VPC indicates greater
#' between-stratum inequality after accounting for individual-level
#' Bernoulli variance.
#'
#' @param stratum_variance Non-negative numeric scalar. Estimated
#'   between-stratum variance (usually on the logit scale).
#'
#' @return Numeric scalar: VPC as a percentage (0--100).
#'
#' @references
#' Evans CR, Williams DR, Onnela JP, Subramanian SV (2018).
#' "A multilevel approach to modeling health inequalities at the
#' intersection of multiple social identities."
#' \emph{SSM - Population Health}, 6, 149--157.
#'
#' O'Sullivan JL, Alonso-Perez E, et al. (2024).
#' "Onset of Type 2 diabetes in adults aged 50 and older in Europe:
#' an intersectional multilevel analysis of individual heterogeneity
#' and discriminatory accuracy."
#' \emph{Diabetology & Metabolic Syndrome}, 16:293.
#'
#' @examples
#' vpc_latent(0.5)   # ~13.2%
#' vpc_latent(0)     # 0%
#' vpc_latent(3.29)  # 50% (stratum variance equals pi^2/3)
#'
#' @export
vpc_latent <- function(stratum_variance) {
  if (!is.numeric(stratum_variance) || length(stratum_variance) != 1) {
    stop("stratum_variance must be a single numeric value")
  }
  if (stratum_variance < 0) {
    stop("stratum_variance must be non-negative")
  }
  100 * stratum_variance / (stratum_variance + LOGISTIC_L1_VARIANCE)
}

#' Proportional Change in Variance
#'
#' Computes the proportional reduction in between-stratum variance when
#' moving from the null (intersectional-only) model to the main-effects
#' (additive) model. This is the key metric for assessing how much of
#' the intersectional inequality can be explained by additive social
#' determinants.
#'
#' \deqn{PCV = \frac{\sigma^2_{null} - \sigma^2_{main}}{\sigma^2_{null}} \times 100\%}
#'
#' A PCV near 100\% suggests the intersectional structure is largely
#' additive (i.e., social determinants act independently). A low PCV
#' suggests residual intersectional heterogeneity not reducible to
#' additive main effects.
#'
#' @param var_null Numeric scalar. Between-stratum variance from the
#'   null model (intersectional strata only, no covariates).
#' @param var_main Numeric scalar. Between-stratum variance from the
#'   main-effects model (strata + additive social determinants).
#'
#' @return Numeric scalar: PCV as a percentage. Returns \code{NaN} if
#'   \code{var_null} is zero or negative (undefined).
#'
#' @references
#' Evans CR, Williams DR, Onnela JP, Subramanian SV (2018).
#' "A multilevel approach to modeling health inequalities at the
#' intersection of multiple social identities."
#' \emph{SSM - Population Health}, 6, 149--157.
#'
#' @examples
#' pcv(1.0, 0.25)  # 75% — most variance explained by additive effects
#' pcv(0.5, 0.4)   # 20% — substantial residual interaction
#' pcv(0.0, 0.0)   # NaN — undefined when null variance is zero
#'
#' @export
pcv <- function(var_null, var_main) {
  if (!is.numeric(var_null) || length(var_null) != 1) {
    stop("var_null must be a single numeric value")
  }
  if (!is.numeric(var_main) || length(var_main) != 1) {
    stop("var_main must be a single numeric value")
  }
  if (var_null <= 0) {
    return(NaN)
  }
  100 * (var_null - var_main) / var_null
}
