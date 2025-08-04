#' Compute t-Statistic from a Bayes Factor
#'
#' This function calculates the **t-statistic** corresponding to a given Bayes Factor (BF)
#' and sample size. The method follows the g-prior framework. The transformation used is
#' derived from eq.5 in Nikolakopoulos, S., & Ntzoufras, I. (2021). Meta Analysis of Bayes Factors.
#'
#' @param bayes_factor A numeric value representing the Bayes Factor.
#' @param sample_size A numeric value representing the sample size.
#' @param g A numeric value representing the prior scaling parameter.
#'          Defaults to `g = sample_size`, following Zellner's g-prior.
#'
#' @return A numeric value representing the computed **t-statistic**.
#' @export
#'
#' @examples
#' # Example: Compute t-statistic for a given Bayes Factor and sample size
#' bf2tstat(bayes_factor = 5, sample_size = 30)
#'
#' # Example: Compute t-statistic with a custom g parameter
#' bf2tstat(bayes_factor = 2.5, sample_size = 30, g = 15)
bf2tstat <- function(bayes_factor, sample_size, g = sample_size) {
  if (sample_size <= 2) stop("Sample size must be greater than 2.")
  if (bayes_factor <= 0) stop("Bayes factor must be > 0.")

  log_BF <- log(bayes_factor)
  exp_part <- exp((-2 * log_BF + (sample_size - 2) * log(1 + g)) / (sample_size - 1))
  denom <- exp_part - 1
  if (denom <= 0) {
    warning("Bayes factor is too low to map to a real t-statistic (denominator <= 0). Returning t = 0.")
    return(0)
  }

  inner <- (sample_size - 2) * (g / denom - 1)
  if (inner < 0 || !is.finite(inner)) {
    warning("Bayes factor is too low to map to a real t-statistic (inner < 0 or not finite). Returning t = 0.")
    return(0)
  }

  t_stat <- sqrt(inner)
  return(t_stat)
}


