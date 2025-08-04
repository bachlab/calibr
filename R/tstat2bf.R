#' Compute Bayes Factor from a t-Statistic
#'
#' This function calculates the Bayes Factor (BF) based on a given t-statistic
#' and sample size, using the g-prior framework. The Bayes Factor quantifies
#' the strength of evidence against a ground measure (a constant). This is the
#' equivalent of null hypothesis testing.
#'
#' @param t_statistic A numeric value representing the t-statistic.
#' @param sample_size A numeric value representing the sample size.
#' @param g A numeric value representing the prior scaling parameter.
#'          Defaults to `g = sample_size`.
#'
#' @return A numeric value representing the computed Bayes Factor.
#' @export
#'
#' @examples
#' # Example: Compute BF for t-statistic of 2.5 with a sample size of 30
#' tstat2bf(2.5, 30)
#'
#' # Example: Compute BF with a custom g parameter
#' tstat2bf(2.5, 30, g = 15)
tstat2bf <- function(t_statistic, sample_size, g = sample_size) {

  # Check if all inputs are scalars (single numeric values)
  if (!is.numeric(t_statistic) || length(t_statistic) != 1) {
    stop("Error: 't_statistic' must be a single numeric value.")
  }

  if (!is.numeric(sample_size) || length(sample_size) != 1) {
    stop("Error: 'sample_size' must be a single numeric value.")
  }

  if (!is.numeric(g) || length(g) != 1) {
    stop("Error: 'g' must be a single numeric value.")
  }

  # Compute BF using the g-prior relationship
  bayes_factor <- exp(((sample_size - 2) * log(1 + g) - (sample_size - 1) * log(1 + g/(t_statistic^2/(sample_size - 2) + 1)))/2)

  return(bayes_factor)
}
