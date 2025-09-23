#' Convert t-Statistic to Pearson Correlation (r-Score)
#'
#' This function converts a t-statistic into a Pearson correlation coefficient (r-score),
#' given the sample size. The conversion is based on the relationship between the t-statistic
#' and the correlation coefficient in linear regression.
#'
#' @param t_statistic A numeric value representing the t-statistic.
#' @param sample_size An integer representing the total sample size.
#'
#' @return A numeric value representing the Pearson correlation coefficient (r-score).
#'
#' @details
#' The conversion formula used is:
#' \deqn{ r = \frac{t}{\sqrt{(N - 2) + t^2}} }
#' where:
#' \itemize{
#'   \item  r  is the Pearson correlation coefficient.
#'   \item  t  is the t-statistic.
#'   \item  N  is the total sample size.
#' }
#'
#' This formula follows from the equivalence of t-statistics and correlation coefficients
#' in simple linear regression.
#'
#' @note The function requires that `sample_size > 2`, since degrees of freedom must be positive.
#'
#' @examples
#' # Example usage:
#' tstat2rscore(t_statistic = 2.5, sample_size = 30)
#'
#' @export
tstat2rscore <- function(t_statistic, sample_size) {
  if (sample_size <= 2) {
    stop("Sample size must be greater than 2.")
  }
  r <- t_statistic / sqrt(sample_size - 2 + t_statistic^2)
  return(r)
}
