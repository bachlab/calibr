#' Compute t-statistic from a log Bayes Factor
#'
#' Inverts the g-prior Bayes factor relationship (Eq. 5) to recover the
#' corresponding \emph{t}-statistic magnitude. The log Bayes factor is
#' assumed to be in natural-log units, as returned by \code{tstat2logbf()}.
#'
#' @param log_BF Numeric, log Bayes factor (natural log).
#' @param sample_size Numeric, sample size (must be > 2).
#' @param g Numeric, prior scaling parameter (default = \code{sample_size}).
#'
#' @return Numeric t-statistic magnitude implied by the Bayes factor.
#'   Returns 0 if the mapping is not real (too small BF).
#'
#' @examples
#' t <- 2.5; n <- 30
#' lb <- tstat2logbf(t, n)
#' logbf2tstat(lb, n)
#'
#' @export
logbf2tstat <- function(log_BF, sample_size, g = sample_size) {
  if (sample_size <= 2) stop("Sample size must be greater than 2.")

  # stable denominator: expm1(x) = exp(x) - 1
  denom <- expm1(((sample_size - 2) * log1p(g) - 2 * log_BF) / (sample_size - 1))
  if (denom <= 0) return(0)

  inner <- (sample_size - 2) * (g / denom - 1)
  if (inner <= 0) return(0)

  sqrt(inner)
}
