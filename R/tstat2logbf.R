#' Compute Log Bayes Factor from a t-Statistic
#'
#' Computes the **log Bayes Factor** (natural logarithm) for a given \emph{t}-statistic
#' and sample size under the Zellner–Siow \(g\)-prior framework.
#' This version is numerically stable and returns results on the log scale,
#' which avoids overflow for very large Bayes Factors.
#'
#' The model corresponds to Equation (5) in the calibration paper:
#' \deqn{\log BF_{10} = \frac{1}{2}\left[(n-2)\log(1+g) - (n-1)\log\left(1+\frac{g}{t^2/(n-2)+1}\right)\right]}
#'
#' @param t_statistic Numeric vector of t-statistics.
#' @param sample_size Numeric vector of sample sizes (must have same length as \code{t_statistic},
#'        or length 1 for recycling).
#' @param g Prior scaling parameter (defaults to \code{g = sample_size}).
#'
#' @details
#' This function returns \eqn{\log BF_{10}} (natural log).
#' To convert to base-2 logarithms, divide by \code{log(2)}.
#' To obtain the raw Bayes Factor, exponentiate the result, i.e. \code{exp(log_bf)}.
#'
#' @return Numeric vector of log Bayes Factors (\eqn{\log BF_{10}}).
#'
#' @examples
#' # Example: Compute log BF for t = 2.5, n = 30
#' tstat2bf(2.5, 30)
#'
#' # Example: Compute log BF with custom g
#' tstat2bf(2.5, 30, g = 15)
#'
#' # Example: Vectorized usage
#' tstat2bf(c(2, 3, 5), c(30, 30, 30))
#'
#' @export
tstat2logbf <- function(t_statistic, sample_size, g = sample_size) {
  # --- input checks
  if (!is.numeric(t_statistic) || !is.numeric(sample_size) || !is.numeric(g)) {
    stop("Error: all inputs must be numeric.")
  }

  # vector recycling
  len <- max(length(t_statistic), length(sample_size), length(g))
  if (length(t_statistic) != len) t_statistic <- rep(t_statistic, length.out = len)
  if (length(sample_size) != len) sample_size <- rep(sample_size, length.out = len)
  if (length(g) != len) g <- rep(g, length.out = len)

  # validity check
  valid <- is.finite(t_statistic) & is.finite(sample_size) & (sample_size > 2)
  out <- rep(NA_real_, len)

  if (any(valid)) {
    t2 <- t_statistic[valid]^2
    n <- sample_size[valid]
    gg <- g[valid]

    # numerically stable computation
    log_bf <- 0.5 * ((n - 2) * log1p(gg) -
                       (n - 1) * log1p(gg / (t2 / (n - 2) + 1)))

    out[valid] <- log_bf
  }

  return(out)
}
