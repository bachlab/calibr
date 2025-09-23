#' Compute Weights for Meta-Analysis
#'
#' This function computes weights for combining test statistics in a meta-analytic
#' framework. The weights are based on sample sizes, degrees of freedom, and
#' optionally provided sum of squares. If sum of squares are not provided,
#' the function estimates weights assuming equal variance.
#'
#' @param t_statistics A numeric vector of t-statistics for each dataset.
#' @param sample_sizes A numeric vector of sample sizes corresponding to each dataset.
#' @param sum_squares A numeric vector of sum of squares for each dataset.
#'                    Defaults to `NA`, in which case weights are approximated
#'                    assuming equal variance.
#'
#' @return A numeric vector of computed weights, normalized such that their sum of squares equals 1.
#' @export
#'
#' @examples
#' # Example 1: Compute weights with sum_squares provided
#' t_statistics <- c(2.5, 1.8, 2.0)
#' sample_sizes <- c(30, 40, 35)
#' sum_squares <- c(5, 6, 5.5)
#' compute_weights(t_statistics, sample_sizes, sum_squares)
#'
#' # Example 2: Compute weights without sum_squares (approximate method)
#' compute_weights(t_statistics, sample_sizes)
compute_weights <- function(t_statistics, sample_sizes, sum_squares = rep(NA, length(sample_sizes))) {
  # Check if t_statistics and sample_sizes are of equal length
  if (length(t_statistics) != length(sample_sizes)) {
    stop("Error: 't_statistics' and 'sample_sizes' must be of equal length.")
  }
  if (length(sum_squares) != length(sample_sizes)) {
    stop("Error: Length of 'sum_squares' must be equal to length of 'sample_sizes'.")
  }

  # Degrees of freedom
  nu_k <- sample_sizes - 2

  # H function
  H <- function(z) as.numeric(sqrt(z) * Rmpfr::igamma(z - 0.5, 0) / Rmpfr::igamma(z, 0))

  # Approximate sum_squares if not provided (partial information case)
  if (any(is.na(sum_squares))) {
    weights <- sqrt(sample_sizes / sum(sample_sizes)) # Approximate ss_k
  } else {
    # Calculate dk
    dk <- t_statistics / (H(nu_k / 2) * sqrt(sum_squares))

    # Calculate variance component vk
    vk <- (H(nu_k / 2)^(-2) * nu_k / (sum_squares * (nu_k - 2))) +
      ((nu_k / (nu_k - 2)) * H(nu_k / 2)^(-2) - 1) * dk^2

    # Compute weights
    weights <- sqrt(1 / vk) / sqrt(sum(1 / vk))
  }

  return(weights)
}
