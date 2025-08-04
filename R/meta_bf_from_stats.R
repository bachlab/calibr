#' Compute Meta-Analytic Bayes Factor From Stats
#'
#' This function computes a **meta-analytic Bayes Factor (BF)** by combining test statistics
#' from multiple datasets using a weighted approach. The function allows input of either
#' t-statistics or precomputed Bayes Factors, transforming them when necessary.
#' It ensures valid statistical combination by weighting individual test statistics appropriately.
#'
#' @param n_datasets An integer specifying the number of datasets.
#' @param t_statistics A numeric vector of t-statistics for each dataset.
#'                     Defaults to `NA` for all entries, meaning Bayes Factors
#'                     must be provided if t-statistics are missing.
#' @param sample_sizes A numeric vector of sample sizes for each dataset.
#' @param sum_squares A numeric vector of sum of squares for each dataset.
#' @param bayes_factors A numeric vector of Bayes Factors for each dataset.
#'                      Defaults to `NA`. If `t_statistics` are provided,
#'                      these are ignored; otherwise, the function will
#'                      transform Bayes Factors into t-statistics.
#'
#' @return A list with two elements:
#'   \item{combined_bayes_factor}{The computed **meta-analytic Bayes Factor**.}
#'   \item{combined_t_statistic}{The combined test statistic used to derive the Bayes Factor.}
#' @export
#'
#' @examples
#' # Example: Compute meta-analytic Bayes Factor from t-statistics
#' n_datasets <- 3
#' t_statistics <- c(2.5, 1.8, 2.0)
#' sample_sizes <- c(30, 40, 35)
#' sum_squares <- c(5, 6, 5.5)
#' meta_bf_from_stats(n_datasets, t_statistics, sample_sizes, sum_squares)
#'
#' # Example: Compute meta-analytic Bayes Factor using Bayes Factors as input
#' bayes_factors <- c(10, 5, 8)  # BF values instead of t-statistics
#' meta_bf_from_stats(n_datasets, bayes_factors = bayes_factors, sample_sizes, sum_squares)
meta_bf_from_stats <- function(n_datasets, t_statistics = rep(NA, n_datasets), sample_sizes, sum_squares = rep(NA, n_datasets), bayes_factors = rep(NA, n_datasets)) {

  # Ensure all input vectors have the same length
  if (length(t_statistics) != length(bayes_factors) ||
      length(t_statistics) != length(sample_sizes) ||
      length(t_statistics) != length(sum_squares)) {
    stop("Error: All input vectors (t_statistics, sample_sizes, sum_squares, bayes_factors) must have the same length.")
  }

  # Transform Bayes Factors into test statistics where needed
  for (i in seq_along(t_statistics)) {
    if (is.na(t_statistics[i]) && is.na(bayes_factors[i])) {
      stop(paste("Error: Both t_statistic and bayes_factor are missing for index", i))
    }
    if (is.na(t_statistics[i]) && !is.na(bayes_factors[i])) {
      t_statistics[i] <- bf2tstat(bayes_factors[i], sample_sizes[i], g = sample_sizes[i])
    }
  }

  # Compute weights
  weights <- compute_weights(t_statistics, sample_sizes, sum_squares)

  # Compute the combined meta-analytic test statistic
  combined_t_statistic <- sum(weights * t_statistics)

  # Convert combined statistic to a Bayes Factor
  combined_bayes_factor <- tstat2bf(combined_t_statistic, sum(sample_sizes), g = sum(sample_sizes))

  return(list(combined_bayes_factor = combined_bayes_factor, combined_t_statistic = combined_t_statistic))
}
