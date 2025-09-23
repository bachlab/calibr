#' Compute Meta-Analytic Bayes Factor From Stats
#'
#' This function computes a **meta-analytic Bayes Factor (BF)** by combining test statistics
#' from multiple datasets using a weighted approach. You may supply per-dataset
#' \(t\)-statistics directly, or provide Bayes Factors to be transformed into
#' \(t\)-statistics where needed. The combined statistic is then converted into a Bayes Factor.
#'
#' To keep outputs uniform with \code{meta_bf_from_data}, the return value includes:
#' \itemize{
#'   \item \code{meta_results}: a 1-row data frame with columns
#'         \code{measurement}, \code{combined_bayes_factor}, \code{log_combined_bayes_factor},
#'         \code{retrodiction_score} (Pearson \(r\) implied by the combined \(t\)),
#'         \code{sample_size} (sum of per-dataset \(n\)), and \code{combined_t_statistic}.
#'   \item \code{subject_overlap}: always \code{NULL} (not applicable when only stats are provided).
#'   \item \code{weights}: a named list (single element, named by \code{measurement_name}) containing
#'         a per-dataset table with columns \code{dataset}, \code{weight}, \code{t_statistic},
#'         \code{sample_size}, \code{sum_squares}.
#' }
#'
#' This function relies on helper routines available in your environment:
#' \code{compute_weights()}, \code{tstat2bf()}, \code{bf2tstat()}, and \code{tstat2rscore()}.
#'
#' @param n_datasets Integer; number of datasets.
#' @param t_statistics Numeric vector of per-dataset \(t\)-statistics. Defaults to \code{rep(NA, n_datasets)}.
#'   Where an entry is \code{NA}, the corresponding \code{bayes_factors} value must be provided and will be
#'   transformed to a \(t\)-statistic.
#' @param sample_sizes Numeric vector of per-dataset sample sizes (\(n\)).
#' @param sum_squares Numeric vector of per-dataset sums of squares (or analogous quantity used by your weighting scheme).
#'   Defaults to \code{rep(NA, n_datasets)} if not used by \code{compute_weights()}.
#' @param bayes_factors Numeric vector of per-dataset Bayes Factors. Defaults to \code{rep(NA, n_datasets)}.
#'   Ignored where \code{t_statistics} is non-missing; otherwise each BF is converted to a \(t\)-statistic via
#'   \code{bf2tstat(bf, n, g = n)}.
#' @param measurement_name Optional single string used to label the output row and the weights list name.
#'   Defaults to \code{"overall"}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{meta_results}{Data frame with one row and columns:
#'     \code{measurement}, \code{combined_bayes_factor}, \code{log_combined_bayes_factor},
#'     \code{retrodiction_score}, \code{sample_size}, \code{combined_t_statistic}.}
#'   \item{subject_overlap}{Always \code{NULL}.}
#'   \item{weights}{Named list of length 1 (name = \code{measurement_name}) containing a per-dataset
#'     data frame with columns \code{dataset}, \code{weight}, \code{t_statistic}, \code{sample_size}, \code{sum_squares}.}
#' }
#' @export
#'
#' @examples
#' # Example 1: Compute from t-statistics
#' n_datasets   <- 3
#' t_statistics <- c(2.5, 1.8, 2.0)
#' sample_sizes <- c(30, 40, 35)
#' sum_squares  <- c(5, 6, 5.5)
#' out1 <- meta_bf_from_stats(
#'     n_datasets,
#'     t_statistics,
#'     sample_sizes,
#'     sum_squares,
#'     measurement_name = "A"
#' )
#' out1$meta_results
#' out1$weights$A
#'
#' # Example 2: Compute using Bayes Factors (converted to t where t is NA)
#' bayes_factors <- c(10, 5, 8)
#' out2 <- meta_bf_from_stats(n_datasets, t_statistics = c(NA, NA, NA),
#'                            sample_sizes = sample_sizes, sum_squares = sum_squares,
#'                            bayes_factors = bayes_factors, measurement_name = "overall")
#' out2$meta_results
#'
meta_bf_from_stats <- function(
    n_datasets,
    t_statistics   = rep(NA, n_datasets),
    sample_sizes,
    sum_squares    = rep(NA, n_datasets),
    bayes_factors  = rep(NA, n_datasets),
    measurement_name = "overall"
) {

  # --- Input checks -----------------------------------------------------------
  # lengths consistent with n_datasets
  lens <- c(length(t_statistics), length(sample_sizes), length(sum_squares), length(bayes_factors))
  if (any(lens != n_datasets)) {
    stop("Error: All input vectors (t_statistics, sample_sizes, sum_squares, bayes_factors) must have length n_datasets.")
  }
  # pairwise length equality as in original
  if (length(t_statistics) != length(bayes_factors) ||
      length(t_statistics) != length(sample_sizes) ||
      length(t_statistics) != length(sum_squares)) {
    stop("Error: All input vectors (t_statistics, sample_sizes, sum_squares, bayes_factors) must have the same length.")
  }
  if (!is.character(measurement_name) || length(measurement_name) != 1) {
    stop("Error: measurement_name must be a single string.")
  }

  # --- Fill missing t's using BF -> t where needed ---------------------------
  for (i in seq_along(t_statistics)) {
    if (is.na(t_statistics[i]) && is.na(bayes_factors[i])) {
      stop(paste("Error: Both t_statistic and bayes_factor are missing for index", i))
    }
    if (is.na(t_statistics[i]) && !is.na(bayes_factors[i])) {
      t_statistics[i] <- bf2tstat(bayes_factors[i], sample_sizes[i], g = sample_sizes[i])
    }
  }

  # --- Weights and combined statistics --------------------------------------
  weights <- compute_weights(t_statistics, sample_sizes, sum_squares)
  combined_t_statistic <- sum(weights * t_statistics)

  total_n <- sum(sample_sizes)
  combined_bayes_factor <- tstat2bf(combined_t_statistic, total_n, g = total_n)
  retrodiction_score <- tstat2rscore(combined_t_statistic, total_n)

  # --- Assemble outputs to mirror meta_bf_from_data --------------------------
  weights_df <- data.frame(
    dataset      = seq_len(n_datasets),
    weight       = as.numeric(weights),
    t_statistic  = as.numeric(t_statistics),
    sample_size  = as.numeric(sample_sizes),
    sum_squares  = as.numeric(sum_squares),
    stringsAsFactors = FALSE
  )

  meta_results <- data.frame(
    measurement               = measurement_name,
    combined_bayes_factor     = combined_bayes_factor,
    log_combined_bayes_factor = log2(combined_bayes_factor),
    retrodiction_score        = retrodiction_score,
    sample_size               = total_n,
    combined_t_statistic      = combined_t_statistic,
    stringsAsFactors = FALSE
  )

  list(
    meta_results   = meta_results,
    subject_overlap = NULL,
    weights        = stats::setNames(list(weights_df), measurement_name)
  )
}
