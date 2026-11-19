#' Compute Strict Pairwise Meta-Analytic Log Bayes Factors
#'
#' Computes **strict, pairwise meta-analytic log Bayes Factors** (natural log scale)
#' for multiple measurement types across datasets.
#' Only datasets containing both measurements are included for each pair.
#' For each pair, the meta-analytic log Bayes factor of each measurement against the null
#' (intercept-only) model is computed on these shared datasets, exactly as in
#' \code{meta_bf_from_data()}, and the two are subtracted:
#' \eqn{\log BF_{12} = \log BF_{10} - \log BF_{20}}.
#'
#' @param data A data frame containing measurements, standard scores, and optionally subject IDs.
#' @param dataset_col Column name identifying datasets.
#' @param measurement_name_col Column name identifying measurement types.
#' @param standard_scores_col Column name for the standard scores (predictor / grouping variable).
#' @param measurement_values_col Column name for measurement values.
#' @param subject_id_col (Optional) Column name identifying unique subjects.
#' @param verbose Logical; if TRUE, prints information on skipped datasets.
#'
#' @return A list with:
#' \describe{
#'   \item{log_bf_matrix}{Antisymmetric matrix of pairwise meta-analytic log Bayes factors (natural log),
#'         \eqn{logBF_{ij} = -logBF_{ji}}. Positive values favour the row measurement over the column one.}
#'   \item{n_datasets_used}{Matrix: number of datasets contributing to each pairwise meta-analysis.}
#'   \item{n_samples_used}{Matrix: per pair, the sum over shared datasets of the smaller of the two
#'         measurements' effective sample sizes.}
#' }
#'
#' @details
#' - Works entirely in **log-space** to prevent overflow and maintain numerical stability.
#' - Datasets where either measurement yields a non-finite t-statistic or log Bayes factor are skipped
#'   for that pair (reported when \code{verbose = TRUE}).
#' - The result does not depend on the order of rows or measurements in \code{data}.
#' - With complete overlap (every dataset contains both measurements), the result equals the
#'   corresponding entry of \code{calibration_summary()$log_bayes_factor_diffs}.
#'
#' @seealso \code{\link{meta_bf_from_data}}, \code{\link{tstat2logbf}}, \code{\link{compute_weights}}
#'
#' @examples
#' # Example with synthetic data (see full package docs)
#' @export
strict_meta_bf_from_data <- function(
    data,
    dataset_col,
    measurement_name_col,
    standard_scores_col,
    measurement_values_col,
    subject_id_col = NULL,
    verbose = FALSE
) {
  required_cols <- c(dataset_col, measurement_name_col, standard_scores_col, measurement_values_col)
  if (!all(required_cols %in% names(data))) {
    stop("Error: one or more required columns do not exist in the dataframe.")
  }

  measurements <- unique(data[[measurement_name_col]])
  datasets     <- unique(data[[dataset_col]])
  n_measurements <- length(measurements)

  log_bf_matrix    <- matrix(NA, n_measurements, n_measurements, dimnames = list(measurements, measurements))
  n_datasets_matrix<- matrix(0,  n_measurements, n_measurements, dimnames = list(measurements, measurements))
  n_samples_matrix <- matrix(0,  n_measurements, n_measurements, dimnames = list(measurements, measurements))

  # Weighted combination of per-dataset t-statistics into a meta-analytic log BF against the
  # null (identical to the combination step in meta_bf_from_data)
  combine_logbf <- function(t_stats, sample_sizes, sum_squares) {
    weights <- compute_weights(t_stats, sample_sizes, sum_squares)
    combined_n <- sum(sample_sizes)
    tstat2logbf(sum(weights * t_stats), combined_n, g = combined_n)
  }

  for (i in seq_len(n_measurements - 1)) {
    for (j in (i + 1):n_measurements) {

      m1 <- measurements[i]
      m2 <- measurements[j]

      # Per-dataset statistics of each measurement, restricted to datasets containing both
      t1s <- n1s <- ss1s <- numeric(0)
      t2s <- n2s <- ss2s <- numeric(0)
      datasets_with_both <- 0
      samples_with_both  <- 0

      for (ds in datasets) {
        subset <- data[data[[dataset_col]] == ds, ]
        if (!all(c(m1, m2) %in% subset[[measurement_name_col]])) next

        # Measurement 1
        s1 <- subset[subset[[measurement_name_col]] == m1, ]
        stats1 <- extract_stats(
          s1[[standard_scores_col]],
          s1[[measurement_values_col]],
          if (!is.null(subject_id_col)) s1[[subject_id_col]]
        )

        # Measurement 2
        s2 <- subset[subset[[measurement_name_col]] == m2, ]
        stats2 <- extract_stats(
          s2[[standard_scores_col]],
          s2[[measurement_values_col]],
          if (!is.null(subject_id_col)) s2[[subject_id_col]]
        )

        t1 <- stats1$t_statistic; n1 <- stats1$sample_size; ss1 <- stats1$sum_squares
        t2 <- stats2$t_statistic; n2 <- stats2$sample_size; ss2 <- stats2$sum_squares
        if (any(!is.finite(c(t1, t2, n1, n2)))) next

        if (!is.finite(tstat2logbf(t1, n1, g = n1)) || !is.finite(tstat2logbf(t2, n2, g = n2))) {
          if (verbose)
            message(sprintf("Skipping dataset '%s' for pair (%s, %s): invalid log BFs.", ds, m1, m2))
          next
        }

        t1s <- c(t1s, t1); n1s <- c(n1s, n1); ss1s <- c(ss1s, ss1)
        t2s <- c(t2s, t2); n2s <- c(n2s, n2); ss2s <- c(ss2s, ss2)
        datasets_with_both <- datasets_with_both + 1
        samples_with_both  <- samples_with_both + min(n1, n2)
      }

      if (datasets_with_both == 0) next

      # Meta-analytic log BF of each measurement against the null (as in meta_bf_from_data),
      # computed on the shared datasets only. Since both are relative to the same null model,
      # their difference is the log BF of m1 versus m2: log BF_12 = log BF_10 - log BF_20.
      meta_logbf <- combine_logbf(t1s, n1s, ss1s) - combine_logbf(t2s, n2s, ss2s)

      log_bf_matrix[i, j]     <- meta_logbf
      log_bf_matrix[j, i]     <- -meta_logbf
      n_datasets_matrix[i, j] <- datasets_with_both
      n_datasets_matrix[j, i] <- datasets_with_both
      n_samples_matrix[i, j]  <- samples_with_both
      n_samples_matrix[j, i]  <- samples_with_both
    }
  }

  if (all(n_datasets_matrix[upper.tri(n_datasets_matrix)] == 0))
    stop("No pair of measurements has overlapping data.")

  list(
    log_bf_matrix    = log_bf_matrix,
    n_datasets_used  = n_datasets_matrix,
    n_samples_used   = n_samples_matrix
  )
}
