#' Compute Strict Pairwise Meta-Analytic Bayes Factors for Multiple Measurements
#'
#' This function computes **strict, pairwise meta-analytic Bayes Factors (BFs)** for measurement types across multiple datasets,
#' using only those datasets in which \emph{both} measures of a pair are present and valid. For each pair of measures,
#' it meta-analytically combines the within-dataset Bayes factor contrasts, properly handling cases where the BF
#' cannot be mapped to a real t-statistic.
#'
#' Returns symmetric matrices of pairwise Bayes factors, log Bayes factors, number of datasets used per pair,
#' and total sample size used per pair.
#'
#' @param data A dataframe containing the measurements, standard scores, and optionally subject IDs across datasets.
#' @param dataset_col A string specifying the column name that identifies datasets.
#' @param measurement_name_col A string specifying the column name that identifies different measurements.
#' @param standard_scores_col A string specifying the column name containing the standard scores (predictor or grouping variable).
#' @param measurement_values_col A string specifying the column name containing the measurement values.
#' @param subject_id_col (Optional) A string specifying the column name that identifies unique subjects.
#'        If provided, subject-level statistics will be extracted if your extract_stats supports it.
#' @param verbose Logical; if TRUE, messages about dataset inclusion/exclusion are printed.
#'
#' @return A list containing:
#'   \describe{
#'     \item{bf_matrix}{A symmetric matrix of pairwise meta-analytic Bayes factors (BFs) for all measurement types.}
#'     \item{log_bf_matrix}{A matrix of pairwise log2 Bayes factors. log(BF_ij) = -log(BF_ji).}
#'     \item{n_datasets_used}{A matrix: number of datasets contributing to each pairwise meta-analysis.}
#'     \item{n_samples_used}{A matrix: total sample size contributing to each pairwise meta-analysis (sum of n for each dataset in which both are present and valid).}
#'   }
#'
#' @details
#' Only datasets in which both measurement types are present are used for each pair. For each such dataset,
#' Bayes factors for each measure (against the null) are computed and their ratio is mapped back to a t-statistic.
#' Datasets in which this mapping fails (i.e., where the Bayes factor ratio is too low to produce a real t-statistic)
#' are counted as "no evidence" (t = 0). The function then performs a meta-analytic weighted combination of the
#' resulting t-statistics, and reports the combined Bayes factor.
#'
#' @seealso \code{\link{meta_bf_from_data}}, \code{\link{tstat2bf}}, \code{\link{bf2tstat}}
#'
#' @examples
#' # Load required library
#' require(dplyr)
#'
#' # Set seed for reproducibility
#' set.seed(123)
#'
#' # Function to generate a study dataset
#' generate_study <- function(study_id, n_per_group, within = TRUE) {
#'   df <- expand.grid(
#'     study_id = study_id,           # Study ID
#'     subject_id = 1:n_per_group,    # Unique subject ID for within-subject design
#'     standard_score = c(0, 1),      # Binary standard scores
#'     measurement_type = c("A", "B") # Two measurement types
#'   )
#'   # In a between-subject design, each subject only has one score
#'   if (!within) {
#'     df$subject_id <- 1:nrow(df)  # Ensure unique subject IDs
#'   }
#'
#'   # Generate measurement values based on standard score and measurement type
#'   df$measurement <- with(df, rnorm(nrow(df),
#'                                    mean = ifelse(measurement_type == "A",
#'                                                  ifelse(standard_score == 1, 10, 5),
#'                                                  ifelse(standard_score == 1, 7, 6)),
#'                                    sd = 5))
#'   return(df)
#' }
#'
#' # Generate example studies
#' all_studies <- do.call(rbind, list(
#'   generate_study("Study1", 30, within = TRUE),   # Within-subject study
#'   generate_study("Study2", 30, within = TRUE),   # Another within-subject study
#'   generate_study("Study3", 30, within = FALSE)   # Between-subject study
#' ))
#'
#' # Compute the strict pairwise meta-analytic Bayes Factor summary
#' result <- strict_meta_bf_from_data(
#'   data = all_studies,
#'   dataset_col = "study_id",
#'   standard_scores_col = "standard_score",
#'   measurement_name_col = "measurement_type",
#'   measurement_values_col = "measurement",
#'   subject_id_col = "subject_id"
#' )
#'
#' print(result$bf_matrix)        # Pairwise Bayes Factor matrix
#' print(result$log_bf_matrix)    # Log2 Bayes Factor matrix
#' print(result$n_datasets_used)  # Datasets used per pair
#' print(result$n_samples_used)   # Sample size used per pair
#'
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
  datasets <- unique(data[[dataset_col]])
  n_measurements <- length(measurements)

  bf_matrix        <- matrix(NA, nrow = n_measurements, ncol = n_measurements, dimnames = list(measurements, measurements))
  log_bf_matrix    <- matrix(NA, nrow = n_measurements, ncol = n_measurements, dimnames = list(measurements, measurements))
  n_datasets_matrix<- matrix(0,  nrow = n_measurements, ncol = n_measurements, dimnames = list(measurements, measurements))
  n_samples_matrix <- matrix(0,  nrow = n_measurements, ncol = n_measurements, dimnames = list(measurements, measurements))

  for (i in seq_len(n_measurements - 1)) {
    for (j in (i + 1):n_measurements) {
      m1 <- measurements[i]
      m2 <- measurements[j]
      t_stats <- c()
      sample_sizes <- c()
      sum_squares <- c()
      datasets_with_both <- 0
      samples_with_both <- 0

      for (ds in datasets) {
        subset <- data[data[[dataset_col]] == ds, ]
        if (!all(c(m1, m2) %in% subset[[measurement_name_col]])) next

        # m1
        s1 <- subset[subset[[measurement_name_col]] == m1, ]
        standard_scores1 <- s1[[standard_scores_col]]
        values1 <- s1[[measurement_values_col]]
        ids1 <- if (!is.null(subject_id_col)) s1[[subject_id_col]] else NULL
        stats1 <- extract_stats(standard_scores1, values1, ids1)
        t1 <- stats1$t_statistic
        n1 <- stats1$sample_size
        ss1 <- stats1$sum_squares

        # m2
        s2 <- subset[subset[[measurement_name_col]] == m2, ]
        standard_scores2 <- s2[[standard_scores_col]]
        values2 <- s2[[measurement_values_col]]
        ids2 <- if (!is.null(subject_id_col)) s2[[subject_id_col]] else NULL
        stats2 <- extract_stats(standard_scores2, values2, ids2)
        t2 <- stats2$t_statistic
        n2 <- stats2$sample_size
        ss2 <- stats2$sum_squares

        n_used <- min(n1, n2)

        bf1 <- tstat2bf(t1, n1, g = n1)
        bf2 <- tstat2bf(t2, n2, g = n2)
        if (!is.finite(bf1) || !is.finite(bf2) || bf1 <= 0 || bf2 <= 0) {
          if (verbose) message(sprintf("Skipping dataset '%s' for pair (%s, %s) due to invalid BFs: bf1=%.4f, bf2=%.4f", ds, m1, m2, bf1, bf2))
          next
        }

        bf_ratio <- bf1 / bf2
        t_contrast <- suppressWarnings(bf2tstat(bf_ratio, n_used, g = n_used))
        if (!is.finite(t_contrast) || is.nan(t_contrast)) {
          if (verbose) message(sprintf("Skipping dataset '%s' for pair (%s, %s): bf2tstat returned NaN/Inf", ds, m1, m2))
          next
        }

        t_stats <- c(t_stats, t_contrast)
        sample_sizes <- c(sample_sizes, n_used)
        sum_squares <- c(sum_squares, ss1 + ss2)
        datasets_with_both <- datasets_with_both + 1
        samples_with_both <- samples_with_both + n_used
      }

      if (datasets_with_both == 0) next

      weights <- compute_weights(t_stats, sample_sizes, sum_squares)
      combined_t_stat <- sum(weights * t_stats)
      combined_sample_size <- sum(sample_sizes)
      meta_bf <- tstat2bf(combined_t_stat, combined_sample_size, g = combined_sample_size)
      meta_logbf <- log2(meta_bf)

      # Fill both upper and lower triangle to enforce symmetry
      bf_matrix[i, j]         <- meta_bf
      bf_matrix[j, i]         <- ifelse(is.finite(meta_bf) && meta_bf > 0, 1 / meta_bf, NA)
      log_bf_matrix[i, j]     <- meta_logbf
      log_bf_matrix[j, i]     <- ifelse(is.finite(meta_logbf), -meta_logbf, NA)
      n_datasets_matrix[i, j] <- datasets_with_both
      n_datasets_matrix[j, i] <- datasets_with_both
      n_samples_matrix[i, j]  <- samples_with_both
      n_samples_matrix[j, i]  <- samples_with_both
    }
  }
  # Fail only if NO pair of measurements has overlap
  if (all(n_datasets_matrix[upper.tri(n_datasets_matrix)] == 0)) {
    stop("No pair of measurements has overlapping data. Cannot compute any meta-analytic Bayes factor.")
  }
  return(list(
    bf_matrix        = bf_matrix,
    log_bf_matrix    = log_bf_matrix,
    n_datasets_used  = n_datasets_matrix,
    n_samples_used   = n_samples_matrix
  ))
}
