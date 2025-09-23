#' Meta-analytic Bayes Factors by Measurement (with optional subject overlap & per-dataset weights)
#'
#' Computes **meta-analytic Bayes Factors (BFs)** for each measurement across datasets by:
#' (i) extracting a per-dataset t-statistic and ancillary stats via `extract_stats()`,
#' (ii) combining those t-statistics with weights from `compute_weights()`,
#' and (iii) turning the combined t into a Bayes factor via `tstat2bf()`.
#' It also returns a retrodiction score (a Pearson correlation implied by the combined t via `tstat2rscore()`).
#' If `subject_id_col` is provided, the function additionally reports the **count of shared subjects** between
#' measurement pairs (based on the dataset–subject combinations), which can help assess the comparability of BFs.
#'
#' This function relies on helper routines that must exist in your environment:
#' `extract_stats()`, `compute_weights()`, `tstat2bf()`, and `tstat2rscore()`.
#'
#' @param data A data frame containing at least the measurement names, standard scores, and measurement values.
#' @param dataset_col A string naming the column identifying datasets/studies. If `NULL` or missing in `data`,
#'   all rows are treated as belonging to a single dataset (a new column `"dataset_col"` is created with
#'   value `"all_data"`), and a warning is emitted.
#' @param measurement_name_col A string naming the column with the measurement identifier (e.g., task/metric name).
#' @param standard_scores_col A string naming the column with the standard scores (often a binary 0/1 score).
#' @param measurement_values_col A string naming the column with the measurement values to be related to the standard scores.
#' @param subject_id_col Optional string naming the column with unique subject identifiers. If supplied, the function
#'   will compute a subject-overlap matrix across measurements (counts of shared subjects across datasets).
#'
#' @details
#' For each measurement:
#' * The function loops over datasets and calls `extract_stats(standard_scores, measurement_values, subject_ids)`
#'   to obtain a per-dataset `t_statistic`, `sample_size`, and `sum_squares`.
#' * It then obtains per-dataset weights via `compute_weights(t_statistics, sample_sizes, sum_squares)` and forms
#'   a **combined t-statistic** as the weighted sum `sum(weights * t_statistics)`.
#' * The Bayes factor is computed by `tstat2bf(combined_t_statistic, sum(sample_sizes), g = sum(sample_sizes))`
#'   (i.e., using total sample size for both the sample-size argument and the `g` hyperparameter).
#' * The **retrodiction score** (a Pearson correlation) is computed via
#'   `tstat2rscore(combined_t_statistic, sum(sample_sizes))`.
#'
#' When `subject_id_col` is provided, subjects are indexed by the pair \{dataset, subject\}; the returned
#' `subject_overlap` matrix reports the **number of shared subjects** between each pair of measurements
#' (NA if undefined because a measurement has no subjects).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{meta_results}{A data frame with one row per measurement containing:
#'     \itemize{
#'       \item `measurement`: Measurement identifier.
#'       \item `combined_bayes_factor`: Meta-analytic Bayes factor computed from the combined t.
#'       \item `log_combined_bayes_factor`: Base-2 logarithm of the Bayes factor.
#'       \item `retrodiction_score`: Pearson correlation implied by the combined t and total sample size.
#'       \item `sample_size`: Total sample size across datasets for the measurement.
#'       \item `combined_t_statistic`: Weighted sum of per-dataset t-statistics.
#'     }
#'   }
#'   \item{subject_overlap}{If `subject_id_col` is supplied, a square matrix whose \eqn{(i,j)} entry is the
#'     **count of shared subjects** between measurement \eqn{i} and \eqn{j} (NA if undefined).}
#'   \item{weights}{A named list (one element per measurement). Each element is a data frame with columns:
#'     `dataset`, `weight`, `t_statistic`, `sample_size`, and `sum_squares`, giving per-dataset inputs and weights.}
#' }
#'
#' @seealso \code{\link{extract_stats}}, \code{\link{compute_weights}}, \code{\link{tstat2bf}}, \code{\link{tstat2rscore}}
#' @export
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#'
#' # Toy helper to fabricate two studies with two measurements (A,B) and a binary standard score
#' make_study <- function(id, n = 80) {
#'   subj <- data.frame(
#'     study_id       = id,
#'     subject_id     = 1:n,
#'     standard_score = rbinom(n, 1, 0.5)
#'   )
#'   # Each subject contributes both measurements; effects differ by measurement
#'   A <- transform(subj,
#'                  measurement_type = "A",
#'                  measurement_value = rnorm(n, mean = 0.8 * standard_score, sd = 1))
#'   B <- transform(subj,
#'                  measurement_type = "B",
#'                  measurement_value = rnorm(n, mean = 0.4 * standard_score, sd = 1))
#'   rbind(A, B)
#' }
#'
#' all_studies <- rbind(make_study("S1", 80), make_study("S2", 100))
#'
#' # Requires user-defined helpers: extract_stats(), compute_weights(), tstat2bf(), tstat2rscore()
#' out <- meta_bf_from_data(
#'   data = all_studies,
#'   dataset_col = "study_id",
#'   measurement_name_col = "measurement_type",
#'   standard_scores_col = "standard_score",
#'   measurement_values_col = "measurement_value",
#'   subject_id_col = "subject_id"
#' )
#'
#' out$meta_results      # main summary table
#' out$subject_overlap   # shared-subject counts between measurements
#' out$weights$A         # per-dataset weights/inputs for measurement A
#' }
meta_bf_from_data <- function(data, dataset_col = NULL, measurement_name_col,
                              standard_scores_col, measurement_values_col,
                              subject_id_col = NULL) {

  required_cols <- c(measurement_name_col, standard_scores_col, measurement_values_col)
  if (!all(required_cols %in% names(data))) {
    stop("Error: one or more columns do not exist in the dataframe.")
  }

  if (is.null(dataset_col) || !(dataset_col %in% names(data))) {
    warning("Warning: all data will be treated as coming from one dataset.")
    dataset_col <- "dataset_col"
    data[[dataset_col]] <- "all_data"
  }

  if (!is.null(subject_id_col) && !(subject_id_col %in% names(data))) {
    stop("Error: subject id column does not exist in the dataframe.")
  }

  unique_measurements <- unique(data[[measurement_name_col]])

  results <- data.frame(
    measurement = character(),
    combined_bayes_factor = numeric(),
    combined_t_statistic = numeric(),
    stringsAsFactors = FALSE
  )

  weights_out <- stats::setNames(vector("list", length(unique_measurements)), unique_measurements)

  for (measurement in unique_measurements) {

    measurement_data <- data[data[[measurement_name_col]] == measurement, ]

    unique_datasets <- unique(measurement_data[[dataset_col]])
    n_datasets <- length(unique_datasets)

    t_statistics <- numeric(n_datasets)
    sample_sizes <- numeric(n_datasets)
    sum_squares <- numeric(n_datasets)

    for (i in seq_along(unique_datasets)) {
      subset_data <- measurement_data[measurement_data[[dataset_col]] == unique_datasets[i], ]

      standard_scores <- subset_data[[standard_scores_col]]
      measurement_values <- subset_data[[measurement_values_col]]
      subject_ids <- if (!is.null(subject_id_col)) subset_data[[subject_id_col]] else NULL

      extracted_stats <- tryCatch({
        extract_stats(standard_scores, measurement_values, subject_ids)
      }, error = function(e) {
        stop(sprintf("extract_stats failed for measurement '%s' in '%s': %s",
                     measurement, unique_datasets[i], e$message))
      })

      t_statistics[i] <- extracted_stats$t_statistic
      sample_sizes[i] <- extracted_stats$sample_size
      sum_squares[i] <- extracted_stats$sum_squares
    }

    weights <- compute_weights(t_statistics, sample_sizes, sum_squares)

    combined_t_statistic <- sum(weights * t_statistics)

    combined_bayes_factor <- tstat2bf(combined_t_statistic, sum(sample_sizes), g = sum(sample_sizes))

    results <- rbind(results, data.frame(
      measurement = measurement,
      combined_bayes_factor = combined_bayes_factor,
      log_combined_bayes_factor = log2(combined_bayes_factor),
      retrodiction_score = tstat2rscore(combined_t_statistic, sum(sample_sizes)),
      sample_size = sum(sample_sizes),
      combined_t_statistic = combined_t_statistic,
      stringsAsFactors = FALSE
    ))

    weights_out[[measurement]] <- data.frame(
      dataset = unique_datasets,
      weight = as.numeric(weights),
      t_statistic = as.numeric(t_statistics),
      sample_size = as.numeric(sample_sizes),
      sum_squares = as.numeric(sum_squares),
      stringsAsFactors = FALSE
    )
  }

  subject_overlap <- NULL

  if (!is.null(subject_id_col) &&
      subject_id_col %in% names(data) &&
      !is.null(dataset_col) &&
      dataset_col %in% names(data)) {

    data$combined_subject_id <- paste(data[[dataset_col]], data[[subject_id_col]], sep = "_")

    subjects_by_measurement <- split(data$combined_subject_id, data[[measurement_name_col]])
    subjects_by_measurement <- lapply(subjects_by_measurement, unique)

    measurement_names <- names(subjects_by_measurement)
    n <- length(measurement_names)
    subject_overlap <- matrix(NA, nrow = n, ncol = n,
                              dimnames = list(measurement_names, measurement_names))

    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        subjects_i <- subjects_by_measurement[[i]]
        subjects_j <- subjects_by_measurement[[j]]
        shared <- length(intersect(subjects_i, subjects_j))
        total <- length(union(subjects_i, subjects_j))
        subject_overlap[i, j] <- ifelse(total == 0, NA, shared)
      }
    }

    data$combined_subject_id <- NULL
  }

  list(
    meta_results = results[order(results$log_combined_bayes_factor, decreasing = TRUE), ],
    subject_overlap = subject_overlap,
    weights = weights_out
  )
}
