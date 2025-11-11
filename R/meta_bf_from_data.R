#' Meta-analytic Bayes Factors by Measurement (with optional subject overlap & per-dataset weights)
#'
#' Computes **meta-analytic Bayes Factors (BFs)** and a **meta-analytic retrodiction score**
#' (a correlation measure) for each measurement across datasets.
#'
#' For each measurement:
#' * Extracts per-dataset t-statistics and ancillary stats via `extract_stats()`.
#' * Combines t-statistics using `compute_weights()` to form a weighted t (`sum(weights * t_i)`).
#' * Computes the Bayes factor from that combined t using `tstat2bf()`.
#' * Computes the **meta-analytic retrodiction score** by converting each t to its
#'   corresponding correlation (\eqn{r_i = t_i / \sqrt{t_i^2 + \mathrm{df}_i}}),
#'   Fisher-z transforming each, weighting by their residual degrees of freedom
#'   (\eqn{w_i = \mathrm{df}_i - 1}), averaging, and back-transforming with \eqn{tanh()}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{meta_results}{Data frame with one row per measurement containing:
#'     `measurement`, `log_combined_bayes_factor`,
#'     `retrodiction_score`, `sample_size`, and `combined_t_statistic`.}
#'   \item{subject_overlap}{Matrix of shared subjects between measurements, if applicable.}
#'   \item{weights}{List of per-dataset weight tables for each measurement.}
#' }
#' @seealso \code{\link{extract_stats}}, \code{\link{compute_weights}}, \code{\link{tstat2bf}}
#' @export
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
  results <- data.frame()
  weights_out <- stats::setNames(vector("list", length(unique_measurements)), unique_measurements)

  for (measurement in unique_measurements) {
    measurement_data <- data[data[[measurement_name_col]] == measurement, ]
    unique_datasets <- unique(measurement_data[[dataset_col]])
    n_datasets <- length(unique_datasets)

    t_statistics <- numeric(n_datasets)
    sample_sizes <- numeric(n_datasets)
    sum_squares  <- numeric(n_datasets)

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
      sum_squares[i]  <- extracted_stats$sum_squares
    }

    weights <- compute_weights(t_statistics, sample_sizes, sum_squares)
    combined_t_statistic <- sum(weights * t_statistics)
    log_combined_bayes_factor <- tstat2logbf(combined_t_statistic, sum(sample_sizes), g = sum(sample_sizes))

    # --- Proper meta-analytic retrodiction score using Fisher-z ---
    df_i <- sample_sizes - 2
    r_i  <- t_statistics / sqrt(t_statistics^2 + df_i)
    z_i  <- atanh(pmax(pmin(r_i, 0.999999), -0.999999))
    w_i  <- pmax(df_i - 1, 1)  # effective Fisher-z weights from residual df
    z_bar <- sum(w_i * z_i) / sum(w_i)
    r_meta <- tanh(z_bar)

    results <- rbind(results, data.frame(
      measurement = measurement,
      log_combined_bayes_factor = log_combined_bayes_factor,
      retrodiction_score = r_meta,
      sample_size = sum(sample_sizes),
      combined_t_statistic = combined_t_statistic,
      stringsAsFactors = FALSE
    ))

    weights_out[[measurement]] <- data.frame(
      dataset = unique_datasets,
      weight = as.numeric(weights),
      t_statistic = as.numeric(t_statistics),
      sample_size = as.numeric(sample_sizes),
      df = as.numeric(df_i),
      r_i = as.numeric(r_i),
      sum_squares = as.numeric(sum_squares),
      stringsAsFactors = FALSE
    )
  }

  subject_overlap <- NULL
  if (!is.null(subject_id_col) &&
      subject_id_col %in% names(data) &&
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
        total  <- length(union(subjects_i, subjects_j))
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
