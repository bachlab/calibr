#' Compute Meta-Analytic Bayes Factors for Multiple Measurements with Subject-Level Effects
#'
#' This function computes **meta-analytic Bayes Factors (BFs)** for multiple
#' measurements across different datasets while optionally accounting for subject-level effects.
#' It extracts statistics from each dataset for each measurement type and applies a weighted combination of
#' t-statistics. If subject IDs are provided and there are multiple datasets, the function also computes
#' the percentage of shared subjects between each pair of measurement types, serving as a proxy for the comparability
#' of Bayes Factors across measurements.
#'
#' @param data A dataframe containing measurements, standard scores, and optionally subject IDs across datasets.
#' @param dataset_col A string specifying the column name that identifies datasets.
#' @param measurement_name_col A string specifying the column name that identifies different measurements.
#' @param standard_scores_col A string specifying the column name containing the standard scores.
#' @param measurement_values_col A string specifying the column name containing the measurement values.
#' @param subject_id_col (Optional) A string specifying the column name that identifies unique subjects.
#'        If provided, and multiple datasets are present, the function computes subject-level overlap between measurements.
#'
#' @return A list containing:
#'   \describe{
#'     \item{meta_results}{A dataframe with the following columns for each measurement:
#'       \itemize{
#'         \item `measurement`: The measurement identifier.
#'         \item `combined_bayes_factor`: The computed meta-analytic Bayes Factor.
#'         \item `log_combined_bayes_factor`: Base 2 logarithm of the Bayes Factor.
#'         \item `retrodiction_score`: Estimated Pearson correlation between the measurement and standard scores.
#'         \item `sample_size`: Total sample size used for the measurement.
#'         \item `combined_t_statistic`: The combined t-statistic used to compute the Bayes Factor.
#'       }}
#'     \item{subject_overlap}{A square matrix indicating the amount of shared subjects (based on combined dataset and subject ID)
#'       between each pair of measurements. Returned only if `subject_id_col` is provided.}
#'   }
#' @export
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
#' # Compute the meta-analytic Bayes Factor summary
#' result <- meta_bf_from_data(
#'   data = all_studies,
#'   dataset_col = "study_id",
#'   standard_scores_col = "standard_score",
#'   measurement_name_col = "measurement_type",
#'   measurement_values_col = "measurement",
#'   subject_id_col = "subject_id"
#' )
#'
#' print(result$meta_results)        # Main Bayes Factor summary
#' print(result$subject_overlap)     # Optional subject overlap matrix

meta_bf_from_data <- function(data, dataset_col = NULL, measurement_name_col,
                              standard_scores_col, measurement_values_col,
                              subject_id_col = NULL) {

  # Ensure input columns exist in the dataframe
  required_cols <- c(measurement_name_col, standard_scores_col, measurement_values_col)
  if (!all(required_cols %in% names(data))) {
    stop("Error: one or more columns do not exist in the dataframe.")
  }

  # Clarify dataset situation
  if (is.null(dataset_col) || !(dataset_col %in% names(data))) {
    warning("Warning: all data will be treated as coming from one dataset.")
    dataset_col <- "dataset_col"
    data[[dataset_col]] <- "all_data"
  }

  # If subject_id_col is provided, ensure it exists
  if (!is.null(subject_id_col) && !(subject_id_col %in% names(data))) {
    stop("Error: subject id column does not exist in the dataframe.")
  }

  # Identify unique measurements
  unique_measurements <- unique(data[[measurement_name_col]])

  # Initialize results storage
  results <- data.frame(
    measurement = character(),
    combined_bayes_factor = numeric(),
    combined_t_statistic = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop through each measurement type
  for (measurement in unique_measurements) {

    # Subset data for the current measurement
    measurement_data <- data[data[[measurement_name_col]] == measurement, ]

    # Identify unique datasets
    unique_datasets <- unique(measurement_data[[dataset_col]])
    n_datasets <- length(unique_datasets)

    # Initialize vectors for extracted statistics
    t_statistics <- numeric(n_datasets)
    sample_sizes <- numeric(n_datasets)
    sum_squares <- numeric(n_datasets)

    # Compute statistics for each dataset
    for (i in seq_along(unique_datasets)) {
      subset_data <- measurement_data[measurement_data[[dataset_col]] == unique_datasets[i], ]

      # Extract variables
      standard_scores <- subset_data[[standard_scores_col]]
      measurement_values <- subset_data[[measurement_values_col]]

      # If subject_id_col is provided, extract it; otherwise, pass NULL
      subject_ids <- if (!is.null(subject_id_col)) subset_data[[subject_id_col]] else NULL

      # Compute statistics using extract_stats
      extracted_stats <- tryCatch({
        extract_stats(standard_scores, measurement_values, subject_ids)
      }, error = function(e) {
        stop(sprintf("extract_stats failed for measurement '%s' in '%s': %s",
                     measurement, unique_datasets[i] , e$message))
      })

      # Store extracted statistics
      t_statistics[i] <- extracted_stats$t_statistic
      sample_sizes[i] <- extracted_stats$sample_size
      sum_squares[i] <- extracted_stats$sum_squares
    }

    # Compute weights
    weights <- compute_weights(t_statistics, sample_sizes, sum_squares)

    # Compute the combined meta-analytic test statistic
    combined_t_statistic <- sum(weights * t_statistics)

    # Convert combined statistic to a Bayes Factor
    combined_bayes_factor <- tstat2bf(combined_t_statistic, sum(sample_sizes), g = sum(sample_sizes))

    # Store results
    results <- rbind(results, data.frame(
      measurement = measurement,
      combined_bayes_factor = combined_bayes_factor,
      log_combined_bayes_factor = log2(combined_bayes_factor),
      retrodiction_score = tstat2rscore(combined_t_statistic, sum(sample_sizes)),
      sample_size = sum(sample_sizes),
      combined_t_statistic = combined_t_statistic,
      stringsAsFactors = FALSE
    ))
  }

  # ---- Compute subject overlap matrix (only if subject IDs) ----
  subject_overlap <- NULL  # default

  if (!is.null(subject_id_col) &&
      subject_id_col %in% names(data) &&
      !is.null(dataset_col) &&
      dataset_col %in% names(data)) {

    # Create a unique subject identifier by combining dataset and subject ID
    data$combined_subject_id <- paste(data[[dataset_col]], data[[subject_id_col]], sep = "_")

    # Get unique subject IDs per measurement type
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

    # Optional cleanup
    data$combined_subject_id <- NULL
  }

  # Return both meta results and overlap matrix
  return(list(
    meta_results = results[order(results$log_combined_bayes_factor, decreasing = TRUE), ],
    subject_overlap = subject_overlap
  ))
}
