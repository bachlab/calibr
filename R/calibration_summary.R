#' Compute a Summary of Meta-Analytic Bayes Factors and Effect Sizes
#'
#' This function computes a calibration summary by:
#' - Running `meta_bf_from_data()` for **Bayes Factors**.
#' - Running `meta_cohensd()` for **Cohen's d and Hedges' g** (only if `standard_scores_col` is binary).
#' - Returning a **matrix of log Bayes Factor differences** between measurements.
#' - Identifying **the best measurement** based on Bayes Factor.
#' - Optionally returning the **subject overlap matrix** (percentage of shared data between measurements).
#'
#' @param data A dataframe containing multiple datasets, measurement types, and measurement values.
#' @param dataset_col A string specifying the column name that identifies datasets - if NULL, all data pertains to one dataset.
#' @param measurement_name_col A string specifying the column name that identifies different measurements.
#' @param measurement_values_col A string specifying the column name containing measurement values.
#' @param standard_scores_col A string specifying the column name containing standard scores.
#' @param subject_id_col A string specifying the column name containing subject IDs - can be NULL.
#' @param attenuation_factor A numeric value (0 < x ≤ 1) indicating the expected proportional transfer
#'        of the calibration effect to the target context. For example, `attenuation_factor = 0.9`
#'        assumes real-world effects are expected to be 90% as strong as those observed in calibration.
#'        Only applicable to binary designs.
#' @param power Power of test (1 minus Type II error probability); used to compute sample size. See pwr::pwr.t.test.
#' @param sig.level Significance level (Type I error probability); used to compute sample size. See pwr::pwr.t.test.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default),
#'        "greater" or "less"; used to compute sample size. See pwr::pwr.t.test.
#'
#' @return A list containing:
#'   \describe{
#'     \item{summary_table}{A dataframe summarizing results from Bayes Factor and effect size computations.}
#'     \item{log_bayes_factor_diffs}{A matrix of Bayes Factor log-ratio differences between all measurements.}
#'     \item{best_measurement}{A descriptive summary of the best-performing measurement.}
#'     \item{subject_overlap}{A matrix showing the percentage of shared subjects between measurement types,
#'         based on combined dataset-subject identifiers. Returned only if subject_id_col is provided and multiple datasets exist.}
#'   }
#'
#' @details
#' - If `standard_scores_col` is **binary (0/1)**, effect sizes (Cohen's d, Hedges' g) will be computed. If not only Bayes Factors will be computed.
#' - The `attenuation_factor` rescales the inferred meta-analytic Cohen’s d to reflect the expected
#'   proportion of the calibration effect that generalizes to a new or target context.
#'   For example, setting `attenuation_factor = 0.8` means that only 80% of the calibration effect
#'   is assumed to hold in the substantive experiment, leading to larger required sample sizes.
#'
#' @examples
#' require(dplyr)
#' set.seed(123)
#'
#' generate_study <- function(study_id, n_per_group, within = TRUE) {
#'   df <- expand.grid(
#'     study_id = study_id,
#'     subject_id = 1:n_per_group,
#'     standard_score = c(0, 1),
#'     measurement_type = c("A", "B")
#'   )
#'   if (!within) df$subject_id <- 1:nrow(df)
#'   df$measurement <- with(df, rnorm(nrow(df),
#'     mean = ifelse(measurement_type == "A",
#'                   ifelse(standard_score == 1, 10, 5),
#'                   ifelse(standard_score == 1, 7, 6)),
#'     sd = 5))
#'   return(df)
#' }
#'
#' all_studies <- do.call(rbind, list(
#'   generate_study("Study1", 30, within = TRUE),
#'   generate_study("Study2", 30, within = TRUE),
#'   generate_study("Study3", 30, within = FALSE)
#' ))
#'
#' result <- calibration_summary(
#'   data = all_studies,
#'   dataset_col = "study_id",
#'   standard_scores_col = "standard_score",
#'   measurement_name_col = "measurement_type",
#'   measurement_values_col = "measurement",
#'   subject_id_col = "subject_id",
#'   attenuation_factor = 0.9
#' )
#'
#' print(result$summary_table)
#' print(result$log_bayes_factor_diffs)
#' print(result$subject_overlap)
#' print(result$best_measurement)
#'
#' @export
calibration_summary <- function(
    data, dataset_col = NULL, measurement_name_col, measurement_values_col,
    standard_scores_col, subject_id_col = NULL,
    attenuation_factor = 0.9, power = 0.8, sig.level = 0.05,
    alternative = "two.sided") {

  if (is.null(standard_scores_col) || !(standard_scores_col %in% names(data)))
    stop("Error: standard scores do not exist in the dataframe.")

  if (is.null(dataset_col) || !(dataset_col %in% names(data))) {
    warning("Warning: all data will be treated as coming from one dataset.")
    dataset_col <- "dataset_col"
    data[[dataset_col]] <- "all_data"
  }

  unique_standard_scores <- unique(data[[standard_scores_col]])
  is_binary <- length(unique_standard_scores) == 2 && all(unique_standard_scores %in% c(0, 1))

  # Meta-analytic Bayes Factors
  bf_output <- meta_bf_from_data(
    data = data,
    dataset_col = dataset_col,
    measurement_name_col = measurement_name_col,
    standard_scores_col = standard_scores_col,
    measurement_values_col = measurement_values_col,
    subject_id_col = subject_id_col
  )
  bayes_results <- bf_output$meta_results
  subject_overlap <- bf_output$subject_overlap

  # Meta-analytic effect sizes (if binary)
  if (is_binary) {
    effect_size_results <- meta_cohensd(
      data, dataset_col, standard_scores_col, measurement_name_col,
      measurement_values_col, subject_id_col)
    effect_size_table <- effect_size_results$meta_results
  } else {
    effect_size_table <- data.frame(
      measurement = unique(data[[measurement_name_col]]),
      meta_cohen_d = NA, meta_hedge_g = NA,
      stringsAsFactors = FALSE
    )
  }

  # Merge results
  summary_table <- merge(bayes_results, effect_size_table, by = "measurement")
  summary_table$sample_size_required <- NA

  if (is_binary) {
    for (i in seq_len(nrow(summary_table))) {
      meta_d <- summary_table$meta_cohen_d[i]
      if (is.finite(meta_d)) {
        # Expected standardized effect in target context
        d_target <- abs(meta_d) * attenuation_factor

        power_test <- pwr::pwr.t.test(
          d = d_target,
          power = power,
          sig.level = sig.level,
          type = "two.sample",
          alternative = alternative
        )
        summary_table$sample_size_required[i] <- ceiling(power_test$n)
      }
    }
  }

  # Log Bayes factor difference matrix
  measurements <- summary_table$measurement
  log_bayes_factors <- summary_table$log_combined_bayes_factor
  log_bayes_factor_diffs <- outer(log_bayes_factors, log_bayes_factors, "-")
  dimnames(log_bayes_factor_diffs) <- list(measurements, measurements)

  # Identify best measurement
  log_max_bf <- max(log_bayes_factors)
  best_measurement <- measurements[which.max(log_bayes_factors)]
  best_retrodiction <- measurements[which.max(summary_table$retrodiction_score)]
  second_best <- sort(log_bayes_factors, decreasing = TRUE)[2]
  margin <- log_max_bf - second_best

  best_statement <- paste(
    "The highest log Bayes Factor is achieved by", best_measurement,
    "with a value of", round(log_max_bf, 2), ".",
    "When comparing methods, consider the degree of subject overlap — greater overlap strengthens the validity of the comparison.",
    "The highest retrodiction score is observed for", best_retrodiction, ".",
    "If the top Bayes Factor and retrodiction score point to different methods, examine differences in sample sizes."
  )

  list(
    summary_table = summary_table[order(summary_table$log_combined_bayes_factor, decreasing = TRUE),],
    log_bayes_factor_diffs = log_bayes_factor_diffs,
    best_measurement = best_statement,
    subject_overlap = subject_overlap
  )
}
