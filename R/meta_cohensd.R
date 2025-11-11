#' Compute Meta-Analytic Cohen's d and Hedges' g for Multiple Measurements
#'
#' This function calculates the meta-analytic **Cohen's d** and **Hedges' g**
#' effect sizes for multiple measurement types across different datasets.
#' It automatically determines whether each dataset follows a **paired (within-subject)**
#' or **unpaired (between-subject) design** based on the presence of repeated subject IDs.
#' Effect sizes are aggregated using inverse variance weighting. For calculations ,we referred to
#' Goulet-Pelletier, J. C., & Cousineau, D. (2018). A review of effect sizes and their confidence intervals, Part I: The Cohen’sd family.
#' The Quantitative Methods for Psychology, 14(4), 242-265.
#'
#' @param data A dataframe containing multiple datasets and measurement types.
#' @param dataset_col A string specifying the column name that identifies datasets.
#' @param standard_scores_col A string specifying the column name that indicates the standard_score variable (binary: 0 or 1). The direction of the difference
#' taken in the computation of Cohen's d is *higher* minus *lower* standard score.
#' @param measurement_name_col A string specifying the column name that identifies different measurement types.
#' @param measurement_values_col A string specifying the column name that contains measurement values (continuous variable).
#' @param subject_id_col (Optional) A string specifying the column name that identifies unique subjects.
#'        If provided, the function determines whether each dataset is paired (within-subject) or unpaired.
#'        Paired analysis only applies if ALL subjects have measurements for both standard_scores.
#'
#' @return A list containing:
#'   \describe{
#'     \item{effect_sizes}{A dataframe containing Cohen’s d and Hedges’ g for each dataset and measurement.}
#'     \item{meta_results}{A dataframe summarizing the meta-analytic effect sizes for each measurement.}
#'   }
#'
#' @details
#' - If `subject_id_col` is **not provided**, all datasets are treated as **unpaired** (between-subject).
#' - If `subject_id_col` **is provided**, the function checks whether subjects appear in both standard_scores:
#'   - If a dataset contains all repeated subjects across standard_scores, it is treated as **paired (within-subject)**.
#'   - Otherwise, it is treated as **unpaired (between-subject)**.
#'
#' - **For both paired and unpaired designs**, Cohen’s d is computed as:
#'   \deqn{ d = \frac{\mathrm{mean_1} - \mathrm{mean_2}}{\mathrm{pooled\_sd}} }
#'   where pooled_sd is the pooled standard deviation. This is to ensure comparability between effect sizes
#'   coming from within and between subjects designs.
#'
#'
#' @examples
#'
#' require(dplyr)
#'
#' set.seed(123)
#'
#' # Function to generate a study with a better (A) measurement in terms of effect size.
#'generate_study <- function(study_id, n_per_group, within = TRUE) {
#'  df <- expand.grid(
#'    study_id = study_id,
#'    subject_id = 1:n_per_group,
#'    standard_score = c(0, 1),
#'    measurement_type = c("A", "B")
#'  )
#'  if(!within){
#'    df$subject_id <- 1:nrow(df)
#'  }
#'  df$measurement <- with(df, rnorm(nrow(df),
#'                                   mean = ifelse(measurement_type == "A",
#'                                   ifelse(standard_score == 1, 10, 5),
#'                                   ifelse(standard_score == 1, 7, 6)),
#'                                   sd = 5))
#'  return(df)
#'}
#'
#'all_studies <- do.call(rbind, list(
#'  generate_study("Study1", 30, within = TRUE),
#'  generate_study("Study2", 30, within = TRUE),
#'  generate_study("Study3", 30, within = FALSE)
#'))
#'
#'result <- meta_cohensd(
#'  data = all_studies, dataset_col = "study_id", standard_scores_col = "standard_score",
#'  measurement_name_col = "measurement_type", measurement_values_col = "measurement",
#'  subject_id_col = "subject_id"
#')
#'
#'print(result)
#'
#' @export
meta_cohensd <- function(data, dataset_col = NULL, standard_scores_col, measurement_name_col,
                         measurement_values_col, subject_id_col = NULL) {


  # Ensure input columns exist in the dataframe
  required_cols <- c(measurement_name_col, standard_scores_col, measurement_values_col)
  if (!all(required_cols %in% names(data)) || (!is.null(dataset_col) && !(dataset_col %in% names(data)))) {
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

  # Extract unique standard_score values
  standard_score_levels <- sort(unique(data[[standard_scores_col]]), decreasing = TRUE)

  # Ensure binary standard_scores
  if (length(standard_score_levels) != 2) {
    stop(sprintf("Error: the standard score variable '%s' is not binary.", standard_scores_col))
  }

  # Identify unique measurement types
  unique_measurements <- unique(data[[measurement_name_col]])

  # Initialize results storage
  effect_sizes <- data.frame()
  meta_results <- data.frame()

  # Loop through each measurement type
  for (measurement in unique_measurements) {

    # Subset data for the current measurement
    measurement_data <- data[data[[measurement_name_col]] == measurement,]

    # Identify unique datasets
    unique_datasets <- unique(measurement_data[[dataset_col]])

    # Loop through each dataset
    for (dataset in unique_datasets) {

      subset_data <- measurement_data[measurement_data[[dataset_col]] == dataset,]

      # Sample sizes
      n1 <- sum(subset_data[[standard_scores_col]] == standard_score_levels[1])
      n2 <- sum(subset_data[[standard_scores_col]] == standard_score_levels[2])

      # Degrees of freedom
      nu <- n1 + n2 - 2

      # Hedge's correction factor (approximation)
      Jh <- 1 - (3 / (4 * nu - 1))

      # Standard deviations
      sd1 <- stats::sd(subset_data[[measurement_values_col]][subset_data[[standard_scores_col]] == standard_score_levels[1]])
      sd2 <- stats::sd(subset_data[[measurement_values_col]][subset_data[[standard_scores_col]] == standard_score_levels[2]])

      # Pooled standard deviation
      pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / nu)


      # Determine if the dataset should be analyzed as paired or unpaired
      is_paired <- !is.null(subject_id_col) && all(table(subset_data[[subject_id_col]]) == 2)

      subset_score_1 <- subset_data[subset_data[[standard_scores_col]] == standard_score_levels[1],]
      subset_score_2 <- subset_data[subset_data[[standard_scores_col]] == standard_score_levels[2],]

      # Check for constant (zero variance) vectors
      sd1_val <- stats::sd(subset_score_1[[measurement_values_col]])
      sd2_val <- stats::sd(subset_score_2[[measurement_values_col]])

      # If either group is constant, handle as unpaired and issue warning
      if (sd1_val == 0 || sd2_val == 0) {
        warning(sprintf(
          "Treated as unpaired for dataset '%s', measurement '%s' due to sd = 0 for either group.",
          as.character(dataset), as.character(measurement)
        ))
        is_paired <- FALSE
      }

      mean1 <- mean(subset_data[[measurement_values_col]][subset_data[[standard_scores_col]] == standard_score_levels[1]])
      mean2 <- mean(subset_data[[measurement_values_col]][subset_data[[standard_scores_col]] == standard_score_levels[2]])

      cohen_d <- (mean1 - mean2) / pooled_sd

      if (is_paired) {
        # Paired (within-subject) case
        rho <- stats::cor(subset_score_1[[measurement_values_col]], subset_score_2[[measurement_values_col]])
        var_d <- (nu / (nu - 2)) * ((2 * (1 - rho)) / n1) * (1 + ((cohen_d^2) * n1)/(2 * (1 - rho))) - cohen_d^2/Jh^2
      } else {
        # Unpaired (between-subject) case
        var_d <- (nu / (nu - 2)) * ((n1 + n2)/(n1*n2)) * (1 + (cohen_d^2)/((n1 + n2)/(n1*n2))) - cohen_d^2/Jh^2
      }

      hedge_g = Jh * cohen_d
      sample_size <- n1 + n2

      # Store effect sizes
      effect_sizes <- rbind(effect_sizes, data.frame(
        measurement = measurement, dataset = dataset, cohen_d = cohen_d, hedge_g = hedge_g, var_d = as.numeric(var_d), is_paired = is_paired,
        sample_size = sample_size,
        stringsAsFactors = FALSE
      ))

    }

    # Here, subset the cohen d's relative to the current measurement
    subset_effect_sizes <- effect_sizes[effect_sizes[['measurement']] == measurement,]

    # Compute meta-analytic estimates
    meta_cohen_d <- sum(subset_effect_sizes$cohen_d / subset_effect_sizes$var_d) / sum(1 / subset_effect_sizes$var_d)
    meta_hedge_g <- sum(subset_effect_sizes$hedge_g / subset_effect_sizes$var_d) / sum(1 / subset_effect_sizes$var_d)

    # Add measurement to meta results record
    meta_results <- rbind(meta_results, data.frame(
      measurement = measurement,
      meta_cohen_d = meta_cohen_d,
      meta_hedge_g = meta_hedge_g,
      stringsAsFactors = FALSE
    ))
  }

  rownames(effect_sizes) <- NULL

  return(list(effect_sizes = effect_sizes, meta_results = meta_results[order(meta_results$meta_cohen_d, decreasing = TRUE), ]))
}
