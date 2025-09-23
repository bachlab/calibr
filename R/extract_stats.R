#' Extract Statistical Information from a Regression Model
#'
#' This function extracts key statistics from a regression model where `standard_scores`
#' are predicted from `measurement_values`. It supports both simple linear regression
#' (when no `subject_ids` is provided) and mixed-effects modeling with subject-specific
#' intercepts (when `subject_ids` is provided). The t-statistics are computed using the
#' number of subejcts minus 2 degrees of freedom.
#'
#' @param standard_scores A numeric vector containing the dependent variable (outcome).
#' @param measurement_values A numeric vector containing the independent variable (predictor).
#' @param subject_ids (Optional) A vector indicating subject IDs. If provided, the function
#'        fits a mixed-effects model with a **random intercept** for each subject.
#'
#' @return A list containing:
#'   \describe{
#'     \item{t_statistic}{The t-statistic for `measurement_values` (slope of the predictor).}
#'     \item{sum_squares}{The sum of squares for `measurement_values`, used for weighting.}
#'     \item{sample_size}{The total number of observations included in the model.}
#'   }
#'
#' @details
#' - If `subject_ids` is **not provided**, the function runs a **simple linear regression**:#'
#' - If `subject_ids` is **provided**, the function subtracts the mean of measurement_values for a given subject_ids.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   subject_ids = rep(1:10, each = 5),
#'   standard_scores = rnorm(50),
#'   measurement_values = rnorm(50, mean = 5, sd = 2)
#' )
#'
#' # Example 1: Simple regression (no subject effects)
#' extract_stats(df$standard_scores, df$measurement_values)
#'
#' # Example 2: Mixed-effects model with subject-specific intercepts
#' extract_stats(df$standard_scores, df$measurement_values, df$subject_ids)
#'
#' @export
extract_stats <- function(standard_scores, measurement_values, subject_ids = NULL) {

  # Check for NA values in measurement_values
  if (any(is.na(measurement_values)) || any(is.na(standard_scores)) ) {
    stop("Error: 'standard_scores' or 'measurement_values' contain NA values.")
  }

  # Convert inputs into a data frame
  data <- data.frame(standard_scores, measurement_values)

  # Boolean for data pairedness
  repeated_measures <- !is.null(subject_ids) && length(unique(table(subject_ids))) == 1 && all(table(subject_ids) > 1)

  if (!is.null(subject_ids)) {
    data$subject_ids <- factor(subject_ids)
  } else {
    data$subject_ids <- factor(seq(1, nrow(data)))
  }

  unique_subjects <- unique(data$subject_ids)

  # Subtract mean if there are multiple subject measurement
  if (repeated_measures) {
    data$measurement_values_centered <- data$measurement_values - stats::ave(data$measurement_values, data$subject_ids, FUN = mean)
  } else {
    data$measurement_values_centered <- data$measurement_values  # No centering if not repeated measures
  }

  # - - - - Module for computing the regression t-statistic

  x <- data$measurement_values_centered
  y <- data$standard_scores

  N <- length(x)  # or N <- nrow(data)

  # Pearson r
  r <- stats::cor(x, y)  # add use="complete.obs" if you ever allow NAs

  # this is identical (up to floating-point rounding) to summary(lm(y ~ x)) slope t
  t_statistic <- r * sqrt((N - 2) / (1 - r^2))

  # Compute sum of squares and sample size
  sum_squares <- sum((measurement_values - mean(measurement_values))^2)
  sample_size <- length(unique_subjects)

  return(list(t_statistic = t_statistic, sum_squares = sum_squares, sample_size = sample_size))

}



