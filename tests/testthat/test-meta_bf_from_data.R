# Mock helper functions required for testing
extract_stats <- function(standard_scores, measurement_values) {
  return(list(
    t_statistic = mean(standard_scores) / (sd(standard_scores) / sqrt(length(standard_scores))),
    sample_size = length(standard_scores),
    sum_squares = sum(standard_scores^2)
  ))
}

compute_weights <- function(t_statistics, sample_sizes, sum_squares) {
  return(sample_sizes / sum(sample_sizes))  # Simple weighting scheme
}

bayes_factor_from_t_statistic <- function(t_statistic, total_sample_size, g) {
  return(exp(t_statistic^2 / (2 * g)))  # Approximate BF calculation
}

# Define a test dataset
test_data <- data.frame(
  dataset_id = rep(1:3, each = 10),
  measurement_name = rep(c("A", "B"), each = 15),
  standard_scores = rnorm(30),
  measurement_values = rnorm(30, mean = 5, sd = 2)
)

# Start testing
test_that("meta_bf_from_data produces correct output structure", {
  result <- meta_bf_from_data(
    data = test_data,
    dataset_col = "dataset_id",
    measurement_name_col = "measurement_name",
    standard_scores_col = "standard_scores",
    measurement_values_col = "measurement_values"
  )

  expect_true(is.list(result))  # Output should be a dataframe
  expect_true(is.data.frame(result$meta_results))
  expect_equal(nrow(result$meta_results), length(unique(test_data$measurement_name)))  # Number of rows should match unique measurements
})

test_that("meta_bf_from_data throws an error when missing columns", {
  expect_warning(
    meta_bf_from_data(
      data = test_data,
      dataset_col = "missing_col",  # Invalid column
      measurement_name_col = "measurement_name",
      standard_scores_col = "standard_scores",
      measurement_values_col = "measurement_values"
    ),
    "Warning: all data will be treated as coming from one dataset."
  )
})

test_that("meta_bf_from_data handles different numbers of datasets per measurement", {
  # Create an imbalanced dataset where one measurement has fewer datasets
  imbalanced_data <- test_data[test_data$dataset_id != 3 | test_data$measurement_name != "A", ]

  result <- meta_bf_from_data(
    data = imbalanced_data,
    dataset_col = "dataset_id",
    measurement_name_col = "measurement_name",
    standard_scores_col = "standard_scores",
    measurement_values_col = "measurement_values"
  )

  expect_equal(nrow(result$meta_results), length(unique(imbalanced_data$measurement_name)))  # Should still match unique measurements
  expect_true(all(result$combined_bayes_factor > 0))  # Ensure Bayes factors are valid numbers
})

test_that("meta_bf_from_data returns non-negative Bayes factors", {
  result <- meta_bf_from_data(
    data = test_data,
    dataset_col = "dataset_id",
    measurement_name_col = "measurement_name",
    standard_scores_col = "standard_scores",
    measurement_values_col = "measurement_values"
  )

  expect_true(all(result$combined_bayes_factor >= 0))  # BF should not be negative
})
