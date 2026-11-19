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
  expect_true(all(is.finite(result$meta_results$log_combined_bayes_factor)))  # Ensure log Bayes factors are valid numbers
})

test_that("meta_bf_from_data returns finite log Bayes factors", {
  result <- meta_bf_from_data(
    data = test_data,
    dataset_col = "dataset_id",
    measurement_name_col = "measurement_name",
    standard_scores_col = "standard_scores",
    measurement_values_col = "measurement_values"
  )

  expect_true(all(is.finite(result$meta_results$log_combined_bayes_factor)))  # log BF must be a finite number
})

test_that("meta_bf_from_data log Bayes factor is in natural log (Eq. 5 of the combined t)", {
  result <- meta_bf_from_data(
    data = test_data,
    dataset_col = "dataset_id",
    measurement_name_col = "measurement_name",
    standard_scores_col = "standard_scores",
    measurement_values_col = "measurement_values"
  )
  w <- result$weights$A
  combined_t <- sum(w$weight * w$t_statistic)
  N <- sum(w$sample_size)
  raw_bf <- (1 + N)^((N - 2) / 2) * (1 + N / (combined_t^2 / (N - 2) + 1))^(-(N - 1) / 2)
  expect_equal(result$meta_results$log_combined_bayes_factor[result$meta_results$measurement == "A"],
               log(raw_bf), tolerance = 1e-8)
})
