# Define a test dataset with multiple measurements
set.seed(123)
test_data <- data.frame(
  dataset = rep(c("Study1", "Study2", "Study3"), each = 60),  # Three datasets
  measurement_name = rep(c("Measure_A", "Measure_B"), each = 30, times = 6),  # Two measurement types
  group = rep(c(0, 1), each = 15, times = 6),  # Binary group: 0 = control, 1 = treatment
  measurement_value = rep(c(rnorm(15, mean = 35, sd = 10),rnorm(15, mean = 65, sd = 10)), each = 1,times = 6) + rnorm(180, mean = 0, sd = 5)
)

test_that("meta_cohensd produces correct output structure", {
  result <- meta_cohensd(
    test_data,
    "dataset", "group", "measurement_name", "measurement_value"
  )

  expect_true(is.list(result))  # Output should be a list
  expect_true(all(c("effect_sizes", "meta_results") %in% names(result)))

  # Check effect_sizes dataframe
  expect_true(is.data.frame(result$effect_sizes))
  expect_true(all(c("measurement", "dataset", "cohen_d", "hedge_g") %in% names(result$effect_sizes)))

  # Check meta_results dataframe
  expect_true(is.data.frame(result$meta_results))
  expect_true(all(c("measurement", "meta_cohen_d", "meta_hedge_g") %in% names(result$meta_results)))
})

test_that("meta_cohensd computes reasonable effect sizes", {
  result <- meta_cohensd(
    test_data,
    "dataset", "group", "measurement_name", "measurement_value"
  )

  # Check that each measurement's meta-analytic Cohen's d is positive
  expect_true(all(result$meta_results$meta_cohen_d > 0))
  expect_true(all(result$meta_results$meta_hedge_g > 0))
})

test_that("meta_cohensd throws an error when group column is not binary", {
  test_data_invalid <- test_data
  test_data_invalid$group <- rep(0:2, length.out = nrow(test_data))  # Introduce a third group

  expect_error(
    meta_cohensd(test_data_invalid, "dataset", "group", "measurement_name", "measurement_value"),
    "Error: the standard score variable group is not binary."
  )
})

test_that("meta_cohensd handles missing columns correctly", {
  expect_error(
    meta_cohensd(test_data, "dataset", "wrong_group", "measurement_name", "measurement_value"),
    "Error: one or more columns do not exist in the dataframe."
  )

  expect_error(
    meta_cohensd(test_data, "wrong_dataset", "group", "measurement_name", "measurement_value"),
    "Error: one or more columns do not exist in the dataframe."
  )
})

test_that("meta_cohensd returns non-negative variances", {
  result <- meta_cohensd(
    test_data,
    "dataset", "group", "measurement_name", "measurement_value"
  )

  expect_true(all(result$effect_sizes$var_d >= 0))  # Variances should not be negative
  expect_true(all(result$effect_sizes$var_g >= 0))  # Variances should not be negative
})

test_that("meta_cohensd handles multiple measurements correctly", {
  result <- meta_cohensd(
    test_data,
    "dataset", "group", "measurement_name", "measurement_value"
  )

  # Ensure all unique measurements are in the results
  unique_measurements <- unique(test_data$measurement_name)
  expect_true(all(unique_measurements %in% result$meta_results$measurement))
})

