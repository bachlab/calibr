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
    "Error: the standard score variable 'group' is not binary."
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

# Paired dataset: 30 subjects measured once at each standard score level
set.seed(42)
paired_base <- rnorm(30)
paired_data <- data.frame(
  dataset = "D1", measurement_name = "A", subject = rep(1:30, 2), group = rep(c(0, 1), each = 30),
  measurement_value = c(paired_base + rnorm(30, 0, 0.3), paired_base + 0.5 + rnorm(30, 0, 0.3))
)

test_that("meta_cohensd paired variance does not depend on row order (subjects matched by ID)", {
  ref <- meta_cohensd(paired_data, "dataset", "group", "measurement_name", "measurement_value", "subject")
  shuffled <- paired_data[sample(nrow(paired_data)), ]
  res <- meta_cohensd(shuffled, "dataset", "group", "measurement_name", "measurement_value", "subject")
  expect_true(ref$effect_sizes$is_paired)
  expect_equal(res$effect_sizes$var_d, ref$effect_sizes$var_d, tolerance = 1e-12)
  expect_equal(res$meta_results$meta_cohen_d, ref$meta_results$meta_cohen_d, tolerance = 1e-12)
})

test_that("meta_cohensd treats repeated rows in the same condition as unpaired", {
  same_cond <- data.frame(
    dataset = "D1", measurement_name = "A", subject = rep(1:30, each = 2),
    group = rep(rep(c(0, 1), each = 15), each = 2), measurement_value = rnorm(60)
  )
  res <- suppressWarnings(
    meta_cohensd(same_cond, "dataset", "group", "measurement_name", "measurement_value", "subject")
  )
  expect_false(res$effect_sizes$is_paired)
})

test_that("meta_cohensd warns when the BF treats a dataset as within-subject but d does not", {
  incomplete <- paired_data[!(paired_data$subject %in% 26:30 & paired_data$group == 0), ]
  expect_warning(
    meta_cohensd(incomplete, "dataset", "group", "measurement_name", "measurement_value", "subject"),
    "treats this dataset as within-subject"
  )
  # No such warning (nor any other output) for perfectly paired data
  expect_silent(
    meta_cohensd(paired_data, "dataset", "group", "measurement_name", "measurement_value", "subject")
  )
})

