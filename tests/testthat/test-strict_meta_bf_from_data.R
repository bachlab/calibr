test_that("strict_meta_bf_from_data stops when there is no overlap for any pair", {
  # Create data where A and B are never in the same dataset
  df <- data.frame(
    dataset_col = rep(1:4, each = 10),
    measurement_name_col = rep(c("A", "B"), each = 20),
    standard_scores_col = rnorm(40),
    measurement_values_col = rnorm(40)
  )
  expect_error(
    strict_meta_bf_from_data(
      data = df,
      dataset_col = "dataset_col",
      measurement_name_col = "measurement_name_col",
      standard_scores_col = "standard_scores_col",
      measurement_values_col = "measurement_values_col"
    ),
    "No pair of measurements has overlapping data."
  )
})

test_that("strict_meta_bf_from_data returns NA for pairs with no overlap, but works for others", {
  # Dataset 1: A and B
  df1 <- data.frame(
    dataset_col = 1,
    measurement_name_col = rep(c("A", "B"), each = 5),
    standard_scores_col = rnorm(10),
    measurement_values_col = rnorm(10)
  )
  # Dataset 2: B and C
  df2 <- data.frame(
    dataset_col = 2,
    measurement_name_col = rep(c("B", "C"), each = 5),
    standard_scores_col = rnorm(10),
    measurement_values_col = rnorm(10)
  )
  # Dataset 3: only C
  df3 <- data.frame(
    dataset_col = 3,
    measurement_name_col = rep("C", 5),
    standard_scores_col = rnorm(5),
    measurement_values_col = rnorm(5)
  )

  df <- rbind(df1, df2, df3)

  # Only "A" and "B" will have overlap if you manipulate accordingly
  result <- strict_meta_bf_from_data(
    data = df,
    dataset_col = "dataset_col",
    measurement_name_col = "measurement_name_col",
    standard_scores_col = "standard_scores_col",
    measurement_values_col = "measurement_values_col"
  )
  expect_true(is.list(result))
  expect_true(any(is.na(result$log_bf_matrix)))
})
