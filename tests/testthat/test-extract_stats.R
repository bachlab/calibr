test_that("extract_stats throws an error when measurements contain NA values", {

  # Case where measurements contain NA values
  standard_scores <- c(1.2, 2.3, 1.8, 2.7, 3.1)
  measurements_with_na <- c(5.1, 6.2, NA, 7.0, 7.5)

  expect_error(
    extract_stats(standard_scores, measurements_with_na),
    "Error: 'standard_scores' or 'measurement_values' contain NA values.",
    fixed = TRUE
  )

  # Case where standard_scores contain NA values
  standard_scores_with_na <- c(1.2, 2.3, NA, 2.7, 3.1)
  measurements <- c(5.1, 6.2, 2, 7.0, 7.5)

  expect_error(
    extract_stats(standard_scores_with_na, measurements),
    "Error: 'standard_scores' or 'measurement_values' contain NA values.",
    fixed = TRUE
  )

})
