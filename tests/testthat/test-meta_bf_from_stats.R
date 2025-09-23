test_that("meta_bf_from_stats throws expected errors", {

  # Error: Input vectors must have the same length
  expect_error(
    meta_bf_from_stats(
      n_datasets = 3,
      t_statistics = c(2.5, 1.8),  # Length 2 instead of 3
      sample_sizes = c(30, 40, 35),
      sum_squares = c(5, 6, 5.5),
      bayes_factors = c(NA, NA, NA)
    ),
    "Error: All input vectors (t_statistics, sample_sizes, sum_squares, bayes_factors) must have length n_datasets.",
    fixed = TRUE
  )

  # Error: Both t_statistic and bayes_factor are missing for an index
  expect_error(
    meta_bf_from_stats(
      n_datasets = 2,
      t_statistics = c(NA, 2.1),
      sample_sizes = c(20, 25),
      sum_squares = c(4.5, 5),
      bayes_factors = c(NA, NA)  # First index has both missing
    ),
    "Error: Both t_statistic and bayes_factor are missing for index 1",
    fixed = TRUE
  )
})
