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

  # Error: Bayes factor supplied without effect direction (BFs carry no sign)
  expect_error(
    meta_bf_from_stats(
      n_datasets = 2,
      t_statistics = c(NA, 2.1),
      sample_sizes = c(20, 25),
      sum_squares = c(4.5, 5),
      bayes_factors = c(3, NA)
    ),
    "effect_directions[1] must be 1 or -1",
    fixed = TRUE
  )
})

test_that("meta_bf_from_stats returns a natural-log BF consistent with tstat2logbf", {
  t_statistics <- c(2.5, 1.8, 2.0)
  sample_sizes <- c(30, 40, 35)
  sum_squares  <- c(5, 6, 5.5)
  out <- meta_bf_from_stats(3, t_statistics, sample_sizes, sum_squares)
  w <- compute_weights(t_statistics, sample_sizes, sum_squares)
  N <- sum(sample_sizes)
  expect_equal(out$meta_results$log_combined_bayes_factor,
               tstat2logbf(sum(w * t_statistics), N, g = N), tolerance = 1e-12)
  expect_false("combined_bayes_factor" %in% names(out$meta_results))
})

test_that("meta_bf_from_stats applies effect directions to Bayes-factor inputs", {
  t_statistics <- c(2.5, 1.8, 2.0)
  sample_sizes <- c(30, 40, 35)
  sum_squares  <- c(5, 6, 5.5)
  bayes_factors <- exp(tstat2logbf(t_statistics, sample_sizes))
  from_t  <- meta_bf_from_stats(3, t_statistics, sample_sizes, sum_squares)
  from_bf <- meta_bf_from_stats(3, rep(NA, 3), sample_sizes, sum_squares,
                                bayes_factors = bayes_factors, effect_directions = c(1, 1, 1))
  flipped <- meta_bf_from_stats(3, rep(NA, 3), sample_sizes, sum_squares,
                                bayes_factors = bayes_factors, effect_directions = c(1, -1, 1))
  expect_equal(from_bf$meta_results$log_combined_bayes_factor,
               from_t$meta_results$log_combined_bayes_factor, tolerance = 1e-8)
  expect_lt(flipped$meta_results$combined_t_statistic, from_bf$meta_results$combined_t_statistic)
})
