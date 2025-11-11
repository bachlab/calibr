test_that("t-statistic is preserved through Bayes Factor conversion", {

  # Define input parameters
  t_stat <- 2.5  # Original t-statistic
  sample_size <- 30  # Sample size

  # Step 1: Convert t-statistic to Bayes Factor
  log_bayes_factor <- tstat2logbf(t_stat, sample_size)

  # Step 2: Convert Bayes Factor back to t-statistic
  recovered_t_stat <- logbf2tstat(log_bayes_factor, sample_size)

  # Step 3: Check that the recovered t-statistic is approximately equal to the original
  expect_equal(recovered_t_stat, t_stat, tolerance = 1e-6)
})
