test_that("weights sum to one when squared", {
  t_stats <- c(1.5, 3)
  sample_sizes <- c(50, 50)
  weights <- compute_weights(t_stats, sample_sizes)

  expect_equal(sum(weights^2), 1, tolerance = 1e-6)
})

test_that("function throws an error when t_statistics and sample_sizes have different lengths", {
  expect_error(compute_weights(c(2, 2, 2), c(10, 10)),
               "'t_statistics' and 'sample_sizes' must be of equal length" , fixed = TRUE)
})

test_that("equal sample sizes return equal weights", {
  weights <- compute_weights(c(2, 2), c(10, 10))

  expect_equal(weights[1], weights[2], tolerance = 1e-6)
})

test_that("larger sample size produces a larger weight", {
  weights <- compute_weights(c(2, 2), c(10, 5))

  expect_gt(weights[1], weights[2])
})
