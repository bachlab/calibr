# tests/testthat/test-tstat2logbf.R

test_that("tstat2logbf computes correct log Bayes factors and handles inputs safely", {
  # 1. Basic expected values (compared to manual computation)
  t <- 2.5
  n <- 30
  g <- 30
  manual_logbf <- 0.5 * ((n - 2) * log1p(g) -
                           (n - 1) * log1p(g / (t^2 / (n - 2) + 1)))

  expect_equal(tstat2logbf(t, n, g), manual_logbf, tolerance = 1e-10)

  # 2. Default g = n should match explicit g = n
  expect_equal(tstat2logbf(t, n), tstat2logbf(t, n, g = n), tolerance = 1e-10)

  # 3. Vectorized input should work
  t_vec <- c(2, 3, 5)
  n_vec <- c(30, 30, 30)
  out_vec <- tstat2logbf(t_vec, n_vec)
  expect_length(out_vec, 3)
  expect_true(all(is.finite(out_vec)))

  # 5. Large t-statistics do not overflow (log BF stays finite)
  large_t <- 1e3
  expect_true(is.finite(tstat2logbf(large_t, 100)))

  # 6. Sample size < 3 should yield NA
  expect_true(is.na(tstat2logbf(2, 2)))

  # 7. Consistency check with known relationship: larger t → larger log BF
  bf_small <- tstat2logbf(2, 30)
  bf_large <- tstat2logbf(5, 30)
  expect_gt(bf_large, bf_small)
})
