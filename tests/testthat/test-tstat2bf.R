test_that("function throws an error for non-scalar inputs", {

  # Expect an error when all inputs are non-scalar
  expect_error(
    tstat2bf(c(2, 3), c(30, 40), c(15, 20)),
    "Error: 't_statistic' must be a single numeric value.",
    fixed = TRUE
  )

  # Check that valid scalar inputs do NOT throw an error
  expect_silent(tstat2bf(2, 30))
  expect_silent(tstat2bf(2, 30, g = 15))

})
