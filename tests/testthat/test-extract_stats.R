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

test_that("extract_stats errors when no subject varies in standard score", {
  # Between-subject design with two trials per subject: no within-subject information
  df <- data.frame(id = rep(1:10, each = 2), s = rep(rep(0:1, each = 5), each = 2), v = rnorm(20))
  expect_error(extract_stats(df$s, df$v, df$id), "No subject has more than one standard score level")
})

test_that("extract_stats within-subject t matches lm with subject fixed effects", {
  df <- data.frame(subject = c("A", "A", "B", "B", "C", "C"), score = c(0, 1, 0, 1, 0, 1),
                   measurement = c(1.2, 1.8, 0.7, 1.1, 2.0, 2.4))
  t_lm <- summary(stats::lm(measurement ~ score + subject, data = df))$coefficients["score", "t value"]
  expect_equal(extract_stats(df$score, df$measurement, df$subject)$t_statistic, t_lm, tolerance = 1e-10)
})
