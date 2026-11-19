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

# Three datasets, each containing two measurements of different quality (B better than A)
set.seed(1)
make_ds <- function(ds) {
  s <- rep(c(0, 1), each = 20)
  rbind(data.frame(dataset_col = ds, measurement_name_col = "A", standard_scores_col = s,
                   measurement_values_col = rnorm(40, 0.1 * s)),
        data.frame(dataset_col = ds, measurement_name_col = "B", standard_scores_col = s,
                   measurement_values_col = rnorm(40, 1.5 * s)))
}
full_overlap <- rbind(make_ds("D1"), make_ds("D2"), make_ds("D3"))
run_strict <- function(d) {
  m <- strict_meta_bf_from_data(d, "dataset_col", "measurement_name_col",
                                "standard_scores_col", "measurement_values_col")$log_bf_matrix
  m[order(rownames(m)), order(colnames(m))]
}

test_that("strict_meta_bf_from_data does not depend on the order of rows or measurements", {
  ref <- run_strict(full_overlap)
  reversed <- full_overlap[order(full_overlap$measurement_name_col, decreasing = TRUE), ]
  expect_equal(run_strict(reversed), ref, tolerance = 1e-12)
  expect_equal(run_strict(full_overlap[sample(nrow(full_overlap)), ]), ref, tolerance = 1e-12)
  expect_gt(ref["B", "A"], 0)                    # better measurement favoured
  expect_equal(ref["A", "B"], -ref["B", "A"])    # antisymmetric
})

test_that("strict_meta_bf_from_data equals calibration_summary differences under full overlap", {
  ref <- run_strict(full_overlap)
  cs <- suppressWarnings(calibration_summary(full_overlap, "dataset_col", "measurement_name_col",
                                             "measurement_values_col", "standard_scores_col"))
  expect_equal(ref["B", "A"], cs$log_bayes_factor_diffs["B", "A"], tolerance = 1e-10)
})
