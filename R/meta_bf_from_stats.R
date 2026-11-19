#' Compute Meta-Analytic Log Bayes Factor From Stats
#'
#' This function computes a **meta-analytic log Bayes Factor** (natural log) by combining
#' test statistics from multiple datasets using the weighted approach of
#' Nikolakopoulos & Ntzoufras (2021). You may supply per-dataset \(t\)-statistics directly,
#' or provide Bayes Factors to be transformed into \(t\)-statistics where needed. Since
#' Bayes Factors carry no information on the sign of the effect, the effect direction must
#' then be supplied via \code{effect_directions}. The combined statistic is then converted
#' into a log Bayes Factor.
#'
#' To keep outputs uniform with \code{meta_bf_from_data}, the return value includes:
#' \itemize{
#'   \item \code{meta_results}: a 1-row data frame with columns
#'         \code{measurement}, \code{log_combined_bayes_factor} (natural log),
#'         \code{retrodiction_score}, \code{sample_size} (sum of per-dataset \(n\)),
#'         and \code{combined_t_statistic}.
#'   \item \code{subject_overlap}: always \code{NULL} (not applicable when only stats are provided).
#'   \item \code{weights}: a named list (single element, named by \code{measurement_name}) containing
#'         a per-dataset table with columns \code{dataset}, \code{weight}, \code{t_statistic},
#'         \code{sample_size}, \code{df}, \code{r_i}, \code{sum_squares}.
#' }
#'
#' @param n_datasets Integer; number of datasets.
#' @param t_statistics Numeric vector of per-dataset (signed) \(t\)-statistics. Defaults to
#'   \code{rep(NA, n_datasets)}. Where an entry is \code{NA}, the corresponding \code{bayes_factors}
#'   and \code{effect_directions} values must be provided.
#' @param sample_sizes Numeric vector of per-dataset sample sizes (\(n\)).
#' @param sum_squares Numeric vector of per-dataset sums of squares used by \code{compute_weights()}.
#'   Defaults to \code{rep(NA, n_datasets)}; if any entry is \code{NA}, weights are based on
#'   sample sizes only.
#' @param bayes_factors Numeric vector of per-dataset Bayes Factors \eqn{BF_{10}} against the null
#'   model, on the raw (not log) scale and > 0. Defaults to \code{rep(NA, n_datasets)}.
#'   Ignored where \code{t_statistics} is non-missing; otherwise each BF is converted to the
#'   magnitude of a \(t\)-statistic via \code{logbf2tstat(log(bf), n, g = n)}.
#' @param effect_directions Numeric vector of effect directions (\code{1} for positive, \code{-1}
#'   for negative). Defaults to \code{rep(NA, n_datasets)}. Required wherever a \(t\)-statistic is
#'   derived from a Bayes Factor, otherwise an error is thrown; ignored where \code{t_statistics}
#'   is supplied (the sign of the \(t\)-statistic is used).
#' @param measurement_name Optional single string used to label the output row and the weights list name.
#'   Defaults to \code{"overall"}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{meta_results}{Data frame with one row and columns:
#'     \code{measurement}, \code{log_combined_bayes_factor} (natural log),
#'     \code{retrodiction_score}, \code{sample_size}, \code{combined_t_statistic}.}
#'   \item{subject_overlap}{Always \code{NULL}.}
#'   \item{weights}{Named list of length 1 (name = \code{measurement_name}) containing a per-dataset
#'     data frame with columns \code{dataset}, \code{weight}, \code{t_statistic}, \code{sample_size},
#'     \code{df}, \code{r_i}, \code{sum_squares}.}
#' }
#'
#' @details
#' The \code{retrodiction_score} is computed exactly as in \code{meta_bf_from_data()}: each
#' \(t\)-statistic is converted to a correlation, Fisher-z transformed, averaged with weights
#' equal to the residual degrees of freedom minus one, and back-transformed.
#'
#' @seealso \code{\link{meta_bf_from_data}}, \code{\link{compute_weights}},
#'   \code{\link{tstat2logbf}}, \code{\link{logbf2tstat}}
#' @export
#'
#' @examples
#' # Example 1: Compute from t-statistics
#' n_datasets   <- 3
#' t_statistics <- c(2.5, 1.8, 2.0)
#' sample_sizes <- c(30, 40, 35)
#' sum_squares  <- c(5, 6, 5.5)
#' out1 <- meta_bf_from_stats(
#'     n_datasets,
#'     t_statistics,
#'     sample_sizes,
#'     sum_squares,
#'     measurement_name = "A"
#' )
#' out1$meta_results
#' out1$weights$A
#'
#' # Example 2: Compute using Bayes Factors (converted to t where t is NA);
#' # the effect directions must be supplied, since Bayes Factors carry no sign
#' bayes_factors <- c(10, 5, 8)
#' out2 <- meta_bf_from_stats(n_datasets, t_statistics = c(NA, NA, NA),
#'                            sample_sizes = sample_sizes, sum_squares = sum_squares,
#'                            bayes_factors = bayes_factors, effect_directions = c(1, 1, 1),
#'                            measurement_name = "overall")
#' out2$meta_results
#'
meta_bf_from_stats <- function(
    n_datasets,
    t_statistics      = rep(NA, n_datasets),
    sample_sizes,
    sum_squares       = rep(NA, n_datasets),
    bayes_factors     = rep(NA, n_datasets),
    effect_directions = rep(NA, n_datasets),
    measurement_name  = "overall"
) {

  # --- Input checks -----------------------------------------------------------
  # lengths consistent with n_datasets
  lens <- c(length(t_statistics), length(sample_sizes), length(sum_squares),
            length(bayes_factors), length(effect_directions))
  if (any(lens != n_datasets)) {
    stop("Error: All input vectors (t_statistics, sample_sizes, sum_squares, bayes_factors) must have length n_datasets.")
  }
  if (!is.character(measurement_name) || length(measurement_name) != 1) {
    stop("Error: measurement_name must be a single string.")
  }

  # --- Fill missing t's from BFs: a BF carries no sign, so the direction must be supplied
  for (i in seq_along(t_statistics)) {
    if (is.na(t_statistics[i]) && is.na(bayes_factors[i])) {
      stop(paste("Error: Both t_statistic and bayes_factor are missing for index", i))
    }
    if (is.na(t_statistics[i])) {
      if (bayes_factors[i] <= 0) {
        stop(sprintf("Error: bayes_factors[%d] must be > 0.", i))
      }
      if (!(effect_directions[i] %in% c(-1, 1))) {
        stop(sprintf(paste0("Error: effect_directions[%d] must be 1 or -1 when the t-statistic ",
                            "is derived from a Bayes factor (Bayes factors carry no sign)."), i))
      }
      t_statistics[i] <- effect_directions[i] *
        logbf2tstat(log(bayes_factors[i]), sample_sizes[i], g = sample_sizes[i])
    }
  }

  # --- Weights and combined statistics (natural log) --------------------------
  weights <- compute_weights(t_statistics, sample_sizes, sum_squares)
  combined_t_statistic <- sum(weights * t_statistics)

  total_n <- sum(sample_sizes)
  log_combined_bayes_factor <- tstat2logbf(combined_t_statistic, total_n, g = total_n)

  # --- Meta-analytic retrodiction score using Fisher-z (as in meta_bf_from_data) ---
  df_i  <- sample_sizes - 2
  r_i   <- t_statistics / sqrt(t_statistics^2 + df_i)
  z_i   <- atanh(pmax(pmin(r_i, 0.999999), -0.999999))
  w_i   <- pmax(df_i - 1, 1)  # effective Fisher-z weights from residual df
  retrodiction_score <- tanh(sum(w_i * z_i) / sum(w_i))

  # --- Assemble outputs to mirror meta_bf_from_data --------------------------
  weights_df <- data.frame(
    dataset      = seq_len(n_datasets),
    weight       = as.numeric(weights),
    t_statistic  = as.numeric(t_statistics),
    sample_size  = as.numeric(sample_sizes),
    df           = as.numeric(df_i),
    r_i          = as.numeric(r_i),
    sum_squares  = as.numeric(sum_squares),
    stringsAsFactors = FALSE
  )

  meta_results <- data.frame(
    measurement               = measurement_name,
    log_combined_bayes_factor = log_combined_bayes_factor,
    retrodiction_score        = retrodiction_score,
    sample_size               = total_n,
    combined_t_statistic      = combined_t_statistic,
    stringsAsFactors = FALSE
  )

  list(
    meta_results    = meta_results,
    subject_overlap = NULL,
    weights         = stats::setNames(list(weights_df), measurement_name)
  )
}
