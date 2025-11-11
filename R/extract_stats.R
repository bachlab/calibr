#' Extract Statistical Information for Bayes Factor Computation
#'
#' This function extracts the key statistics needed for computing a correlation-based
#' Bayes factor (Eq. 5). It computes the Pearson correlation between `measurement_values`
#' and `standard_scores`, and converts it to a t-statistic. When `subject_ids` are provided,
#' it computes the **within-subject partial correlation** (i.e., the correlation of
#' subject-centered residuals) to account for repeated or paired measurements.
#'
#' @param standard_scores A numeric (typically binary 0/1) vector representing the
#'        standard or outcome scores.
#' @param measurement_values A numeric vector of measurement values to be correlated
#'        with `standard_scores`.
#' @param subject_ids (Optional) A vector of subject identifiers. If provided, the
#'        correlation is computed within subjects by partialling out subject effects
#'        (equivalent to including subject fixed effects in a regression model).
#'
#' @return A list containing:
#'   \describe{
#'     \item{t_statistic}{The t-statistic corresponding to the (partial) correlation
#'           between `measurement_values` and `standard_scores`. This statistic is
#'           directly compatible with Eq. 5 for Bayes factor computation.}
#'     \item{sum_squares}{The sum of squares of `measurement_values` at the appropriate
#'           level (within-subject if repeated, across all data if unpaired), used for weighting.}
#'     \item{sample_size}{An effective sample size such that \code{sample_size - 2}
#'           equals the residual degrees of freedom associated with the t-statistic.}
#'   }
#'
#' @details
#' * If `subject_ids` is **not provided**, the function computes a standard Pearson correlation
#'   and converts it to a t-statistic using \eqn{t = r\sqrt{(n-2)/(1-r^2)}} with \eqn{df = n-2}.
#' * If `subject_ids` is **provided**, both variables are centered within each subject,
#'   yielding a within-subject (partial) correlation. The corresponding t-statistic is computed
#'   using \eqn{df = N - S - 1}, where \eqn{N} is the number of rows and \eqn{S} is the number
#'   of unique subjects. This t is identical to that from \code{lm(measurement ~ score + subject)}.
#'
#' The resulting \code{t_statistic} and \code{sample_size} are suitable for direct use in
#' \code{tstat2bf()} and then \code{meta_bf_from_data()}.
#'
#' @examples
#' # Example 1: Unpaired design (between subjects)
#' df1 <- data.frame(
#'   subject = paste0("S", 1:4),
#'   score = c(0, 0, 1, 1),
#'   measurement = c(1.0, 0.6, 1.9, 1.5)
#' )
#' extract_stats(df1$score, df1$measurement)
#'
#' # Example 2: Paired / repeated design (within subjects)
#' df2 <- data.frame(
#'   subject = c("A","A","B","B","C","C"),
#'   score = c(0,1,0,1,0,1),
#'   measurement = c(1.2,1.8,0.7,1.1,2.0,2.4)
#' )
#' extract_stats(df2$score, df2$measurement, df2$subject)
#'
#' @seealso \code{\link{tstat2bf}}, \code{\link{meta_bf_from_data}}
#' @export
extract_stats <- function(standard_scores, measurement_values, subject_ids = NULL) {
  # --- checks
  if (anyNA(standard_scores) || anyNA(measurement_values)) {
    stop("Error: 'standard_scores' or 'measurement_values' contain NA values.")
  }
  if (length(unique(standard_scores)) < 2) {
    stop("Error: 'standard_scores' must contain at least two unique values.")
  }

  df <- data.frame(
    score = as.numeric(standard_scores),
    meas  = as.numeric(measurement_values),
    subj  = if (is.null(subject_ids)) factor(seq_along(standard_scores)) else factor(subject_ids)
  )

  # design diagnostics
  tab <- table(df$subj)
  has_repeats <- any(tab > 1)

  if (!has_repeats) {
    # ---------------- Unpaired: plain correlation ----------------
    n <- nrow(df)
    if (n < 3) stop("Not enough observations (need at least 3).")

    r <- stats::cor(df$meas, df$score)
    df_t <- n - 2
    t_statistic <- r * sqrt(df_t / (1 - r^2))
    sample_size <- n                                   # so (n - 2) = df_t
    sum_squares <- sum((df$meas - mean(df$meas))^2)    # SS at independent (subject) level

  } else {
    # ---------------- Repeated/paired: partial correlation via FWL ----------------
    # Residualize BOTH meas and score on subject (within-subject centering)
    meas_ws  <- df$meas  - ave(df$meas,  df$subj, FUN = mean)
    score_ws <- df$score - ave(df$score, df$subj, FUN = mean)

    # If every subject has both 0 and 1 exactly once, this equals the paired t case.
    # More generally, it matches lm(meas ~ score + subj) by FWL.
    N  <- nrow(df)
    S  <- nlevels(df$subj)
    df_t <- N - S - 1
    if (df_t < 1) stop("Not enough residual degrees of freedom after removing subject effects.")

    r <- stats::cor(meas_ws, score_ws)
    t_statistic <- r * sqrt(df_t / (1 - r^2))

    # Make BF's (n - 2) equal df_t  -> n = df_t + 2
    sample_size <- df_t + 2

    # SS at the independent (within-subject) level
    sum_squares <- sum(meas_ws^2)
  }

  list(
    t_statistic = as.numeric(t_statistic),
    sum_squares = as.numeric(sum_squares),
    sample_size = as.numeric(sample_size)
  )
}






