#' critical_r
#'
#' Calculates the critical value of a correlation coefficient for a selected
#' alpha and number of observations.
#'
#' For Pearson correlation the critical value follows from the t distribution
#' with n - 2 degrees of freedom. Spearman's rho uses the same transformation
#' applied to ranks, which is the standard approximation and is accurate for
#' the sample sizes typical of tree-ring analysis.
#'
#' Kendall's tau has a different null distribution and requires a different
#' critical value. Under the null hypothesis tau is approximately normal with
#' variance 2(2n + 5) / 9n(n - 1), the approximation used by
#' \code{\link[stats]{cor.test}} when \code{exact = FALSE}. Applying the
#' Pearson formula to tau is substantially too strict: at n = 24 it gives
#' 0.404 where the correct value is 0.286, so genuine responses are discarded.
#'
#' @param n number of observations
#' @param alpha significance level
#' @param cor_method the correlation coefficient the critical value is for.
#' One of "pearson" (the default), "spearman" or "kendall". "pearson" and
#' "spearman" share the same formula.
#'
#' @return the calculated critical value of the correlation coefficient
#'
#' @export
#'
#' @examples
#' threshold_1 <- critical_r(n = 55, alpha = 0.01)
#' threshold_2 <- critical_r(n = 55, alpha = 0.05)
#'
#' # Kendall's tau needs its own critical value
#' threshold_3 <- critical_r(n = 55, alpha = 0.05, cor_method = "kendall")

critical_r <- function(n, alpha = .05, cor_method = "pearson") {

  allowed_cor_methods <- c("pearson", "spearman", "kendall")

  if (!(cor_method %in% allowed_cor_methods)) {
    stop(paste0("cor_method is '", cor_method, "'. Instead it should be one ",
                "of: ", paste(allowed_cor_methods, collapse = ", "), "."))
  }

  if (cor_method == "kendall") {

    # Kendall (1938). Under the null, tau is approximately normal with
    # variance 2(2n + 5) / 9n(n - 1).
    tau_variance <- 2 * (2 * n + 5) / (9 * n * (n - 1))
    critical_value <- stats::qnorm(1 - alpha / 2) * sqrt(tau_variance)

  } else {

    # Pearson, and Spearman through the same transformation applied to ranks
    df <- n - 2
    critical_t <- stats::qt(alpha / 2, df, lower.tail = FALSE)
    critical_value <- sqrt((critical_t ^ 2) / ((critical_t ^ 2) + df))
  }

  return(critical_value)
}


#' p_from_r
#'
#' Two-sided p values for a matrix of correlation coefficients, using the null
#' distribution appropriate to the coefficient. Used by the false discovery
#' rate criteria, which need a p value for every screened window.
#'
#' @param r matrix or vector of correlation coefficients
#' @param n number of observations
#' @param cor_method "pearson", "spearman" or "kendall"
#' @param df_adjust degrees of freedom to subtract beyond the usual two, used
#' for partial correlations where a control variable is removed
#'
#' @return an object of the same shape as r, containing p values
#' @keywords internal
#'
p_from_r <- function(r, n, cor_method = "pearson", df_adjust = 0) {

  if (cor_method == "kendall") {

    # The normal approximation does not have a degrees-of-freedom parameter,
    # so a control variable is accounted for by reducing the effective number
    # of observations. This is an approximation; for partial correlations with
    # Kendall's tau the permutation criteria are preferable.
    n_eff <- n - df_adjust
    tau_variance <- 2 * (2 * n_eff + 5) / (9 * n_eff * (n_eff - 1))

    p <- 2 * stats::pnorm(-abs(r) / sqrt(tau_variance))

  } else {

    df <- n - 2 - df_adjust
    t_stat <- r * sqrt(df / (1 - r ^ 2))
    p <- 2 * stats::pt(-abs(t_stat), df = df)
  }

  p
}
