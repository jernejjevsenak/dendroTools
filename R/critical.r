#' critical_r
#'
#' Calculates critical value of Pearson correlation coefficient
#' for a selected alpha.
#'
#' @param n number of observations
#' @param alpha significance level
#'
#' @return calculated critical value of Pearson correlation coefficient
#'
#' @export
#'
#' @examples
#' threshold_1 <- critical_r(n = 55, alpha = 0.01)
#' threshold_2 <- critical_r(n = 55, alpha = 0.05)

critical_r <- function(n, alpha = .05) {
  df <- n - 2
  critical_t <- qt(alpha / 2, df, lower.tail = FALSE)
  critical_r_value <- sqrt((critical_t ^ 2) / ((critical_t ^ 2) + df))
  return(critical_r_value)
}
