#' Climate windows used to estimate the effective number of tests
#'
#' Internal helper for remove_insignificant_method = "sidak" in
#' daily_response(), monthly_response(), daily_response_seascorr() and
#' monthly_response_seascorr().
#'
#' The effective number of tests is obtained from the eigenvalues of the
#' correlation matrix among the screened windows, which requires complete
#' data. Earlier versions dropped every window with a missing year. After
#' dc_method = "Ar", which removes the first years of nearly every window,
#' that left only the few windows whose fitted autoregressive order was zero,
#' and the effective number of tests was underestimated, i.e. the filter was
#' too permissive. A year missing in all windows dropped every window and
#' switched the filter off.
#'
#' Now windows without a single value, which are not screened, are dropped
#' first. Years that are still missing in some window are then left out, so
#' that the estimate describes every screened window. Only if fewer than
#' min_years complete years would remain are the incomplete windows dropped
#' instead, as before. Complete data are returned unchanged.
#'
#' @param x matrix of aggregated climate windows, years as rows
#' @param x_control optional matrix of the paired control windows, as in the
#' seascorr functions. A value counts as missing when either series is missing.
#' @param min_years the smallest number of complete years that is used
#'
#' @return x, restricted to complete years and windows
#' @keywords internal
#' @noRd

sidak_window_matrix <- function(x, x_control = NULL, min_years = 10L) {

  missing <- is.na(x)

  if (!is.null(x_control)) {
    missing <- missing | is.na(x_control)
  }

  if (!any(missing)) {
    return(x)
  }

  # windows without a single value are not screened
  has_values <- colSums(!missing) > 0
  x <- x[, has_values, drop = FALSE]
  missing <- missing[, has_values, drop = FALSE]

  complete_years <- rowSums(missing) == 0

  if (sum(complete_years) >= min_years) {
    x[complete_years, , drop = FALSE]
  } else {
    x[, colSums(missing) == 0, drop = FALSE]
  }
}
