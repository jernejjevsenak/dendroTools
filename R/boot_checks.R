#' Warning about too few bootstrap replicates for bca intervals
#'
#' Internal helper of daily_response(), monthly_response(),
#' daily_response_seascorr() and monthly_response_seascorr(). boot::boot.ci()
#' estimates the acceleration of bca intervals by a regression on the
#' bootstrap replicates, which needs at least as many replicates as resampled
#' years. With fewer replicates every bca interval is NA.
#'
#' @param n_years the number of years that are resampled
#'
#' @return NULL, invisibly
#' @keywords internal
#' @noRd

bca_boot_n_warning <- function(boot, boot_ci_type, boot_n, n_years) {

  if (isTRUE(boot) && identical(boot_ci_type, "bca") &&
      is.numeric(boot_n) && length(boot_n) == 1 && !is.na(boot_n) &&
      boot_n < n_years) {

    warning(paste0("boot_ci_type = 'bca' needs at least as many bootstrap ",
                   "replicates as analysed years, but boot_n = ", boot_n,
                   " and ", n_years, " years are analysed. boot::boot.ci() ",
                   "then cannot estimate the bca intervals, so all of them ",
                   "will be NA. Set boot_n to at least ", n_years,
                   " or use boot_ci_type = 'norm', 'perc' or 'basic'."),
            call. = FALSE)
  }

  invisible(NULL)
}
