#' boot_f_lm
#'
#' Generate R bootstrap replicates of a lm function applied to data.
#'
#' @param data data frame with variables for model fitting
#' @param index indices to be used for calculation
#' @param lm.formula object of class formula to be passed into lm function
#'
#' @return A matrix with bootstrapped estimates of r squared and adjusted r squared
#'
#' @references https://www.datacamp.com/community/tutorials/bootstrap-r
#' @keywords internal
#'

boot_f_lm <- function(data, index, lm.formula) {
  d <- data[index, , drop = FALSE]

  # a resample can lose complete cases; guard like boot_f_cor
  if (sum(stats::complete.cases(d)) < 4) {
    return(c(NA_real_, NA_real_))
  }

  bsFit <- try(lm(lm.formula, data = d), silent = TRUE)
  if (inherits(bsFit, "try-error")) {
    return(c(NA_real_, NA_real_))
  }

  s <- summary(bsFit)
  c(s$r.squared, s$adj.r.squared)
}
