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
  bsFit <- lm(lm.formula, data=data[index,])
  c(
    summary(bsFit)$r.squared,
    summary(bsFit)$adj.r.squared
  )
}
