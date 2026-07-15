#' boot_f_brnn
#'
#' Generate R bootstrap replicates of a brnn function applied to data. Function is based on boot() from
#' boot R package
#'
#' @param data data frame with variables for model fitting
#' @param index indices to be used for calculation
#' @param brnn.formula object of class formula to be passed into brnn function
#'
#' @return A matrix with bootstrapped estimates of r squared and adjusted r squared
#'
#' @references https://www.datacamp.com/community/tutorials/bootstrap-r
#' @keywords internal
#'

boot_f_brnn <- function(data, index, brnn.formula, neurons = 3) {

  d <- data[index, , drop = FALSE]

  if (sum(stats::complete.cases(d)) < 4) {
    return(c(NA_real_, NA_real_))
  }

  capture.output(
    bsFit <- try(brnn(as.formula(brnn.formula), data = d,
                      neurons = neurons, tol = 1e-6), silent = TRUE)
  )
  if (inherits(bsFit, "try-error")) {
    return(c(NA_real_, NA_real_))
  }

  preds <- try(predict.brnn(bsFit, d), silent = TRUE)
  if (inherits(preds, "try-error")) {
    return(c(NA_real_, NA_real_))
  }

  temporal_r_squared <- 1 - (sum((d[, 1] - preds)^2) /
                               sum((d[, 1] - mean(d[, 1]))^2))
  temporal_adj_r_squared <- 1 - ((1 - temporal_r_squared) *
                                   (nrow(d) - 1) /
                                   (nrow(d) - ncol(as.data.frame(d)) - 1 + 1))

  c(temporal_r_squared, temporal_adj_r_squared)
}
