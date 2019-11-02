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

  capture.output(bsFit <- try(brnn(as.formula(brnn.formula), data = data[index,], neurons = neurons, tol = 1e-6), silent = TRUE))

  temporal_bsFit_predictions <- try(predict.brnn(bsFit, data[index,]), silent = TRUE)

  if (class(bsFit)[[1]] != "try-error"){

    temporal_r_squared <- 1 - (sum((data[index,][, 1] - temporal_bsFit_predictions) ^ 2) /
                                 sum((data[index,][, 1] - mean(data[index,][, 1])) ^ 2))
    temporal_adj_r_squared <- 1 - ((1 - temporal_r_squared) *
                                     ((nrow(data[index,]) - 1)) /
                                     (nrow(data[index,]) - ncol(as.data.frame(data[index,])) -  1 + 1))

  }

  c(
    temporal_r_squared,
    temporal_adj_r_squared
  )

}
