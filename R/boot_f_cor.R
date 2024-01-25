#' boot_f_cor
#'
#' Generate R bootstrap replicates of a statistic applied to data.
#'
#' @param data A vector, matrix, or data frame containing the input data.
#' If a matrix or data frame is provided, each row is treated as a separate
#' multivariate observation.
#' @param indices A numeric or integer vector specifying the indices to be
#' utilized in the calculation. These indices determine which elements or rows
#' of the 'data' are to be included in the analysis
#' @param cor.type A string specifying the method of correlation to be applied
#' in the bootstrap procedure. Available options are "pearson" (default),
#' "kendall", or "spearman".
#'
#' @return A matrix with bootstrapped estimates of correlation coefficients
#'
#' @references https://www.datacamp.com/community/tutorials/bootstrap-r
#' @keywords internal
#'
boot_f_cor <- function(data, indices, cor.type){
  dt <- data[indices,]
    c(
      cor(dt[,1], dt[,2], method = cor.type),
      median(dt[,1]),
      median(dt[,2])
    )

  }
