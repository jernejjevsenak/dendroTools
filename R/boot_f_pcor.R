#' boot_f_pcor
#'
#' Generate R bootstrap replicates of a statistic applied to data.
#'
#' @param data The data as a vector, matrix or data frame. If it is a matrix or data frame then each
#' row is considered as one multivariate observation.
#' @param indices indices to be used for calculation
#' @param cor.type a character string indicating which correlation
#' @param fun a function to use in bootstrap procedure
#' coefficient is to be computed. One of "pearson" (default), "kendall", or
#' "spearman".
#'
#' @return A matrix with bootstrapped estimates of partial correlation coefficients
#'
#' @references https://www.datacamp.com/community/tutorials/bootstrap-r
#' @keywords internal
#'
boot_f_pcor <- function(data, indices, cor.type){

    partial.r(data=data[indices,], x=c("x","y"), y="z", use= "pairwise" , method = cor.type)[2]

  }



