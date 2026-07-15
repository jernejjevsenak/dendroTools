#' boot_f_pcor
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
#' @return A matrix with bootstrapped estimates of partial correlation coefficients
#'
#' @references https://www.datacamp.com/community/tutorials/bootstrap-r
#' @keywords internal
#'
boot_f_pcor <- function(data, indices, cor.type, pcor.use = "pairwise.complete.obs") {

  d <- data[indices, , drop = FALSE]

  if (sum(stats::complete.cases(d[, c("x", "y", "z")])) < 4) {
    return(NA_real_)
  }

  out <- try(
    psych::partial.r(data = d, x = c("x", "y"), y = "z",
                     use = pcor.use, method = cor.type)[2],
    silent = TRUE
  )
  if (inherits(out, "try-error")) return(NA_real_)
  as.numeric(out)
}



