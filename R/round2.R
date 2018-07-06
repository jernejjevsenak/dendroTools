#' round2
#'
#' Rounding function, where numbers are rounded up, 0.5 = 1
#'
#' @param x a numeric vector
#' @param n integer indicating the number of decimal places
#'
#' @return a numeric vector with rounded numbers
#'
#' @references
#' http://alandgraf.blogspot.com/2012/06/rounding-in-r.html
#'
#' @keywords internal

# This is a small function to round numbers up
round2 = function(x, n = 0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}
