#' count_ones
#'
#' calculates share of integer 1 in a vector
#'
#' @param vector a vector of integers
#'
#' @return an integer of counted ones
#'
#'
#' @examples
#' vector_1 <- seq(1:10)
#' # count_ones(vector_1)
#' @keywords internal

count_ones <- function(vector){
  sum(vector == 1)
}
