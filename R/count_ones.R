#' count_ones
#'
#' calculates share of intiger 1 in a vector
#'
#' @param vector a vector of intigers
#'
#' @return an intiger of counted ones

count_ones <- function(vector){
  sum(vector == 1)
}
