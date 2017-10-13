#' count_ones
#'
#' calculates share of intiger 1 in a vector
#'
#' @param vector a vector of intigers
#'
#' @return an intiger of counted ones
#'
#' @export
#'
#' @examples
#' vector_1 <- seq(1:10)
#' count_ones(vector_1)

count_ones <- function(vector){
  sum(vector == 1)
}
