#' smooth_matrix
#'
#' Removes unrealistic values in a matrix and replace them with an average of
#' values in a window 3 x 3 around the unrealistic value. Unrealistic value is
#' determined by a factor_drop.
#'
#' @param mat a matrix or data.frame
#' @param factor_drop a number that specifies by how many % should a value drop,
#' comparing to two the closest values in a row (i +1 and i -1), to be
#' considered as a unrealistic value.
#' @param repeats an integer that specifies number of repeats of smoothing.
#' Important when there are more unrealistic values one by another.
#'
#' @return a matrix with replaced unrealistic values
#'
#' @export
#'
#' @examples
#' library(dendroTools)
#' data(LJ_daily_temperatures)
#' data(example_proxies_1)
#' Example1 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "brnn",
#' measure = "r.squared", lower = 250, upper = 251,
#' previous_year = FALSE, brnn_smooth = TRUE, alpha = 0.1)
#' smoothed <- smooth_matrix(mat = Example1[[1]])
#'
#' mat_1 <-  matrix(seq(1.01, 2, by = 0.01)  , ncol = 10, byrow = TRUE)
#' mat_1[5 ,5] <- -1
#' mat_2 <- smooth_matrix(mat_1)

smooth_matrix <- function(mat, factor_drop = 0.7, repeats = 3){
for (r in 1:repeats){
  for (k in 1:nrow(mat)){

    for (l in 1:(ncol(mat))){

      temp_x <- mean(as.numeric(mat[k, max(1, l - 1):
          min(ncol(mat), l + 1)]), na.rm = TRUE)
      factor <- temp_x * factor_drop
      if (is.na(temp_x) == FALSE) {
        if (is.na(mat[k, l]) == FALSE) {
          if (mat[k, l] < factor) {
            mat[k, l] <- NA
            mat[k, l] <-
              mean(as.matrix(mat[max(1, k - 1):
                  min(nrow(mat), k + 1), max(1, l - 1):
                  min(ncol(mat), l + 1)]), na.rm = T)
          }
        }
      }
    }
  }
}
  return(as.matrix(mat))
}
