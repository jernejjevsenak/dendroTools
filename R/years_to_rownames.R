#'
#'
#' Function returns a data frame with row names as years
#'
#' @param data a data frame to be manipulated
#' @param column_year string specifying a column with years
#'
#' @return a data frame with years as row names
#'
#' @export
#'
#' @examples
#' data <- data.frame(years = seq(1950, 2015), observations = rnorm(66))
#' new_data <- years_to_rownames(data = data, column_year = "years")
#'
#' data <- data.frame(observations1 = rnorm(66), years = seq(1950, 2015),
#' observations2 = rnorm(66), observations3 = rnorm(66))
#' new_data <- years_to_rownames(data = data, column_year = "years")

years_to_rownames <- function(data, column_year) {
  data <- data.frame(data) # Data needs to of class data.frame!

  year_index <- grep(column_year, colnames(data))

  names <- colnames(data)
  names <- names[-year_index]

  row.names(data) <- data[, year_index]
  data <- as.data.frame(data[, -as.numeric(year_index), F])

  colnames(data) <- names

  return(data)
}
