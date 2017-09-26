#' years_to_rownames
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
#' \dontrun{
#' new_df = years_to_rownames(data = daily_sequences, column_year = "years")
#' }
years_to_rownames <- function(data, column_year) {
  year_index <- grep(column_year, colnames(data))
  row.names(data) <- data[, year_index]
  data <- data[, -year_index]
}
