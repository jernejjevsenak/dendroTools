#' round_df
#'
#' Round all numeric columns in a data frame.
#'
#' @param df a data frame
#' @param digits number of digits for the round function
#'
#' @return data frame with rounded values
#'
#' @export
#'
#' @examples
#' ID <- c("a", "b", "c", "d", "e")
#' Value1 <- as.numeric(c("3.4", "6.4", "8.7", "1.1", "0.1"))
#' Value2 <- as.numeric(c("8.2", "1.7", "6.4", "1.9", "10.3"))
#' df <- data.frame(ID, Value1, Value2, stringsAsFactors = FALSE)
#' round_df(df, digits = 0)
#'
#' @references
#' https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables


round_df <- function(df, digits = 3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}
