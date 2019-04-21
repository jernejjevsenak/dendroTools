#' KNMI_daily_transform
#'
#' transforms daily data obtained from KNMI Climate explorer into data frame suitable for daily_response()
#'
#' @param input typical KNMI daily data format: Data frame with two columns, first column represents date,
#' second column represents variable, such as mean temperature, precipitation, etc.
#'
#' @return env_data suitable for daily_response()
#'
#' @examples
#' data(swit272_daily_temperatures)
#' proper_env_data <- KNMI_daily_transform(swit272_daily_temperatures)

KNMI_daily_transform <- function(input){

  colnames(input) <- c("date", "variable")

  input$date <- ymd(input[,"date"])
  input$year <- year(input$date)
  input$doy <- yday(input$date)
  input$date <- NULL
  daily_matrix <- dcast(year ~ doy, data = input, value.var = "variable")
  row.names(daily_matrix) <- daily_matrix$year
  daily_matrix$year <- NULL
  daily_matrix
}
