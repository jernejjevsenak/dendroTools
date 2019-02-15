#' KNMI_daily_transform
#'
#' transforms daily data obtained from KNMI Climate explorer into data frame suitable for daily_response()
#'
#' @param input typical KNMI daily data format
#'
#' @return env_data suitable for daily_response()
#'
#'
#' @examples
#' data(swit272_daily_temperatures)
#' proper_env_data <- KNMI_daily_transform(swit272_daily_temperatures)

KNMI_daily_transform <- function(input){

  input$date <- ymd(input$V1)
  input$year <- year(input$date)
  input$doy <- yday(input$date)
  input$V1 <- NULL
  input$date <- NULL
  daily_matrix <- dcast(year ~ doy, data = input, value.var = "V2")
  row.names(daily_matrix) <- daily_matrix$year
  daily_matrix$year <- NULL
  daily_matrix
}
