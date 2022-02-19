#' data_transform
#'
#' Transforms daily data with two columns (date and variable) into data frame suitable for daily or
#' monthly analysis with dendroTools.
#'
#' @param input typical daily data format: Data frame with two columns, first column represents date,
#' second column represents variable, such as mean temperature, precipitation, etc. Date should be in
#' format Year-Month-Day (e.g. "2019-05-15")
#' @param format character string indicating the desired output format. Should be "daily" or "monthly".
#' Daily format returns a data frame with 366 columns (days), while monthly format returns data frame
#' with 12 columns (months). Years are indicated as row names.
#' @param monthly_aggregate_function character string indicating, how to aggregate daily into monthly
#' data. It can be "mean" or "sum". Third option is "auto" (default). In this case function will try
#' to guess whether input is temperature or precipitation data. For temperature, it will use "mean",
#' for precipitation "sum".
#' @param date_format Describe the format of date. It should be one of "ymd", "ydm", "myd", "mdy",
#' "dmy", "dym".
#'
#' @return env_data suitable for daily or monthly analysis with dendroTools.
#'
#' @examples
#' data(swit272_daily_temperatures)
#' proper_daily_data <- data_transform(swit272_daily_temperatures, format = "daily",
#'    date_format = "ymd")
#'
#' proper_monthly_data <- data_transform(swit272_daily_temperatures, format = "monthly",
#'    date_format = "ymd")
#'
#' data(swit272_daily_precipitation)
#' proper_daily_data <- data_transform(swit272_daily_precipitation, format = "daily",
#'    date_format = "ymd")
#'
#' proper_monthly_data <- data_transform(swit272_daily_precipitation, format = "monthly",
#'    date_format = "ymd")

data_transform <- function(input, format = "daily",
                           monthly_aggregate_function = "auto",
                            date_format = "ymd"){

  colnames(input) <- c("date", "variable")

  # This is needed to ensure the proper functionality
  input <- data.frame(input)

  if (date_format == "ymd"){
    input$date <- ymd(input[,"date"])
  } else if ( date_format == "ydm"){
    input$date <- ydm(input[,"date"])
  } else if (date_format == "mdy"){
    input$date <- mdy(input[,"date"])
  } else if ( date_format == "myd"){
    input$date <- myd(input[,"date"])
  } else if ( date_format == "dmy"){
    input$date <- dmy(input[,"date"])
  } else if ( date_format == "dym"){
    input$date <- dym(input[,"date"])
  } else {
    stop(paste("The argument date_format is not in one of ydm, ymd, mdy, myd, dym, dmy!"))
  }

    input$year <- year(input$date)

  if (monthly_aggregate_function == "auto"){
    share_na <- sum(input$variable == 0, na.rm = TRUE)/length(input$variable)

    if (share_na > 0.1){

      monthly_aggregate_function <- "sum"

    } else {

      monthly_aggregate_function <- "mean"
    }
  }

  if (format == "daily"){

    input$doy <- yday(input$date)
    input$date <- NULL
    daily_matrix <- dcast(year ~ doy, data = input, value.var = "variable", fill = NA)
    row.names(daily_matrix) <- daily_matrix$year
    daily_matrix$year <- NULL

  } else if (format == "monthly") {

    input$month <- month(input$date)
    input$date <- NULL
    colnames(input)[1] <- "value"


    if (monthly_aggregate_function == "mean"){
      daily_matrix <- dcast(year ~ month, data = input, value.var = "value", fun.aggregate = mean, na.rm = TRUE)

    } else if (monthly_aggregate_function == "sum") {

      # here we create user specified sum, since Primitive sum returns 0 when only NA values are present
      daily_matrix <- dcast(year ~ month, data = input, value.var = "value", fun.aggregate = function(x)

        if (length(x) < 1) NA_real_ else if (sum(!is.na(x)) < 1) NA_real_ else sum(x, na.rm = TRUE))

    } else {
      stop(paste0("monthly_aggregate_function argument should be 'mean' or 'sum'! It is ", monthly_aggregate_function))
    }

    row.names(daily_matrix) <- daily_matrix$year
    daily_matrix$year <- NULL

  } else {
    stop(paste0("format argument should be 'daily' or 'monthly'! It is ", format))
    }

  daily_matrix
}
