#' @method summary ds
#' @export

summary.ds <- function(object, ...){

  # A) Extracting a matrix from a list and converting it into a data frame
  result_daily_response <- object

  result_daily_element1 <- data.frame(object[[1]])

  reference_window <- object[[15]]

  # With the following chunk, overall_maximum and overall_minimum values of
  # result_daily_element1 matrix are calculated.
  overall_max <- max(result_daily_element1, na.rm = TRUE)
  overall_min <- min(result_daily_element1, na.rm = TRUE)

  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  # There are unimportant warnings produced:
  # no non-missing arguments to max; returning -Inf
  # Based on the answer on the StackOverlow site:
  # https://stackoverflow.com/questions/24282550/no-non-missing-arguments-warning-when-using-min-or-max-in-reshape2
  # Those Warnings could be easily ignored
  if ((abs(overall_max) > abs(overall_min)) == TRUE) {

    # maximum value is located. Row indeces are needed to query information
    # about the window width used to calculate the maximum. Column name is
    # needed to query the starting day.
    max_result <- suppressWarnings(which.max(apply(result_daily_element1,
                                                   MARGIN = 2, max, na.rm = TRUE)))
    plot_column <- max_result
    plot_column_source <- plot_column
    max_index <- which.max(result_daily_element1[, names(max_result)])
    row_index <- row.names(result_daily_element1)[max_index]
    temporal_vector <- unlist(result_daily_element1[max_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_metric <- round(max(temporal_vector, na.rm = TRUE), 3)

    # Here we remove missing values at the end of the temporal_vector.
    # It is important to remove missing values only at the end of the
    # temporal_vector!
    row_count <- nrow(temporal_vector)
    delete_rows <- 0
    while (is.na(temporal_vector[row_count, ] == TRUE)){
      delete_rows <- delete_rows + 1
      row_count <-  row_count - 1
    }
    # To check if the last row is a missing value
    if (is.na(temporal_vector[nrow(temporal_vector), ] == TRUE)) {
      temporal_vector <-  temporal_vector[-c((row_count + 1):(row_count +
                                                                delete_rows)), ]
    }
    temporal_vector <- data.frame(temporal_vector)
  }

  if ((abs(overall_max) < abs(overall_min)) == TRUE) {

    # minimum value is located. Row indeces are needed to query information
    # about the window width used to calculate the minimum. Column name is
    # needed to query the starting day.
    min_result <- suppressWarnings(which.min(apply(result_daily_element1,
                                                   MARGIN = 2, min, na.rm = TRUE)))
    plot_column <- min_result
    plot_column_source <- plot_column
    min_index <- which.min(result_daily_element1[, names(min_result)])
    row_index <- row.names(result_daily_element1)[min_index]
    temporal_vector <- unlist(result_daily_element1[min_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_metric <- round(min(temporal_vector, na.rm = TRUE), 3)

    # Here we remove missing values
    # We remove missing values at the end of the temporal_vector.
    # It is important to remove missing values only at the end of the
    # temporal_vector!

    row_count <- nrow(temporal_vector)
    delete_rows <- 0
    while (is.na(temporal_vector[row_count, ] == TRUE)){
      delete_rows <- delete_rows + 1
      row_count <-  row_count - 1
    }
    # To check if the last row is a missing value
    if (is.na(temporal_vector[nrow(temporal_vector), ] == TRUE)) {
      temporal_vector <-  temporal_vector[-c((row_count + 1):(row_count +
                                                                delete_rows)), ]
    }
    temporal_vector <- data.frame(temporal_vector)
  }



  # In case of previous_year == TRUE, we calculate the day of a year
  # (plot_column), considering 366 days of previous year.
  if (nrow(temporal_vector) > 366 & plot_column > 366) {
    plot_column_extra <- plot_column %% 366
  } else {
    plot_column_extra <- plot_column
  }


  # Here we define a data frame of dates and corresponing day of year (doi). Later
  # this dataframe will be used to describe tht optimal sequence od days
  doy <- seq(1:730)
  date <- seq(as.Date('2013-01-01'),as.Date('2014-12-31'), by = "+1 day")
  # date[366] <- as.Date('2015-12-31')
  date <- format(date, "%b %d")
  date_codes <- data.frame(doy = doy, date = date)


  # Here, there is a special check if optimal window width is divisible by 2 or not.
  if (as.numeric(row_index)%%2 == 0){
    adjustment_1 = 0
    adjustment_2 = 1
  } else {
    adjustment_1 = 1
    adjustment_2 = 2
  }



  if (reference_window == "start"){
    Optimal_string <- paste("\nOptimal Selection:",
                            as.character(date_codes[plot_column_extra, 2]),"-",
                            as.character(date_codes[plot_column_extra + as.numeric(row_index) - 1, 2]))
  } else if (reference_window == "end") {
    Optimal_string <- paste("\nOptimal Selection:",
                            as.character(date_codes[plot_column_extra - as.numeric(row_index) + 1, 2]),"-",
                            as.character(date_codes[plot_column_extra, 2]))
  } else if (reference_window == "middle") {
    Optimal_string <- paste("\nOptimal Selection:",
                            as.character(date_codes[(round2((plot_column_extra - as.numeric(row_index)/2)) - adjustment_1), 2]),"-",
                            as.character(date_codes[(round2((plot_column_extra + as.numeric(row_index)/2)) - adjustment_2), 2]))
  }




  # Here we define titles. They differ importantly among methods and arguments
  # in the final output list from daily_response() function
  if (result_daily_response[[2]] == "cor"){
    y_lab <- "Correlation Coefficient"
  } else if (result_daily_response[[2]] == "pcor"){
    y_lab <- "Partial Correlation Coefficient"
  } else if (result_daily_response[[3]] == "r.squared"){
    y_lab <- "Explained Variance"
  } else if (result_daily_response[[3]] == "adj.r.squared"){
    y_lab <- "Adjusted Explained Variance"
  }

  if (nrow(temporal_vector) > 366){
    x_lab <- "Day of Year  (Including Previous Year)"
  } else if (nrow(temporal_vector) <= 366){
    x_lab <- "Day of Year"
  }

  if (reference_window == 'start' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("\nStarting Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}

  if (reference_window == 'start' &&  plot_column <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("\nStarting Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}

  if (reference_window == 'start' &&  plot_column <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("\nStarting Day of Optimal Window Width: Day ",
                               plot_column_extra)}


  if (reference_window == 'end' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("\nEnding Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}

  if (reference_window == 'end' &&  plot_column  <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("\nEnding Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}

  if (reference_window == 'end' &&  plot_column  <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("\nEnding Day of Optimal Window Width: Day ",
                               plot_column_extra)}


  if (reference_window == 'middle' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("\nMiddle Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}

  if (reference_window == 'middle' &&  plot_column  <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("\nMiddle Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}

  if (reference_window == 'middle' &&  plot_column  <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("\nMiddle Day of Optimal Window Width: Day ",
                               plot_column_extra)}

  optimal_window_string <- paste0("\nOptimal Window Width: ", as.numeric(row_index),
                                  " Days")

  optimal_calculation <- paste0("\nThe Highest ", y_lab,": " , calculated_metric)

  period_string <- paste0("Analysed Period: ", result_daily_response[[4]])

  if (result_daily_response[[2]] == 'cor'){
    method_string <- paste0("\nMethod: Pearson Correlation")
  } else if (result_daily_response[[2]] == 'pcor'){
    method_string <- paste0("\nMethod: Partial Pearson Correlation")
  } else if (result_daily_response[[2]] == 'lm'){
    method_string <- paste0("\nMethod: Linear Regression")
  } else if (result_daily_response[[2]] == 'brnn'){
    method_string <- paste0("\nMethod: ANN with Bayesian Regularization")
  }

  cat("The highest calculated metric", calculated_metric)

}

