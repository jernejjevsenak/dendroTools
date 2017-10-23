#' plot_extreme
#'
#' Graphs a line plot of a row with the highest measure in a matrix, produced by
#' \code{\link{daily_response}} function.
#'
#' @param result_daily_response a list with three objects as produced by
#' daily_response function
#' @param title logical, if set to FALSE, no plot title is displayed
#'
#' @return A ggplot2 object containing the plot display
#' @export
#'
#' @examples
#' \dontrun{
#' data(LJ_daily_temperatures)
#' data(example_proxies_1)
#' Example1 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "lm", measure = "r.squared",
#' fixed_width = 90, previous_year = TRUE)
#' plot_extreme(Example1)
#'
#' Example2 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "brnn",
#' measure = "adj.r.squared", lower_limit = 50, upper_limit = 55, neurons = 1,
#' row_names_subset = TRUE)
#' plot_extreme(Example2)
#'
#' # Example with negative correlations
#' data(example_proxies_2)
#' LJ_daily_temperatures_subset = LJ_daily_temperatures[-c(53:55), ]
#' Example3 <- daily_response(response = example_proxies_2,
#' env_data = LJ_daily_temperatures_subset, method = "cor",
#' lower_limit = 35, upper_limit = 40)
#' plot_extreme(Example3)
#' }

plot_extreme <- function(result_daily_response, title = TRUE) {

  # Short description of the function. It
  # - extracts matrix (the frst object of a list)
  # - in case of method == "cor" (second object of a list), calculates the
  # highest minimum and maximum and compare its absolute values. If absolute
  # minimum is higher than maximum, we have to plot minimum correlations.
  # - query the information about windows width (row names of the matrix) and
  # starting day of the mighest (absolute) measure (column names of the matrix).
  # - draws a ggplot line plot.

  # A) Extracting a matrix from a list and converting it into a data frame
  result_daily_element1 <- data.frame(result_daily_response[[1]])

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
    max_index <- which.max(result_daily_element1[, names(max_result)])
    row_index <- row.names(result_daily_element1)[max_index]
    temporal_vector <- unlist(result_daily_element1[max_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_measure <- round(max(temporal_vector, na.rm = TRUE), 3)

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
      temporal_vector <-  temporal_vector[-c(row_count:(row_count +
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
    min_index <- which.min(result_daily_element1[, names(min_result)])
    row_index <- row.names(result_daily_element1)[min_index]
    temporal_vector <- unlist(result_daily_element1[min_index, ])
    temporal_vector <- temporal_vector[!is.na(temporal_vector)]
    temporal_vector <- data.frame(temporal_vector)
    calculated_measure <- round(min(temporal_vector, na.rm = TRUE), 3)

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
      temporal_vector <-  temporal_vector[-c(row_count:(row_count +
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

  # The final plot is being created. The first part of a plot is the same,
  # the second part is different, depending on temporal.vector, plot_column,
  # method and measure string stored in result_daily_response. The second part
  # defines xlabs, xlabs and ggtitles.

  # The definition of theme
  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18), text = element_text(size = 18),
          plot.title = element_text(size = 16,  face = "bold"))

  if (title == FALSE){
    journal_theme <- journal_theme +
      theme(plot.title = element_blank())
  }

  # in the next chunk, warnings are supressed. At the end of the vector,
  # there are always missing values, which are a result of changing window
  # width calclulations. Those warnings are not important and do not affect
  # our results at all
  final_plot <- suppressWarnings(
  ggplot(temporal_vector, aes(y = temporal_vector,
    x = seq(1, length(temporal_vector)))) + geom_line(lwd = 1.2) +
     geom_vline(xintercept = plot_column, col = "red") +
     scale_x_continuous(breaks = sort(c(seq(0, nrow(temporal_vector), 50)), decreasing = FALSE),
       labels = sort(c(seq(0, nrow(temporal_vector), 50)))) +
       annotate("label", label = as.character(calculated_measure),
         y = calculated_measure, x = plot_column + 15) +
    annotate("label", label = paste("Day", as.character(plot_column), sep = " "),
             y = min(temporal_vector, na.rm = TRUE) + 0.2*min(temporal_vector, na.rm = TRUE), x = plot_column + 15) +
    journal_theme)

  # If previous_year = TRUE, we add a vertical line with labels of
  # previous and current years
  if (ncol(result_daily_element1) > 366) {
    final_plot <- final_plot +
      annotate(fontface = "bold", label = 'Previous Year', geom = 'label',
               x = 366 - ncol(result_daily_element1) / 12.8,
               y = calculated_measure - (calculated_measure/5)) +
      annotate(fontface = "bold", label = 'Current Year', geom = 'label',
               x = 366 + ncol(result_daily_element1) / 13.5,
               y = calculated_measure -(calculated_measure/5)) +
      geom_vline(xintercept = 366, size = 1)
  }


  # Here we define titles. They differ importantly among methods and arguments
  # in the final output list from daily_response() function
  if (is.na(result_daily_response[[5]])) {

    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Maximal correlation coefficient:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Maximal correlation coefficient:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Maximal correlation coefficient:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Correlation Coefficient")
    }







    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }

    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN with Bayesian Regularization",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }








    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }


    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN with Bayesian Regularization",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }





    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width:  day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }

    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN with Bayesian Regularization",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width:  day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }






    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "lm") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }


    # plot for lm and brnn method; using adj.r.squared
    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "brnn") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN with Bayesian Regularization",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }









    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width: day", as.numeric(row_index),
                      "days", "\nStarting day of optimal window width:",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }


    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN with Bayesian Regularization",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width: day", as.numeric(row_index),
                      "days", "\nStarting day of optimal window width:",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }











    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN with Bayesian Regularization",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")

    }


  }



  # If there is information about analysed period in the fifth element, then this
  # will be added to the plot information

  if (!is.na(result_daily_response[[5]])) {

    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMaximal correlation coefficient:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMaximal correlation coefficient:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMaximal correlation coefficient:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Correlation Coefficient")
    }







    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: Linear Regression",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }

    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: ANN with Bayesian Regularization",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }








    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: Linear Regression",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }


    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: ANN with Bayesian Regularization",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Explained Variance")
    }





    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: Linear Regression",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width:  day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }

    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: ANN with Bayesian Regularization",
                      "\nMaximal R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width:  day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }






    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "lm") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: Linear Regression",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }


    # plot for lm and brnn method; using adj.r.squared
    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "brnn") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: ANN with Bayesian Regularization",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column_extra, "of current year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }









    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: Linear Regression",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width: day", as.numeric(row_index),
                      "days", "\nStarting day of optimal window width:",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }


    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: ANN with Bayesian Regularization",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width: day", as.numeric(row_index),
                      "days", "\nStarting day of optimal window width:",
                      plot_column_extra, "of previous year")) +
        xlab("Day of Year (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }











    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: Linear Regression",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed period:", result_daily_response[[5]],
                      "\nMethod: ANN with Bayesian Regularization",
                      "\nMaximal Adjusted R squared:", calculated_measure,
                      "\nOptimal window width:", as.numeric(row_index), "days",
                      "\nStarting day of optimal window width: day",
                      plot_column)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")

    }


  }





  final_plot
}
