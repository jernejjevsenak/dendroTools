#' plot_specific
#'
#' Graphs a line plot of a row with a selected window width in a matrix,
#' produced by \code{\link{daily_response}} function.
#'
#' @param result_daily_response a list with three objects as produced by
#' daily_response function
#' @param window_width integer representing window width to be displayed
#' @param title logical, if set to FALSE, no plot title is displayed
#' @param ylimits limit of the y axes. It should be given as ylimits = c(0,1)
#'
#' @return A ggplot2 object containing the plot display
#'
#' @examples
#' }
#' data(LJ_daily_temperatures)
#' data(KRE_daily_temperatures)
#' data(example_proxies_1)
#' Example1 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "lm", metric = "r.squared",
#' lower_limit = 90, upper_limit = 150, row_names_subset = TRUE,
#' previous_year = TRUE)
#' plot_specific(Example1, window_width = 90)
#'
#' Example2 <- daily_response(response = data_TRW_1,
#' env_data = KRE_daily_temperatures, method = "cor",
#' metric = "adj.r.squared", lower_limit = 150, upper_limit = 155,
#' neurons = 1, row_names_subset = TRUE, previous_year = TRUE)
#' plot_specific(Example2, window_width = 153, title = TRUE)
#'
#' Example3 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "brnn",
#' metric = "adj.r.squared", lower_limit = 150, upper_limit = 155,
#' neurons = 1, previous_year = TRUE, row_names_subset = TRUE)
#' plot_specific(Example3, window_width = 153, title = TRUE)
#' }
#'
#' @keywords internal

plot_specific <- function(result_daily_response, window_width, title = TRUE,
                          ylimits = NULL) {

  # Short description of the function. It
  # - extracts matrix (the frst object of a list)
  # - verification of whether we deal with negative correlations. In this case
  # we will expose globail minimum (and not maximum, as in the case of positive
  # correlations, r.squared and adj.r.squared)
  # - subseting extracted matrix to keep only row, as defined by the argument
  # window_width
  # In case of there are more than 366 columns (days), xlabs are properly
  # labeled

  # A) Extracting a matrix from a list and converting it into a data frame
  result_daily_element1 <- data.frame(result_daily_response[[1]])

  # warning msg in case of selected window_width not among row.names.
  # support_string suggests, which window_widths are avaliable.
  support_string <- paste(row.names(result_daily_element1), sep = "",
    collapse = ", ")

  if (as.character(window_width) %in% row.names(result_daily_element1)
    == FALSE) {
    stop(paste("Selected window_width is not avaliable.",
                "Select one among:", support_string, sep = ""))
  }

  # Subseting result_daily_element1 to keep only row with selected window_width.
  # Subset is transposed and converted to a data frame, so data is ready for
  # ggplot2
  temporal_vector <-
    data.frame(t(result_daily_element1[row.names(result_daily_element1)
    == as.character(window_width), ]))

  # Removing missing values at the end of tempora_vector
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
    temporal_vector <- temporal_vector[-c(row_count:(row_count +
                                                         delete_rows)), ]
  }
  temporal_vector <- data.frame(temporal_vector)




  # renaming the first column.
  # (I later use this name in the script for plotting)
  names(temporal_vector) <- c("var1")


  # B) Sometime we have the case of negative correlations, the following code
  # examines minimum and compare it with the maximum, to get the information if
  # we have negative values. In this case, we will not be looking for max, but
  # for min!

  # With the following chunk, overall_maximum and overall_minimum values of
  # subset are calculated and compared.
  overall_max <- max(temporal_vector$var1, na.rm = TRUE)
  overall_min <- min(temporal_vector$var1, na.rm = TRUE)

  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  if ((abs(overall_max) > abs(overall_min)) == TRUE) {

    # maximum value is calculated and index of column is stored as index
    # index represent the starting day (location) in the matrix, which gives
    # the maximum result
    max_result <- max(temporal_vector, na.rm = TRUE)
    calculated_metric <- round(max_result, 3)
    index <- which(temporal_vector$var1 == max_result, arr.ind = TRUE)
    plot_column <- index
  }

  if ((abs(overall_max) < abs(overall_min)) == TRUE) {

    # This is in case of negative values
    # minimum value is calculated and index of column is stored as index
    # index represent the starting day (location) in the matrix, which gives
    # the minimum result
    min_result <- min(temporal_vector, na.rm = TRUE)
    calculated_metric <- round(min_result, 3)
    index <- which(temporal_vector$var1 == min_result, arr.ind = TRUE)
    plot_column <- index
  }

  # In case of we have more than 366 days, we calculate the day of a year
  # (plot_column), considering 366 days of previous year.
  if (nrow(temporal_vector) > 366 & plot_column > 366) {
    plot_column_extra <- plot_column %% 366
  } else {
    plot_column_extra <- plot_column
  }



  # C) The final plot is being created. The first part of a plot is universal,
  # the second part defines xlabs, ylabs and ggtitles.

  # The definition of theme
  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18), text = element_text(size = 18),
          plot.title = element_text(size = 16,  face = "bold"))

  if (title == FALSE){
    journal_theme <- journal_theme +
      theme(plot.title = element_blank())
  }


  # Here we define a data frame of dates and corresponing day of year (doi). Later
  # this dataframe will be used to describe tht optimal sequence od days
  doy <- seq(1:366)
  date <- seq(as.Date('2013-01-01'),as.Date('2013-12-31'), by = "+1 day")
  date[366] <- as.Date('2015-12-31')
  date <- format(date, "%b %d")
  date_codes <- data.frame(doy = doy, date = date)


  Optimal_string <- paste("\nOptimal Selection:",
                          as.character(date_codes[plot_column_extra, 2]),"-",
                          as.character(date_codes[plot_column_extra + as.numeric(window_width) - 1, 2]))

  final_plot <- suppressWarnings(
  ggplot(temporal_vector, aes_(y = ~var1,
    x = ~ seq(1, nrow(temporal_vector)))) + geom_line(lwd = 1.2) +
    geom_vline(xintercept = plot_column, col = "red", size = 1) +
    scale_x_continuous(breaks = sort(c(seq(0, nrow(temporal_vector), 50)), decreasing = FALSE),
      labels = sort(c(seq(0, nrow(temporal_vector), 50)))) +
  scale_y_continuous(limits = ylimits) +
    annotate("label", label = as.character(calculated_metric),
      y = calculated_metric, x = plot_column + 15) +
    annotate("label", label = paste("Day", as.character(plot_column), sep = " "),
             y = min(temporal_vector, na.rm = TRUE) + 0.2*min(temporal_vector, na.rm = TRUE), x = plot_column + 15) +
    journal_theme)

  # If previous_year = TRUE, we add a vertical line with labels of
  # previous and current years
  if (ncol(result_daily_element1) > 366) {
    final_plot <- final_plot +
      annotate(fontface = "bold", label = 'Previous Year', geom = 'label',
               x = 366 - ncol(result_daily_element1) / 12.8,
               y = calculated_metric - (calculated_metric/5)) +
      annotate(fontface = "bold", label = 'Current Year', geom = 'label',
               x = 366 + ncol(result_daily_element1) / 13.5,
               y = calculated_metric -(calculated_metric/5)) +
      geom_vline(xintercept = 366, size = 1)
  }


  # Here we define titles. They differ importantly among methods and arguments
  # in the final output list from daily_response() function
  if (is.na(result_daily_response[[4]])) {


    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Maximal correlation coefficient:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Maximal correlation coefficient:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Previous Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Maximal correlation coefficient:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Correlation Coefficient")
    }







    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }


    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN With Bayesian Regularization",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }







    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Previous Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }

    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN With Bayesian Regularization",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Previous Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }







    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }

    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN With Bayesian Regularization",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }











    # plot for lm and brnn method; using adj.r.squared
    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "lm") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }



    # plot for lm and brnn method; using adj.r.squared
    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "brnn") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN With Bayesian Regularization",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }









    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra,
                      "of Previous Year", Optimal_string, "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }

    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN With Bayesian Regularization",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra,
                      "of Previous Year", Optimal_string, "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }











    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: Linear Regression",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Method: ANN With Bayesian Regularization",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")
    }


  }


  if (!is.na(result_daily_response[[4]])) {


    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMaximal correlation coefficient:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMaximal correlation coefficient:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Previous Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Correlation Coefficient")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "cor")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMaximal correlation coefficient:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Correlation Coefficient")
    }







    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: Linear Regression",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }


    # plot for lm and brnn method; using r.squared
    if ((nrow(temporal_vector) > 366) &&  (plot_column > 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: ANN With Bayesian Regularization",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }







    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: Linear Regression",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Previous Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }

    if ((nrow(temporal_vector) > 366) && (plot_column < 366) &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: ANN With Bayesian Regularization",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Previous Year", Optimal_string,
                      "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Explained Variance")
    }







    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "lm") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: Linear Regression",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }

    if (nrow(temporal_vector) < 366 &&
        (result_daily_response[[2]] == "brnn") &&
        result_daily_response[[3]] == "r.squared") {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: ANN With Bayesian Regularization",
                      "\nMaximal R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Explained Variance")
    }











    # plot for lm and brnn method; using adj.r.squared
    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "lm") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: Linear Regression",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }



    # plot for lm and brnn method; using adj.r.squared
    if ((nrow(temporal_vector) > 366) && (plot_column > 366) &&
        (result_daily_response[[2]] == "brnn") &&
        (result_daily_response[[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: ANN With Bayesian Regularization",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, "of Current Year", Optimal_string)) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }









    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: Linear Regression",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra,
                      "of Previous Year", Optimal_string, "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }

    if ((nrow(temporal_vector) > 366) &&  (plot_column < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: ANN With Bayesian Regularization",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra,
                      "of Previous Year", Optimal_string, "(Previous Year)")) +
        xlab("Day of Year  (Including Previous Year)") +
        ylab("Adjusted Explained Variance")
    }











    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "lm") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: Linear Regression",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")
    }

    if ((nrow(temporal_vector) < 366) &&
        (result_daily_response [[2]] == "brnn") &&
        (result_daily_response [[3]] == "adj.r.squared")) {
      final_plot <- final_plot +
        ggtitle(paste("Analysed Period:", result_daily_response[[4]],
                      "\nMethod: ANN With Bayesian Regularization",
                      "\nMaximal Adjusted R Squared:", calculated_metric,
                      "\nSelected Window Width:", window_width, "Days",
                      "\nStarting Day of Selected Window Width: Day",
                      plot_column_extra, Optimal_string)) +
        xlab("Day of Year") +
        ylab("Adjusted Explained Variance")
    }

  }

  final_plot
}
