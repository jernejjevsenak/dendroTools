#' plot_heatmap
#'
#' Graphs a heatmap of values stored in a matrix, such as produced
#' by \code{\link{daily_response}} function.
#'
#' @param result_daily_response a list with three objects as produced by
#' \code{\link{daily_response}} function
#' @param reference_window character string describing how to relate doy and
#' window of each calculation. There are two options: 'start' (default) and
#' 'end'. If the reference window is set as 'start', then each calculation is
#' related to the starting day of window, otherwise it is related to the ending
#' day of window calculation. For example, if we consider correlations with window
#' from doy 15 to doy 35. If reference window is set to ‘start’, then this calculation
#' will be related to the doy 15. If the reference window is set to ‘end’, then this
#' calculation will be related to the doy 35.
#'
#' @return A ggplot2 object containing the heatmap display
#'
#' @examples
#' \dontrun{
#' data(LJ_daily_temperatures)
#' data(example_proxies_1)
#' Example1 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "lm", metric = "r.squared",
#' fixed_width = 90, previous_year = TRUE, row_names_subset = TRUE)
#' plot_heatmap(Example1)
#'
#' Example2 <- daily_response(response = example_proxies_1,
#' env_data = LJ_daily_temperatures, method = "brnn",
#' metric = "adj.r.squared", lower_limit = 50, upper_limit = 55,
#' row_names_subset = TRUE)
#' plot_heatmap(Example2)
#'
#' library(dplyr)
#' oxygen_isotope <- dplyr::select(example_proxies_1, O18)
#' Example3 <- daily_response(response = oxygen_isotope,
#' env_data = LJ_daily_temperatures, method = "cor", lower_limit = 50,
#' upper_limit = 55, previous_year = TRUE, row_names_subset = TRUE)
#' plot_heatmap(Example3)
#' }
#'
#' @keywords internal

plot_heatmap <- function(result_daily_response, reference_window = "start"){

  # Extracting a matrix from a list and converting it into a data frame
  result_daily_element1 <- data.frame(result_daily_response [[1]])


  # Creating a nice string that will be used to generate ggplot Legend
  if (result_daily_response[[2]] == "cor"){
    temp_string <- "Correlation Coefficient"
  } else if (result_daily_response[[3]] == "r.squared"){
    temp_string <- "Explained Variance"
  } else if (result_daily_response[[3]] == "adj.r.squared"){
    temp_string <- "Adjusted Explained Variance"
  }

  # Data manipulation. The goal of this part is to prepare data for ggplot
  result_daily_element1$temp_row_names <- row.names(result_daily_element1)
  result_daily_element1_melted <- melt(result_daily_element1,
                                       id.vars = "temp_row_names")

  # colname is changed, for a more sufficient plotting
  colnames(result_daily_element1_melted)[3] <- "Value"

  # Calculating parameters for heatmap. Our goal is to expose /
  # point out extreme values.
  min_limit <- min(result_daily_element1_melted$Value, na.rm = TRUE)
  max_limit <- max(result_daily_element1_melted$Value, na.rm = TRUE)
  bounds <- quantile(result_daily_element1_melted$Value,
    probs = seq(0, 1, 0.01), na.rm = TRUE)
  bound1 <- bounds[1]
  bound2 <- bounds[20]
  bound3 <- bounds[50]
  bound4 <- bounds[100]

  # When the matrix in result_daily_element_1 is small, for the conviniece,
  # different solution is needed.
    if (nrow(result_daily_element1) * ncol(result_daily_element1) < 500){
    bounds <- quantile(result_daily_element1_melted$Value,
      probs = seq(0, 1, 0.1), na.rm = TRUE)
    bound1 <- bounds[1]
    bound2 <- bounds[2]
    bound3 <- bounds[5]
    bound4 <- bounds[11]
  }

  # The definition of theme
  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18), text = element_text(size = 18),
          plot.title = element_text(size = 16,  face = "bold"),
          legend.position = "bottom", legend.key.width = unit(3, "line"))


  final_plot <- suppressWarnings(ggplot(result_daily_element1_melted,
    aes_(x = ~as.numeric(variable), y = ~as.numeric(temp_row_names),
    fill = ~Value)) +
    geom_tile() +
    scale_fill_gradientn(temp_string,
                         colours = c("blue", "red", "yellow", "black"),
      values = rescale(c(bound1, bound2, bound3, bound4)),
     guide = "colorbar", limits = c(min_limit, max_limit),
     na.value = "white") +
    xlab("Day of Year") +
    ylab("Window Width") +
    scale_x_continuous(expand = c(0, 0), breaks = sort(c(seq(0, nlevels(result_daily_element1_melted$variable), 50)),
                                                       decreasing = FALSE),
                         labels = sort(c(seq(0, nlevels(result_daily_element1_melted$variable), 50)))) +
    journal_theme)

  # Scale_y_continuous is added separately. When there is only a few  rows
  # e.g. fixed_width = TRUE, breaks are specified separately

  if (nrow(result_daily_element1) < 5) {
    final_plot <- final_plot +
      scale_y_continuous(expand = c(0, 0),
                         breaks = pretty_breaks(n =
                                  nrow(result_daily_element1)))
  } else {
    final_plot <- final_plot +
      scale_y_continuous(expand = c(0, 0),
                         breaks = pretty_breaks())
  }

  # If previous_year == TRUE(function daily_response), different xlab
  # is needed
  if (ncol(result_daily_element1) > 367) {
    final_plot <- suppressWarnings(final_plot +
      annotate(fontface = "bold", label = 'Previous Year', geom = 'label',
              x = 366 - ncol(result_daily_element1) / 12.8,
              y = round(max(as.numeric(row.names(result_daily_element1))) -
                         (max(as.numeric(row.names(result_daily_element1))) -
                            min(as.numeric(row.names(result_daily_element1)))) / 5), 0) +
      annotate(fontface = "bold", label = 'Current Year', geom = 'label',
               x = 366 + ncol(result_daily_element1) / 13.5,
               y = round(max(as.numeric(row.names(result_daily_element1))) -
                           (max(as.numeric(row.names(result_daily_element1))) -
                              min(as.numeric(row.names(result_daily_element1)))) / 5), 0) +
      geom_vline(xintercept = 366, size = 1) +
      xlab("Day of Year  (Including Previous Year)"))
  }


  if (ncol(result_daily_element1) > 367){
    x_lab <- "Day of Year  (Including Previous Year)"
  } else if (ncol(result_daily_element1) <= 367){
    x_lab <- "Day of Year"
  }

  if (result_daily_response[[2]] == 'cor'){
    method_string <- paste0("\nMethod: Pearson Correlation")
  } else if (result_daily_response[[2]] == 'lm'){
    method_string <- paste0("\nMethod: Linear Regression")
  } else if (result_daily_response[[2]] == 'brnn'){
    method_string <- paste0("\nMethod: ANN with Bayesian Regularization")
  }

  period_string <- paste0("\nAnalysed Period: ", result_daily_response[[4]])

  if (reference_window == 'start'){
    reference_string <- "\nDOY Reference of Each Calculation is the Beginning of the Window"
  } else if (reference_window == 'end'){
    reference_string <- "\nDOY Reference of Each Calculation is the End of the Window"
  }

  final_plot <- final_plot +
        xlab(x_lab) +
        ggtitle(paste0(period_string, method_string, reference_string))

  final_plot
}

