#' plot_heatmap
#'
#' Graphs a heatmap of values stored in a matrix, such as produced
#' by \code{\link{daily_response}} function.
#'
#' @param result_daily_response a list with three objects as produced by
#' \code{\link{daily_response}} function
#'
#' @return A ggplot2 object containing the heatmap display
#' @export
#'
#' @examples
#' \dontrun{
#' data(daily_temperatures_example)
#' data(example_proxies_1)
#' Example1 <- daily_response(response = example_proxies_1,
#' env_data = daily_temperatures_example, method = "lm", measure = "r.squared",
#' fixed_width = 90, previous_year = TRUE)
#' plot_heatmap(Example1)
#'
#' Example2 <- daily_response(response = example_proxies_1,
#' env_data = daily_temperatures_example, method = "lm",
#' measure = "adj.r.squared", lower_limit = 50, upper_limit = 55)
#' plot_heatmap(Example2)
#' }

plot_heatmap <- function(result_daily_response){

  # Extracting a matrix from a list and converting it into a data frame
  result_daily_element1 <- data.frame(result_daily_response [[1]])

  # Creating a nice string that will be used to generate ggplot Legend
  if  (result_daily_response[[3]] == "r.squared"){
    temp_string <- "R squared"
  } else if (result_daily_response[[3]] == "adj.r.squared"){
    temp_string <- "Adjusted R squared"
  } else if (result_daily_response[[2]] == "cor") {
    temp_string <- "Correlation coefficient"
  } else {
    stop("Check your method and measures")
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
          legend.position = "bottom", legend.key.width = unit(3, "line"),
          plot.title = element_blank())

  final_plot <- ggplot(result_daily_element1_melted,
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
    scale_x_continuous(expand = c(0, 0)) +
    journal_theme

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
  if (ncol(result_daily_element1) > 366) {
    final_plot <- final_plot +
      xlab("Day of Year  (Including Previous Year)")
  }

  final_plot
}
