#' plot_heatmap
#'
#' Graphs a heatmap of values stored in a matrix, such as produced
#' by \code{\link{daily_response}} or \code{\link{monthly_response}} functions.
#'
#' @param result_daily_response a list with objects as produced by
#' \code{\link{daily_response}}, \code{\link{monthly_response}},
#' \code{\link{daily_response_seascorr}} or
#' \code{\link{monthly_response_seascorr}}.
#' @param reference_window character string, the reference_window argument describes
#' how each calculation is referred. Options are 'start', 'end' and 'middle'.
#' If available, the value stored in result_daily_response$reference_window is used.
#' @param type the character string describing type of analysis: daily or monthly.
#' @param show_title logical. If TRUE, the plot title with analysed period, method
#' and reference-window information is shown. The default is FALSE.
#' @param show_year_separators logical. If TRUE, dashed vertical lines are
#' added between relative-year blocks in plots with previous-year data.
#' The default is TRUE.
#'
#' @return A ggplot2 object containing the heatmap display
#'
#' @keywords internal

plot_heatmap <- function(result_daily_response,
                         reference_window = "start",
                         type = "daily",
                         show_title = FALSE,
                         show_year_separators = TRUE) {

  variable <- NULL
  temp_row_names <- NULL
  Value <- NULL

  # Prefer metadata stored in the result object
  if (!is.null(result_daily_response$reference_window)) {
    reference_window <- result_daily_response$reference_window
  }

  reference_position_label <- if (reference_window == "start") {
    "start of window"
  } else if (reference_window == "end") {
    "end of window"
  } else if (reference_window == "middle") {
    "middle of window"
  } else {
    reference_window
  }

  # Extract matrix and keep original column names
  result_daily_element1 <- data.frame(result_daily_response[[1]],
                                      check.names = FALSE)

  n_matrix_cols <- ncol(result_daily_element1)

  a <- max(result_daily_element1, na.rm = TRUE)
  b <- min(result_daily_element1, na.rm = TRUE)

  if ((a > 0 & b > 0) | (a < 0 & b < 0)) {
    paleta <- "viridis"
  } else {
    paleta <- "no_viridis"
  }

  # Strings for optional title
  period_string <- paste0("\nAnalysed Period: ", result_daily_response[[4]])

  if (result_daily_response[[2]] == "cor") {
    method_string <- paste0("\nMethod: Correlation Coefficient (",
                            result_daily_response[[3]], ")")
  } else if (result_daily_response[[2]] == "pcor") {
    method_string <- paste0("\nMethod: Partial Correlation Coefficient (",
                            result_daily_response[[3]], ")")
  } else if (result_daily_response[[2]] == "lm") {
    method_string <- paste0("\nMethod: Linear Regression")
  } else if (result_daily_response[[2]] == "brnn") {
    method_string <- paste0("\nMethod: ANN with Bayesian Regularization")
  } else {
    method_string <- paste0("\nMethod: ", result_daily_response[[2]])
  }

  if (reference_window == "start") {
    reference_string <- "\nDOY Reference of Each Calculation is the Beginning of the Window"
  } else if (reference_window == "end") {
    reference_string <- "\nDOY Reference of Each Calculation is the End Day of the Window"
  } else if (reference_window == "middle") {
    reference_string <- "\nDOY Reference of Each Calculation is the Middle Day of the Window"
  } else {
    reference_string <- ""
  }

  # Legend label
  if (result_daily_response[[2]] == "cor") {
    temp_string <- "Correlation Coefficient"
  } else if (result_daily_response[[2]] == "pcor") {
    temp_string <- "Partial Correlation Coefficient"
  } else if (result_daily_response[[3]] == "r.squared") {
    temp_string <- "Explained Variance"
  } else if (result_daily_response[[3]] == "adj.r.squared") {
    temp_string <- "Adjusted Explained Variance"
  } else {
    temp_string <- "Value"
  }

  # Prepare data for ggplot
  result_daily_element1$temp_row_names <- row.names(result_daily_element1)

  result_daily_element1_melted <- melt(result_daily_element1,
                                       id.vars = "temp_row_names")

  colnames(result_daily_element1_melted)[3] <- "Value"

  # Ensure x positions are numeric
  result_daily_element1_melted$variable <- as.numeric(gsub(
    "^X", "", as.character(result_daily_element1_melted$variable)
  ))

  result_daily_element1_melted <- result_daily_element1_melted[
    !is.na(result_daily_element1_melted$variable), , drop = FALSE
  ]

  # Colour scaling
  min_limit <- min(result_daily_element1_melted$Value, na.rm = TRUE)
  max_limit <- max(result_daily_element1_melted$Value, na.rm = TRUE)

  bounds <- quantile(result_daily_element1_melted$Value,
                     probs = seq(0, 1, 0.01),
                     na.rm = TRUE)

  bound1 <- bounds[1]
  bound2 <- bounds[20]
  bound3 <- bounds[50]
  bound4 <- bounds[100]

  if (nrow(result_daily_response[[1]]) * n_matrix_cols < 500) {

    bounds <- quantile(result_daily_element1_melted$Value,
                       probs = seq(0, 1, 0.1),
                       na.rm = TRUE)

    bound1 <- bounds[1]
    bound2 <- bounds[2]
    bound3 <- bounds[5]
    bound4 <- bounds[11]
  }

  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18),
          text = element_text(size = 18),
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom",
          legend.key.width = unit(3, "line"))

  if (type == "daily") {

    days_per_year <- 366L

    if (!is.null(result_daily_response$number_previous_years)) {
      number_previous_years <- result_daily_response$number_previous_years
    } else {
      number_previous_years <- ceiling(n_matrix_cols / days_per_year) - 1L
    }

    if (is.na(number_previous_years) || number_previous_years < 0L) {
      number_previous_years <- 0L
    }

    number_previous_years <- as.integer(number_previous_years)

    year_labels <- if (number_previous_years == 0L) {
      "Y"
    } else {
      c(paste0("Y-", number_previous_years:1), "Y")
    }

    year_start <- seq(1,
                      by = days_per_year,
                      length.out = length(year_labels))

    # Breaks: first day of each year block gets only Y-label,
    # then 100, 200 and 300 are repeated within each block.
    x_breaks <- c()
    x_labels <- c()

    for (i in seq_along(year_start)) {

      start_i <- year_start[i]

      candidate_breaks <- c(start_i,
                            start_i + 99,
                            start_i + 199,
                            start_i + 299)

      candidate_breaks <- candidate_breaks[
        candidate_breaks <= n_matrix_cols
      ]

      candidate_labels <- c(year_labels[i], "100", "200", "300")
      candidate_labels <- candidate_labels[seq_along(candidate_breaks)]

      x_breaks <- c(x_breaks, candidate_breaks)
      x_labels <- c(x_labels, candidate_labels)
    }

    year_separator_positions <- year_start[-1] - 0.5
    year_separator_positions <- year_separator_positions[
      year_separator_positions <= n_matrix_cols
    ]

    x_axis_label <- if (number_previous_years > 0L) {
      paste0("Relative year and DOY (", reference_position_label, ")")
    } else {
      paste0("DOY (", reference_position_label, ")")
    }

    final_plot <- suppressWarnings(
      ggplot(result_daily_element1_melted,
             aes(x = variable,
                 y = as.numeric(temp_row_names),
                 fill = Value)) +
        geom_tile() +
        xlab(x_axis_label) +
        ylab("Window Width") +
        scale_x_continuous(
          expand = c(0, 0),
          breaks = x_breaks,
          labels = x_labels
        ) +
        journal_theme +
        theme(axis.text.x.bottom = element_text(size = 13))
    )

    if (show_year_separators && number_previous_years > 0L) {
      final_plot <- final_plot +
        geom_vline(xintercept = year_separator_positions,
                   linewidth = 0.4,
                   linetype = "dashed",
                   colour = "grey40")
    }

    if (paleta == "viridis") {

      final_plot <- final_plot +
        scale_fill_viridis(temp_string,
                           na.value = "white",
                           direction = -1)

    } else {

      final_plot <- final_plot +
        scale_fill_gradientn(temp_string,
                             colours = c("red4", "orange", "cyan", "blue4"),
                             values = rescale(c(bound1, bound2, bound3, bound4)),
                             guide = "colorbar",
                             limits = c(min_limit, max_limit),
                             na.value = "white")
    }

    if (nrow(result_daily_response[[1]]) < 5) {

      final_plot <- final_plot +
        scale_y_continuous(expand = c(0, 0),
                           breaks = pretty_breaks(
                             n = nrow(result_daily_response[[1]])
                           ))

    } else {

      final_plot <- final_plot +
        scale_y_continuous(expand = c(0, 0),
                           breaks = pretty_breaks())
    }

    if (show_title) {
      final_plot <- final_plot +
        ggtitle(paste0(period_string, method_string, reference_string))
    }
  }

  if (type == "monthly") {

    months_per_year <- 12L

    if (!is.null(result_daily_response$number_previous_years)) {
      number_previous_years <- result_daily_response$number_previous_years
    } else {
      number_previous_years <- ceiling(n_matrix_cols / months_per_year) - 1L
    }

    if (is.na(number_previous_years) || number_previous_years < 0L) {
      number_previous_years <- 0L
    }

    number_previous_years <- as.integer(number_previous_years)

    year_labels <- if (number_previous_years == 0L) {
      "Y"
    } else {
      c(paste0("Y-", number_previous_years:1), "Y")
    }

    year_start <- seq(1,
                      by = months_per_year,
                      length.out = length(year_labels))

    month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

    x_breaks <- c()
    x_labels <- c()

    for (i in seq_along(year_start)) {

      start_i <- year_start[i]

      candidate_breaks <- start_i:(start_i + months_per_year - 1)
      candidate_breaks <- candidate_breaks[
        candidate_breaks <= n_matrix_cols
      ]

      candidate_labels <- month_labels[seq_along(candidate_breaks)]

      if (number_previous_years > 0L && length(candidate_labels) > 0) {
        candidate_labels[1] <- year_labels[i]
      }

      x_breaks <- c(x_breaks, candidate_breaks)
      x_labels <- c(x_labels, candidate_labels)
    }

    year_separator_positions <- year_start[-1] - 0.5
    year_separator_positions <- year_separator_positions[
      year_separator_positions <= n_matrix_cols
    ]

    x_axis_label <- if (number_previous_years > 0L) {
      paste0("Relative year and month (", reference_position_label, ")")
    } else {
      paste0("Month (", reference_position_label, ")")
    }

    final_plot <- suppressWarnings(
      ggplot(result_daily_element1_melted,
             aes(x = variable,
                 y = as.numeric(temp_row_names),
                 fill = Value)) +
        geom_tile() +
        ylab("Number of Consecutive Months") +
        xlab(x_axis_label) +
        scale_x_continuous(expand = c(0, 0),
                           breaks = x_breaks,
                           labels = x_labels) +
        scale_y_continuous(expand = c(0, 0),
                           breaks = pretty_breaks()) +
        journal_theme +
        theme(axis.text.x.bottom = element_text(size = 13))
    )

    if (show_year_separators && number_previous_years > 0L) {
      final_plot <- final_plot +
        geom_vline(xintercept = year_separator_positions,
                   linewidth = 0.4,
                   linetype = "dashed",
                   colour = "grey40")
    }

    if (show_title) {
      final_plot <- final_plot +
        ggtitle(paste0(period_string, method_string))
    }

    if (paleta == "viridis") {

      final_plot <- final_plot +
        scale_fill_viridis(temp_string,
                           na.value = "white",
                           direction = -1)

    } else {

      final_plot <- final_plot +
        scale_fill_gradientn(temp_string,
                             colours = c("red4", "orange", "cyan", "blue4"),
                             values = rescale(c(bound1, bound2, bound3, bound4)),
                             guide = "colorbar",
                             limits = c(min_limit, max_limit),
                             na.value = "white")
    }
  }

  final_plot
}
