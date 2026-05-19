#' plot_extreme
#'
#' Graphs a line or bar plot of a row with the highest metric in a matrix,
#' produced by \code{\link{daily_response}} or \code{\link{monthly_response}}
#' functions. Bar plot is drawn for monthly_response(), while for daily_response,
#' line plot is produced.
#'
#' @param result_daily_response a list with objects as produced by
#' daily_response or monthly_response function
#' @param title logical, if set to FALSE, no plot title is displayed
#' @param ylimits limit of the y axes. It should be given as ylimits = c(0,1)
#' @param reference_window character string, the reference_window argument describes,
#' how each calculation is referred. There are three different options: 'start'
#' (default), 'end' and 'middle'.
#' @param type the character string describing type of analysis: daily or monthly
#' @param show_year_separators logical. If TRUE, dashed vertical lines are
#' added between relative-year blocks. The default is TRUE.
#'
#' @return A ggplot2 object containing the plot display
#'
#' @keywords internal

plot_extreme <- function(result_daily_response,
                         title = TRUE,
                         ylimits = NULL,
                         reference_window = "start",
                         type = "daily",
                         show_year_separators = TRUE) {

  # This needs to be set to provide results in English language
  Sys.setlocale("LC_TIME", "C")

  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------

  round2_local <- function(x, digits = 0) {
    if (exists("round2", mode = "function")) {
      round2(x, digits)
    } else {
      round(x, digits)
    }
  }

  make_daily_date_labels <- function() {

    doy <- seq(1:366)
    date <- seq(as.Date("2013-01-01"),
                as.Date("2013-12-31"),
                by = "+1 day")

    # Original dendroTools convention: 366 columns are possible.
    date[366] <- as.Date("2015-12-31")
    date <- format(date, "%b %d")

    data.frame(doy = doy, date = date, stringsAsFactors = FALSE)
  }

  get_daily_position_info <- function(position,
                                      number_previous_years,
                                      days_per_year = 366L) {

    block_index <- ceiling(position / days_per_year)
    doy <- ((position - 1L) %% days_per_year) + 1L

    if (number_previous_years == 0L) {

      relative_year_number <- 0L
      relative_year_label <- "Y"

    } else {

      relative_year_number <- block_index - (number_previous_years + 1L)

      relative_year_label <- if (relative_year_number == 0L) {
        "Y"
      } else {
        paste0("Y", relative_year_number)
      }
    }

    list(
      doy = doy,
      relative_year_number = relative_year_number,
      relative_year_label = relative_year_label
    )
  }

  get_monthly_position_info <- function(position,
                                        number_previous_years,
                                        months_per_year = 12L) {

    block_index <- ceiling(position / months_per_year)
    month_number <- ((position - 1L) %% months_per_year) + 1L

    if (number_previous_years == 0L) {

      relative_year_number <- 0L
      relative_year_label <- "Y"

    } else {

      relative_year_number <- block_index - (number_previous_years + 1L)

      relative_year_label <- if (relative_year_number == 0L) {
        "Y"
      } else {
        paste0("Y", relative_year_number)
      }
    }

    list(
      month_number = month_number,
      relative_year_number = relative_year_number,
      relative_year_label = relative_year_label
    )
  }

  # Function to find the longest non-NA range
  longest_non_na_sequence <- function(x) {

    max_length <- 0
    current_length <- 0

    for (value_temp in x) {

      if (!is.na(value_temp)) {

        current_length <- current_length + 1

      } else {

        if (current_length > max_length) {
          max_length <- current_length
        }

        current_length <- 0
      }
    }

    max(max_length, current_length)
  }

  # ---------------------------------------------------------------------------
  # Extract matrix and find optimal row/column
  # ---------------------------------------------------------------------------

  result_matrix <- as.matrix(result_daily_response[[1]])

  if (!any(is.finite(result_matrix))) {
    return("All calculations are insignificant! No plot output available.")
  }

  result_daily_element1 <- data.frame(result_daily_response[[1]],
                                      check.names = FALSE)

  n_matrix_cols <- ncol(result_matrix)

  overall_max <- max(result_matrix, na.rm = TRUE)
  overall_min <- min(result_matrix, na.rm = TRUE)

  if (abs(overall_max) >= abs(overall_min)) {

    calculated_value <- overall_max
    optimal_position <- which(result_matrix == overall_max, arr.ind = TRUE)[1, ]

  } else {

    calculated_value <- overall_min
    optimal_position <- which(result_matrix == overall_min, arr.ind = TRUE)[1, ]
  }

  row_position <- optimal_position[1]
  plot_column_source <- optimal_position[2]
  plot_column <- plot_column_source

  row_index <- row.names(result_daily_element1)[row_position]

  if (is.null(row_index) || is.na(row_index)) {
    row_index <- as.character(row_position)
  }

  optimal_window_width <- as.numeric(row_index)
  calculated_metric <- round(calculated_value, 3)

  temporal_vector <- as.numeric(result_matrix[row_position, ])

  plot_df <- data.frame(
    Position = seq_along(temporal_vector),
    Value = temporal_vector
  )

  # ---------------------------------------------------------------------------
  # General labels
  # ---------------------------------------------------------------------------

  if (result_daily_response[[2]] == "cor") {

    y_lab <- "Correlation Coefficient"

  } else if (result_daily_response[[2]] == "pcor") {

    y_lab <- "Partial Correlation Coefficient"

  } else if (result_daily_response[[3]] == "r.squared") {

    y_lab <- "Explained Variance"

  } else if (result_daily_response[[3]] == "adj.r.squared") {

    y_lab <- "Adjusted Explained Variance"

  } else {

    y_lab <- "Calculated Metric"
  }

  if (result_daily_response[[2]] == "cor") {

    method_string <- paste0("\nMethod: Correlation Coefficient (",
                            result_daily_response[[3]], ")")

  } else if (result_daily_response[[2]] == "pcor") {

    method_string <- paste0("\nMethod: Partial Correlation Coefficient (",
                            result_daily_response[[3]], ")")

  } else if (result_daily_response[[2]] == "lm") {

    method_string <- "\nMethod: Linear Regression"

  } else if (result_daily_response[[2]] == "brnn") {

    method_string <- "\nMethod: ANN with Bayesian Regularization"

  } else {

    method_string <- paste0("\nMethod: ", result_daily_response[[2]])
  }

  period_string <- paste0("Analysed Period: ", result_daily_response[[4]])
  optimal_window_string <- paste0("\nOptimal Window Width: ",
                                  optimal_window_width)

  optimal_calculation <- paste0("\nOptimal ", y_lab, ": ",
                                calculated_metric)

  # ---------------------------------------------------------------------------
  # Determine optimal start/end positions
  # ---------------------------------------------------------------------------

  if (optimal_window_width %% 2 == 0) {

    adjustment_1 <- 0
    adjustment_2 <- 1

  } else {

    adjustment_1 <- 1
    adjustment_2 <- 2
  }

  if (reference_window == "start") {

    optimal_start_position <- plot_column_source
    optimal_end_position <- plot_column_source + optimal_window_width - 1

  } else if (reference_window == "end") {

    optimal_start_position <- plot_column_source - optimal_window_width + 1
    optimal_end_position <- plot_column_source

  } else if (reference_window == "middle") {

    optimal_start_position <- round2_local(
      plot_column_source - optimal_window_width / 2, 0
    ) - adjustment_1

    optimal_end_position <- round2_local(
      plot_column_source + optimal_window_width / 2, 0
    ) - adjustment_2

  } else {

    stop("reference_window must be 'start', 'end' or 'middle'.")
  }

  # ---------------------------------------------------------------------------
  # Theme
  # ---------------------------------------------------------------------------

  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18),
          text = element_text(size = 18),
          plot.title = element_text(size = 16, face = "bold"))

  if (title == FALSE) {
    journal_theme <- journal_theme +
      theme(plot.title = element_blank())
  }

  # ---------------------------------------------------------------------------
  # Daily plot
  # ---------------------------------------------------------------------------

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

    full_sequence_width <- days_per_year * (number_previous_years + 1L)

    optimal_start_position <- max(1, optimal_start_position)
    optimal_end_position <- min(full_sequence_width, optimal_end_position)
    plot_column_source <- max(1, min(full_sequence_width,
                                     plot_column_source))

    date_codes <- make_daily_date_labels()

    daily_label <- function(position) {

      position_info <- get_daily_position_info(
        position = position,
        number_previous_years = number_previous_years,
        days_per_year = days_per_year
      )

      date_label <- date_codes$date[position_info$doy]

      if (number_previous_years > 0L) {
        paste0(position_info$relative_year_label, " ", date_label)
      } else {
        date_label
      }
    }

    reference_info <- get_daily_position_info(
      position = plot_column_source,
      number_previous_years = number_previous_years,
      days_per_year = days_per_year
    )

    if (reference_window == "start") {

      reference_string <- paste0("\nStarting Day of Optimal Window Width: Day ",
                                 reference_info$doy)

    } else if (reference_window == "end") {

      reference_string <- paste0("\nEnding Day of Optimal Window Width: Day ",
                                 reference_info$doy)

    } else if (reference_window == "middle") {

      reference_string <- paste0("\nMiddle Day of Optimal Window Width: Day ",
                                 reference_info$doy)
    }

    if (number_previous_years > 0L) {
      reference_string <- paste0(reference_string,
                                 " of ",
                                 reference_info$relative_year_label)
    }

    if (optimal_window_width == 1) {

      Optimal_string <- paste0("\nOptimal Selection: ",
                               daily_label(optimal_start_position))

    } else {

      Optimal_string <- paste0("\nOptimal Selection: ",
                               daily_label(optimal_start_position),
                               " - ",
                               daily_label(optimal_end_position))
    }

    # Axis breaks and labels, same logic as plot_heatmap()
    year_labels <- if (number_previous_years == 0L) {
      "Y"
    } else {
      c(paste0("Y-", number_previous_years:1), "Y")
    }

    year_start <- seq(1,
                      by = days_per_year,
                      length.out = length(year_labels))

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

    longest_sequence <- longest_non_na_sequence(plot_df$Value)

    if (longest_sequence > 1) {

      final_plot <- suppressWarnings(
        ggplot(plot_df, aes(x = Position, y = Value)) +
          geom_line(size = 1.2) +
          geom_vline(xintercept = plot_column_source,
                     colour = "red",
                     size = 0.8) +
          scale_x_continuous(expand = c(0, 0),
                             breaks = x_breaks,
                             labels = x_labels) +
          scale_y_continuous(limits = ylimits) +
          journal_theme
      )

    } else {

      final_plot <- suppressWarnings(
        ggplot(plot_df, aes(x = Position, y = Value)) +
          geom_point() +
          geom_vline(xintercept = plot_column_source,
                     colour = "red",
                     size = 0.8) +
          scale_x_continuous(expand = c(0, 0),
                             breaks = x_breaks,
                             labels = x_labels) +
          scale_y_continuous(limits = ylimits) +
          journal_theme
      )
    }

    if (show_year_separators && number_previous_years > 0L) {

      final_plot <- final_plot +
        geom_vline(xintercept = year_separator_positions,
                   size = 0.4,
                   linetype = "dashed",
                   colour = "grey40")
    }

    # Annotation placement
    y_min <- min(plot_df$Value, na.rm = TRUE)
    y_max <- max(plot_df$Value, na.rm = TRUE)
    y_span <- y_max - y_min

    if (!is.finite(y_span) || y_span == 0) {
      y_span <- abs(y_max)

      if (!is.finite(y_span) || y_span == 0) {
        y_span <- 1
      }
    }

    label_x <- plot_column_source + max(10, round(n_matrix_cols * 0.015))
    label_x <- min(label_x, n_matrix_cols)

    reference_label <- paste0(reference_info$relative_year_label,
                              " DOY ",
                              reference_info$doy)

    final_plot <- final_plot +
      annotate("label",
               label = as.character(calculated_metric),
               y = calculated_metric,
               x = label_x) +
      annotate("label",
               label = reference_label,
               y = y_min + 0.15 * y_span,
               x = label_x)

    x_lab <- if (number_previous_years > 0L) {
      "Relative year and day of year"
    } else {
      "Day of year"
    }

    final_plot <- final_plot +
      ggtitle(paste0(period_string,
                     method_string,
                     optimal_calculation,
                     optimal_window_string,
                     " Days",
                     reference_string,
                     Optimal_string)) +
      xlab(x_lab) +
      ylab(y_lab)
  }

  # ---------------------------------------------------------------------------
  # Monthly plot
  # ---------------------------------------------------------------------------

  if (type == "monthly") {

    months_per_year <- 12L

    month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

    if (!is.null(result_daily_response$number_previous_years)) {

      number_previous_years <- result_daily_response$number_previous_years

    } else {

      number_previous_years <- ceiling(n_matrix_cols / months_per_year) - 1L
    }

    if (is.na(number_previous_years) || number_previous_years < 0L) {
      number_previous_years <- 0L
    }

    number_previous_years <- as.integer(number_previous_years)

    full_sequence_width <- months_per_year * (number_previous_years + 1L)

    optimal_start_position <- max(1, optimal_start_position)
    optimal_end_position <- min(full_sequence_width, optimal_end_position)
    plot_column_source <- max(1, min(full_sequence_width,
                                     plot_column_source))

    monthly_label <- function(position) {

      position_info <- get_monthly_position_info(
        position = position,
        number_previous_years = number_previous_years,
        months_per_year = months_per_year
      )

      month_label <- month_names[position_info$month_number]

      if (number_previous_years > 0L) {
        paste0(position_info$relative_year_label, " ", month_label)
      } else {
        month_label
      }
    }

    reference_info <- get_monthly_position_info(
      position = plot_column_source,
      number_previous_years = number_previous_years,
      months_per_year = months_per_year
    )

    if (reference_window == "start") {

      reference_string <- paste0("\nStarting Month of Optimal Window Width: Month ",
                                 reference_info$month_number)

    } else if (reference_window == "end") {

      reference_string <- paste0("\nEnding Month of Optimal Window Width: Month ",
                                 reference_info$month_number)

    } else if (reference_window == "middle") {

      reference_string <- paste0("\nMiddle Month of Optimal Window Width: Month ",
                                 reference_info$month_number)
    }

    if (number_previous_years > 0L) {
      reference_string <- paste0(reference_string,
                                 " of ",
                                 reference_info$relative_year_label)
    }

    if (optimal_window_width == 1) {

      Optimal_string <- paste0("\nOptimal Selection: ",
                               monthly_label(optimal_start_position))

    } else {

      Optimal_string <- paste0("\nOptimal Selection: ",
                               monthly_label(optimal_start_position),
                               " - ",
                               monthly_label(optimal_end_position))
    }

    if (optimal_window_width == 1) {
      month_string <- " Month"
    } else {
      month_string <- " Months"
    }

    optimal_window_string <- paste0("\nOptimal Window Width: ",
                                    optimal_window_width,
                                    month_string)

    year_labels <- if (number_previous_years == 0L) {
      "Y"
    } else {
      c(paste0("Y-", number_previous_years:1), "Y")
    }

    year_start <- seq(1,
                      by = months_per_year,
                      length.out = length(year_labels))

    x_breaks <- c()
    x_labels <- c()

    for (i in seq_along(year_start)) {

      start_i <- year_start[i]

      candidate_breaks <- start_i:(start_i + months_per_year - 1)
      candidate_breaks <- candidate_breaks[
        candidate_breaks <= n_matrix_cols
      ]

      candidate_labels <- month_names[seq_along(candidate_breaks)]
      candidate_labels[1] <- year_labels[i]

      x_breaks <- c(x_breaks, candidate_breaks)
      x_labels <- c(x_labels, candidate_labels)
    }

    year_separator_positions <- year_start[-1] - 0.5
    year_separator_positions <- year_separator_positions[
      year_separator_positions <= n_matrix_cols
    ]

    plot_df$BarColour <- ifelse(plot_df$Position == plot_column_source,
                                "optimal", "other")

    final_plot <- suppressWarnings(
      ggplot(plot_df, aes(x = Position, y = Value, fill = BarColour)) +
        geom_col() +
        geom_hline(yintercept = 0) +
        scale_x_continuous(expand = c(0, 0),
                           breaks = x_breaks,
                           labels = x_labels) +
        scale_y_continuous(limits = ylimits) +
        scale_fill_manual(values = c("other" = "grey50",
                                     "optimal" = "red"),
                          guide = "none") +
        journal_theme
    )

    if (show_year_separators && number_previous_years > 0L) {

      final_plot <- final_plot +
        geom_vline(xintercept = year_separator_positions,
                   size = 0.4,
                   linetype = "dashed",
                   colour = "grey40")
    }

    x_lab <- if (number_previous_years > 0L) {
      "Relative year and month"
    } else {
      "Month"
    }

    final_plot <- final_plot +
      annotate("label",
               label = as.character(calculated_metric),
               y = calculated_metric,
               x = plot_column_source) +
      xlab(x_lab) +
      ylab(y_lab) +
      ggtitle(paste0(period_string,
                     method_string,
                     optimal_calculation,
                     optimal_window_string,
                     reference_string,
                     Optimal_string))
  }

  final_plot
}
