#' @method summary dmrs
#' @export

summary.dmrs <- function(object, ...) {

  result_daily_response <- object
  result_matrix <- as.matrix(result_daily_response[[1]])

  if (!any(is.finite(result_matrix))) {

    return("All calculations are insignificant! No summary output available.")

  } else {

    # This needs to be set to provide results in English language
    Sys.setlocale("LC_TIME", "C")

    # -------------------------------------------------------------------------
    # Helper functions
    # -------------------------------------------------------------------------

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

      # Keep the original dendroTools convention with 366 daily columns.
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

    # -------------------------------------------------------------------------
    # A) Extract basic information
    # -------------------------------------------------------------------------

    type <- as.character(result_daily_response$type)
    result_daily_element1 <- data.frame(result_daily_response[[1]],
                                        check.names = FALSE)

    reference_window <- result_daily_response$reference_window

    # -------------------------------------------------------------------------
    # B) Find the strongest absolute result
    # -------------------------------------------------------------------------

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

    row_index <- row.names(result_daily_element1)[row_position]

    if (is.null(row_index) || is.na(row_index)) {
      row_index <- as.character(row_position)
    }

    optimal_window_width <- as.numeric(row_index)
    calculated_metric <- round(calculated_value, 3)

    if (!is.null(result_daily_response$boot_lower) &&
        !is.null(result_daily_response$boot_upper) &&
        row_position <= nrow(result_daily_response$boot_lower) &&
        plot_column_source <= ncol(result_daily_response$boot_lower)) {

      lower_bound <- result_daily_response$boot_lower[row_position,
                                                      plot_column_source]
      upper_bound <- result_daily_response$boot_upper[row_position,
                                                      plot_column_source]

    } else {

      lower_bound <- NA
      upper_bound <- NA
    }

    # -------------------------------------------------------------------------
    # C) Metric and method labels
    # -------------------------------------------------------------------------

    if (result_daily_response[[2]] == "cor") {

      metric_label <- "Correlation Coefficient"

    } else if (result_daily_response[[2]] == "pcor") {

      metric_label <- "Partial Correlation Coefficient"

    } else if (result_daily_response[[3]] == "r.squared") {

      metric_label <- "Explained Variance"

    } else if (result_daily_response[[3]] == "adj.r.squared") {

      metric_label <- "Adjusted Explained Variance"

    } else {

      metric_label <- "Calculated Metric"
    }

    if (result_daily_response[[2]] == "cor") {

      method_string <- paste0("Correlation Coefficient (",
                              result_daily_response[[3]], ")")

    } else if (result_daily_response[[2]] == "pcor") {

      method_string <- paste0("Partial Correlation Coefficient (",
                              result_daily_response[[3]], ")")

    } else if (result_daily_response[[2]] == "lm") {

      method_string <- "Linear Regression"

    } else if (result_daily_response[[2]] == "brnn") {

      method_string <- "ANN with Bayesian Regularization"

    } else {

      method_string <- as.character(result_daily_response[[2]])
    }

    # -------------------------------------------------------------------------
    # D) Window start and end positions
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # E) Daily output
    # -------------------------------------------------------------------------

    if (type == "daily") {

      days_per_year <- 366L

      if (!is.null(result_daily_response$number_previous_years)) {

        number_previous_years <- result_daily_response$number_previous_years

      } else {

        number_previous_years <- ceiling(ncol(result_daily_element1) /
                                           days_per_year) - 1L
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

        reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                                   reference_info$doy)

      } else if (reference_window == "end") {

        reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                                   reference_info$doy)

      } else if (reference_window == "middle") {

        reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                                   reference_info$doy)
      }

      if (number_previous_years > 0L) {
        reference_string <- paste0(reference_string,
                                   " of ",
                                   reference_info$relative_year_label)
      }

      if (optimal_window_width == 1) {

        Optimal_string <- daily_label(optimal_start_position)

      } else {

        Optimal_string <- paste0(daily_label(optimal_start_position),
                                 " - ",
                                 daily_label(optimal_end_position))
      }

      analysed_previous_year <- number_previous_years > 0L
      optimal_reference_year <- reference_info$relative_year_label
    }

    # -------------------------------------------------------------------------
    # F) Monthly output
    # -------------------------------------------------------------------------

    if (type == "monthly") {

      months_per_year <- 12L
      month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

      if (!is.null(result_daily_response$number_previous_years)) {

        number_previous_years <- result_daily_response$number_previous_years

      } else {

        number_previous_years <- ceiling(ncol(result_daily_element1) /
                                           months_per_year) - 1L
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

        reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                   reference_info$month_number)

      } else if (reference_window == "end") {

        reference_string <- paste0("Ending Month of Optimal Window Width: Month ",
                                   reference_info$month_number)

      } else if (reference_window == "middle") {

        reference_string <- paste0("Middle Month of Optimal Window Width: Month ",
                                   reference_info$month_number)
      }

      if (number_previous_years > 0L) {
        reference_string <- paste0(reference_string,
                                   " of ",
                                   reference_info$relative_year_label)
      }

      if (optimal_window_width == 1) {

        Optimal_string <- monthly_label(optimal_start_position)

      } else {

        Optimal_string <- paste0(monthly_label(optimal_start_position),
                                 " - ",
                                 monthly_label(optimal_end_position))
      }

      analysed_previous_year <- number_previous_years > 0L
      optimal_reference_year <- reference_info$relative_year_label
    }

    # -------------------------------------------------------------------------
    # G) Final output
    # -------------------------------------------------------------------------

    output_df <- data.frame(
      Variable = c("approach",
                   "method",
                   "metric",
                   "analysed_years",
                   "maximal_calculated_metric",
                   "lower_ci",
                   "upper_ci",
                   "reference_window",
                   "analysed_previous_year",
                   "number_previous_years",
                   "optimal_reference_year",
                   "optimal_time_window",
                   "optimal_time_window_length"),

      Value = c(result_daily_response$type,
                method_string,
                metric_label,
                result_daily_response[[4]],
                calculated_metric,
                round(lower_bound, 3),
                round(upper_bound, 3),
                reference_string,
                analysed_previous_year,
                number_previous_years,
                optimal_reference_year,
                Optimal_string,
                optimal_window_width),
      stringsAsFactors = FALSE
    )

    return(output_df)
  }
}
