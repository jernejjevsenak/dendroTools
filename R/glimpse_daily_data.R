#' glimpse_daily_data
#'
#' Visual presentation of daily data to spot missing values.
#'
#' @param env_data a data frame of daily sequences of environmental data as
#' columns and years as row names. Each row represents a year and
#' each column represents a day of a year. Alternatively, env_data could be a
#' tidy data with three columns, i.e. Year, DOY and third column representing
#' values of mean temperatures, sum of precipitation etc. If tidy data is passed
#' to the function, set the argument tidy_env_data to TRUE.
#' @param tidy_env_data if set to TRUE, env_data should be inserted as a data frame with three
#' columns: "Year", "DOY", "Precipitation/Temperature/etc."
#' @param na.color color to use for missing values
#' @param low_color colours for low end of the gradient
#' @param high_color colours for high end of the gradient
#'
#' @examples
#' data("LJ_daily_temperatures")
#' glimpse_daily_data(env_data = LJ_daily_temperatures, tidy_env_data = FALSE, na.color = "white")
#'
#' data("LJ_daily_precipitation")
#' glimpse_daily_data(env_data = LJ_daily_precipitation, tidy_env_data = TRUE, na.color = "white")

glimpse_daily_data <- function(env_data, na.color = "red", low_color = "blue", high_color = "green",
                               tidy_env_data = FALSE){

  if (tidy_env_data == FALSE){

  a <- seq(from = 1, to = ncol(env_data), by = 1)
  colnames(env_data) <- a

  # Data manipulation. The goal of this part is to prepare data for ggplot
  env_data$temp_row_names <- row.names(env_data)
  env_data_melted <- melt(env_data, id.vars = "temp_row_names")
  colnames(env_data_melted)[3] <- "Legend"
  } else if (tidy_env_data == TRUE){
    colnames(env_data) <- c("temp_row_names", "variable", "Legend")
    env_data_melted <- env_data
  }

  min_year = min(as.numeric(row.names(env_data)))
  max_year = max(as.numeric(row.names(env_data)))

  b <- length(unique(env_data_melted$variable))

  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18), text = element_text(size = 18),
          plot.title = element_text(size = 16,  face = "bold"))

  ggplot(env_data_melted, aes_(x = ~as.numeric(variable),
                               y = ~as.numeric(temp_row_names),
                               fill = ~Legend)) + geom_tile() +
    scale_fill_gradient(na.value = na.color, low = low_color, high = high_color) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(from = min_year, to = max_year, by = 10)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(from = 50, to = b, by = 50)) +
    xlab("Day of Year") +
    ylab("Year") +
    journal_theme

}

