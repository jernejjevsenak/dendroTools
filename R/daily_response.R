#' daily_response
#'
#' Function calculates all possible values of a selected statistical measure
#' between one or more response variables and daily sequences of environmental
#' data. Calculations are based on moving window which is defined with two
#' arguments: window width and a location in a matrix of daily sequences of
#' environmental data. Window width could be fixed (use fixed_width) or
#' variable width (use lower_limit and upper_limit arguments). In this case,
#' all window widths between lower and upper limit will be used. All calculated
#' measures are stored in a matrix. The location of stored calculated measure
#' in the matrix is indicating a window width (row names) and a location in a
#' matrix of daily sequences of environmental data (column names).
#'
#' @param response a data frame with tree-ring proxy variables as columns and
#' (optional) years as row names. Row.names should be matched with those from a
#' env_data data frame. If not, set row_names_subset = TRUE.
#' @param env_data a data frame of daily sequences of environmental data as
#' columns and (optional) years as row names. Each row represents a year and
#' each column represents a day of a year. Row.names should be matched with
#' those from a response data frame. If not, set row_names_subset = TRUE.
#' @param method a string specifying which method to use. Current possibilities
#' are "cor", "lm" and "brnn".
#' @param measure a string specifying which measure to use. Current
#' possibilities are "r.squared" and "adj.r.squared". If method = "cor",
#' measure is not relevant.
#' @param lower_limit lower limit of window width
#' @param upper_limit upper limit of window width
#' @param fixed_width fixed width used for calculation. If fixed_width is
#' assigned a value, upper_limit and lower_limit will be ignored
#' @param previous_year if set to TRUE, env_data and response variables will be
#' rearranged in a way, that also previous year will be used for calculations of
#' selected statistical measure.
#' @param neurons positive integer that indicates the number of neurons used
#'  for brnn method
#' @param brnn_smooth if set to TRUE, a smoothing algorithm is applied that
#' removes unrealistic calculations which are a result of neural net failure.
#' @param remove_insignificant if set to TRUE, removes all correlations bellow
#' the significant threshold level, based on a selected alpha. For "lm" and
#' "brnn" method, squared threshold is used, which corresponds to R squared
#' statistics.
#' @param alpha significance level used to remove insignificant calculations.
#' @param row_names_subset if set to TRUE, row.names are used to subset
#' env_data and response data frames. Only years from both data frames are
#' kept.
#'
#' @return a list with four elements:
#'   1. calculations is a matrix with all calculated results,
#'   2. method is a string indicating method that was used
#'   3. measure is a string indicating a calculated measure
#'   4. optimized_result is aggregated daily data, that returned the best
#'   calculated measure
#'
#' @export
#'
#' @examples

#' data(daily_temperatures_example)
#' data(example_proxies_1)
#' library(dplyr)
#' oxygen_isotope <- select(example_proxies_1, O)
#' carbon_isotope <- select(example_proxies_1, C)
#'
#' Example1a <- daily_response(response = carbon_isotope,
#' env_data = daily_temperatures_example, method = "lm", measure = "r.squared",
#' lower_limit = 357, upper_limit = 358)
#'
#' \dontrun{
#' Example1b <- daily_response(response = oxygen_isotope,
#' env_data = daily_temperatures_example, method = "lm", measure = "adj.r.squared",
#' lower_limit = 350, upper_limit = 351, remove_insignificant = FALSE)
#'
#' Example1c <- daily_response(response = carbon_isotope,
#' env_data = daily_temperatures_example, method = "lm", measure = "r.squared",
#' lower_limit = 100, upper_limit = 104)
#'
#' Example2 <- daily_response(response = example_proxies_1,
#' env_data = daily_temperatures_example, method = "brnn",
#' measure = "adj.r.squared", fixed_width = 90)
#'
#' Example3 <- daily_response(response = oxygen_isotope,
#' env_data = daily_temperatures_example, method = "cor", lower_limit = 60,
#' upper_limit = 70, remove_insignificant = TRUE)
#'
#' # Example with negative correlations. Data frames are automatically subset.
#' data(example_proxies_2)
#' Example4 <- daily_response(response = example_proxies_2,
#' env_data = daily_temperatures_example, method = "cor",
#' lower_limit = 30, upper_limit = 40, row_names_subset = TRUE)
#' }

daily_response <- function(response, env_data, method = "lm",
                           measure = "r.squared", lower_limit = 30,
                           upper_limit = 270, fixed_width = 0,
                           previous_year = FALSE, neurons = 2,
                           brnn_smooth = TRUE, remove_insignificant = TRUE,
                           alpha = .05, row_names_subset = FALSE) {

  # PART 1 - general data arrangements, warnings abd stops

  set.seed(neurons * 55)

  # Both bojects (response and env_data) are converted to data frames
  response <- data.frame(response)
  env_data <- data.frame(env_data)

  # For measure calculations, both objects need to have the same length,
  # with the exception, when row_names_subset is set to TRUE
  # Stop message in case both data frames do not have the same length
  if (nrow(response) !=  nrow(env_data) & row_names_subset == FALSE)
    stop("Length of env_data and response records differ")

  # Stop message if fixed_width is not between 0 and 365
  if (fixed_width < 0 | fixed_width > 365)
    stop("fixed_width should be between 0 and 365")

  # Stop in case of method == "cor" and ncol(proxies) > 1
  # Correlations could be calculated only for one variable
  if (method == "cor" & ncol(response) > 1)
    stop(paste("More than 1 variable in response data frame not suitable ",
  "for 'cor' method. Use 'lm' or 'brnn'"))

  if (lower_limit >= upper_limit)
    stop("lower_limit can not be higher than upper_limit!")

  if (lower_limit > 365 | lower_limit < 1)
    stop("lower_limit out of bounds! It should be between 1 and 365")

  if (upper_limit > 365 | upper_limit < 1)
    stop("upper_limit out of bounds! It should be between 1 and 365")

  # If row_names_subset == TRUE, data is subseted and ordered based on matching
  # row.names. Additionally, number of characters in row.names is checked.
  # There should be at least three characters (assuming years before 100 will
  # never be analysed, there is no such environmental data avaliable)
  if (row_names_subset == TRUE & nchar(row.names(env_data)[1]) >= 3){

    ncol_response <- ncol(response)

    colnames_response <- colnames(response)

    env_data$year <- row.names(env_data)
    response$year <- row.names(response)

    temporal_data <- merge(response, env_data, by = "year")

    response <- data.frame(temporal_data[, c(2:(1 + ncol_response))],
                           row.names = temporal_data$year)
    colnames(response) <- colnames_response

    env_data <- data.frame(temporal_data[, c((1 + ncol_response + 1):
                                               ncol(temporal_data))],
                           row.names = temporal_data$year)
  }

  # if row.names of env_data and the response data frames are not equal,
  # warning is given.
  if (sum(row.names(env_data) == row.names(response)) != nrow(env_data)) {
    warning("row.names between env_data and response do not match!")
  }

  # If row_names_subset == TRUE, but row.names does not appear to be years,
  # error is given.
  if (row_names_subset == TRUE & nchar(row.names(env_data)[1]) < 3){
    stop(paste("row.names does not appear to be years!",
                "At least three characters needed!"))
  }

  # Data manipulation
  # If use.previous == TRUE, env_data data has to be rearranged accordingly
  if (previous_year == TRUE) {
    response <- response[-nrow(response), ]
    env_data_previous <- env_data[-1, ]
    env_data_current <- env_data[-nrow(env_data), ]
    env_data <- cbind(env_data_previous, env_data_current)
    response <- data.frame(response)
    env_data <- data.frame(env_data)
  }

  # PART 2 - Based on the selected function arguments, different chunks of code
  # will be used. For demonstration:
  # A) Chunks are used if fixed.withd != 0
    # A.1 method == "cor"
    # A.2 method == "lm"
    # A.3 method == "brnn"

  # B) Chunks are used if fixed.withd == 0
    # B.1 method == "cor"
    # B.2 method == "lm"
    # B.3 method == "brnn"

  # A.1 method = "cor"
  if (fixed_width != 0 & method == "cor") {

      # This is an empty matrix, currently filled with NA's
      # Latter, calculations will be stored in this matrix
      temporal_matrix <- matrix(NA, nrow = 1,
        ncol = (ncol(env_data) - fixed_width) + 1)

      # An iterating loop. In each itteration x is calculated and represents
      # response (dependent) variable. X is a moving average. Window width of
      # a moving window is fixed_width. Next, statistical measure is calculated
      # based on a selected method (cor, lm or brnn). Calculation is stored in
      # temporal matrix.
      for (j in 0: (ncol(env_data) - fixed_width)) {
        x <- rowMeans(env_data[1:nrow(env_data),
         (1 + j): (j + fixed_width)], na.rm = TRUE)


        print(paste(j, fixed_width), sep = "")

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)
        temporal_correlation <- cor(response[, 1], x[, 1])

        # Each calculation is printed. Reason: usually it takes several minutes
        # to go through all loops and therefore, users might think that R is
        # not responding. But if each calculation is printed, user could be
        # confident, that R is responding.


        #print (temporal_correlation)
        temporal_matrix[1, j + 1] <- temporal_correlation
      }

     # temporal_matrix is given rownames and colnames. Rownames represent a
     # window width used fot calculations. Colnames represent the position of
     # moving window in a original env_data data frame.
     row.names(temporal_matrix) <- fixed_width
     temporal_colnames <- as.vector(seq(from = 1,
       to = ncol(temporal_matrix), by = 1))
     colnames(temporal_matrix) <- temporal_colnames
  }

  # A.2 method == "lm"
  # For a description see A.1
  if (fixed_width != 0 & method == "lm") {

    temporal_matrix <- matrix(NA, nrow = 1,
      ncol = (ncol(env_data) - fixed_width) + 1)

    for (j in 0:(ncol(env_data) - fixed_width)) {
      x <- rowMeans(env_data[1:nrow(env_data),
        (1 + j) : (j + fixed_width)], na.rm = TRUE)
      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      temporal_df <- data.frame(cbind(x, response))
      temporal_model <- lm(x ~ ., data = temporal_df)
      temporal_summary <- summary(temporal_model)
      temporal_r_squared <- temporal_summary$r.squared
      temporal_adj_r_squared <- temporal_summary$adj.r.squared

      if (measure == "r.squared"){
        temporal_matrix[1, j + 1] <- temporal_r_squared
        print(temporal_r_squared)
      }

      if (measure == "adj.r.squared"){
        temporal_matrix[1, j + 1] <- temporal_adj_r_squared
        print(temporal_adj_r_squared)
      }
    }

    row.names(temporal_matrix) <- fixed_width
    temporal_colnames <- as.vector(seq(from = 1,
      to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
  }

  # A.3 method == "brnn"
  # For a description see A.1
  if (fixed_width != 0 & method == "brnn") {

    temporal_matrix <- matrix(NA, nrow = 1,
      ncol = (ncol(env_data) - fixed_width) + 1)

    for (j in 0: (ncol(env_data) - fixed_width)) {
      x <- rowMeans(env_data[1:nrow(env_data),
        (1 + j): (j + fixed_width)], na.rm = TRUE)
      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      temporal_df <- data.frame(cbind(x, response))
      temporal_model <- try(brnn(x ~ ., data = temporal_df,
                                 neurons = neurons, tol = 1e-6),
                            silent = TRUE)
      temporal_predictions <- try(predict.brnn(temporal_model, temporal_df),
                                  silent = TRUE)

      if (class(temporal_model)[[1]] != "try-error"){

        temporal_r_squared <- 1 - (sum((x[, 1] - temporal_predictions) ^ 2) /
                                     sum((x[, 1] - mean(x[, 1])) ^ 2))
        temporal_adj_r_squared <- 1 - ((1 - temporal_r_squared) *
                                         ((nrow(x) - 1)) /
                                         (nrow(x) -
                                            ncol(as.data.frame(response[, 1]))
                                          -  1))

        if (measure == "r.squared"){
          temporal_matrix[1, j + 1] <- temporal_r_squared
          print(temporal_r_squared)
        }

        if (measure == "adj.r.squared"){
          temporal_matrix[1, j + 1] <- temporal_adj_r_squared
          print(temporal_adj_r_squared)
        } else {

          temporal_matrix[1, j + 1] <- NA
        }
      }
    }

    row.names(temporal_matrix) <- fixed_width
    temporal_colnames <- as.vector(seq(from = 1,
      to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
  }

  # B fixed_width == 0, in this case, lower_limit and upper_limit arguments
  # will be used to define window width of a moving window.
  # B.1 method == "cor"

  if (fixed_width == 0 & method == "cor") {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be stored in this matrix
  temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
    ncol = (ncol(env_data) - lower_limit) + 1)

  # An iterating double loop: 1 outer loop) iterating from lower_limit :
  # upper_limit defines windo.width used for a moving window. 2) inner loop
  # defines the starting position of a moving window.
  # In each itteration, x is calculated and represents a response (dependent)
  # variable. x is a moving average, based on rowMeans function.
  # Next, statistical measure is calculated based on a selected method (cor,
  # lm or brnn). Calculation is stored in temporal matrix in a proper place.
  # The position of stored calculation is informative later used for
  # indiciating optimal values.

  for (k in lower_limit:upper_limit) {

    for (j in 0: (ncol(env_data) - k)) {
      x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + k)], na.rm = T)
      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      temporal_correlation <- cor(response[, 1], x[, 1])
      print(temporal_correlation)
      temporal_matrix[(k - lower_limit) + 1, j + 1] <- temporal_correlation
    }
  }

  # temporal_matrix is given rownames and colnames. Rownames represent a
  # window width used fot calculations. Colnames represent the position of
  # moving window in a original env_data data frame.
  temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit,
    by = 1))
  row.names(temporal_matrix) <- temporal_rownames

  temporal_colnames <- as.vector(seq(from = 1,
    to = ncol(temporal_matrix), by = 1))
  colnames(temporal_matrix) <- temporal_colnames
  }

  # B.2 method == "lm"
  # For a description see B.1
  if (fixed_width == 0 & method == "lm") {

    temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
      ncol = (ncol(env_data) - lower_limit) + 1)

    for (k in lower_limit:upper_limit) {

      for (j in 0: (ncol(env_data) - k)) {
        x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + k)],
          na.rm = T)
        x <- matrix(x, nrow = nrow(env_data), ncol = 1)
        temporal_df <- data.frame(cbind(x, response))
        temporal_model <- lm(x ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        temporal_r_squared <- temporal_summary$r.squared
        temporal_adj_r_squared <- temporal_summary$adj.r.squared

        if (measure == "r.squared"){
          temporal_matrix[(k - lower_limit) + 1, j + 1]  <-
            temporal_r_squared
          print(temporal_r_squared)
        }

        if (measure == "adj.r.squared"){
          temporal_matrix[(k - lower_limit) + 1, j + 1]  <-
            temporal_adj_r_squared
          print(temporal_adj_r_squared)
        }
      }
    }

    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit,
      by = 1))
    row.names(temporal_matrix) <- temporal_rownames

    temporal_colnames <- as.vector(seq(from = 1,
      to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
  }

  # B.3 method == "brnn"
  # For a description see B.1
  if (fixed_width == 0 & method == "brnn") {

    temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
      ncol = (ncol(env_data) - lower_limit) + 1)

    for (k in lower_limit:upper_limit) {

      for (j in 0: (ncol(env_data) - k)) {
        x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + k)],
          na.rm = T)
        x <- matrix(x, nrow = nrow(env_data), ncol = 1)
        temporal_df <- data.frame(cbind(x, response))
        temporal_model <- try(brnn(x ~ ., data = temporal_df, neurons = neurons,
                                   tol = 1e-6), silent = TRUE)
        temporal_predictions <- try(predict.brnn(temporal_model, temporal_df),
                                    silent = TRUE)

        if (class(temporal_model)[[1]] != "try-error"){

          temporal_r_squared <- 1 - (sum((x[, 1] - temporal_predictions) ^ 2) /
                                       sum((x[, 1] - mean(x[, 1])) ^ 2))
          temporal_adj_r_squared <- 1 - ((1 - temporal_r_squared) *
                                           ((nrow(x) - 1)) /
                                           (nrow(x) -
                                              ncol(as.data.frame(response[, 1]))
                                            - 1))

          if (measure == "r.squared"){
            temporal_matrix[(k - lower_limit) + 1, j + 1]  <- temporal_r_squared
            print(temporal_r_squared)
          }

          if (measure == "adj.r.squared"){
            temporal_matrix[(k - lower_limit) + 1, j + 1]  <-
              temporal_adj_r_squared
            print(temporal_adj_r_squared)
          }

        } else {
          temporal_matrix[(k - lower_limit) + 1, j + 1] <- NA

        }
      }
    }

    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit,
      by = 1))
    row.names(temporal_matrix) <- temporal_rownames

    temporal_colnames <- as.vector(seq(from = 1,
      to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
  }

  # PART 3: smoothing function, if brnn method is used. It turnes out, that
  # brnn sometimes (1 - 3 % of calculations) fails to construct a realistic
  # result. In those cases, much lower r.squared (or adj.r.squared) are
  # calculated. smooth_matrix function removes unrealistic calculations and
  # replace them with an average of values in a window 3 x 3. Maximum value
  # is not affected by any means.
  # [i - 1, j - 1], [i - 1, j], [i - 1, j + 1]
  # [    i, j - 1], [    i, j], [    i, j + 1]
  # [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]
  if (method == "brnn" & brnn_smooth == TRUE){
    temporal_matrix <- smooth_matrix(temporal_matrix, factor_drop = 0.7,
      repeats = 2)
  }

  # To enhance the visualisation, insignificant values
  # are removed if remove_insignificant == TRUE
  if (remove_insignificant == TRUE){
    critical_threshold_cor <- critical_r(nrow(response), alpha = alpha)
    critical_threshold_cor2 <- critical_threshold_cor ^ 2

    # With the following chunk, we need to cosider negative correlations
    # and threat them separately
    overall_max <- max(temporal_matrix, na.rm = TRUE)
    overall_min <- min(temporal_matrix, na.rm = TRUE)

    # 1 positive correlations
    if (method == "cor" & (abs(overall_max) > abs(overall_min))) {
      temporal_matrix[abs(temporal_matrix) < critical_threshold_cor] <- NA
    # 2 negative correlations
    } else if (method == "cor" & (abs(overall_max) < abs(overall_min))) {
      temporal_matrix[temporal_matrix > -critical_threshold_cor] <- NA
    # 3 lm and brnn method
      } else if (method == "lm" | method == "brnn") {
      temporal_matrix[temporal_matrix < critical_threshold_cor2] <- NA
      }
  }

  # PART 4: Final list is being created and returned as a function output
  # When metohod == "cor", different final_list is created
  if (method == "lm" | method == "brnn") {
    final_list <- list(calculations = temporal_matrix, method = method,
      measure = measure)
  }

  if (method == "cor"){
    final_list <- list(calculations = temporal_matrix, method = method,
                        measure = method)
  }


  # Here we add additional, fourth element: the optimal sequence of days
  # that returns the best selected statistical measure

  # In case of negative correlations, different strategy is applied.
  # For more detailed description see plot_extreme()

  overall_max <- max(temporal_matrix, na.rm = TRUE)
  overall_min <- min(temporal_matrix, na.rm = TRUE)

  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  # There are unimportant warnings produced:
  # no non-missing arguments to max; returning -Inf


  if ((abs(overall_max) > abs(overall_min)) == TRUE) {

    # maximum value is located. Row indeces are needed to query information
    # about the window width used to calculate the maximum. Column name is
    # needed to query the starting day.
    max_result <- suppressWarnings(which.max(apply(temporal_matrix,
                                                   MARGIN = 2, max,
                                                   na.rm = TRUE)))
    plot_column <- max_result
    max_index <- which.max(temporal_matrix[, names(max_result)])
    row_index <- row.names(temporal_matrix)[max_index]
  }

  if ((abs(overall_max) < abs(overall_min)) == TRUE) {

    min_result <- suppressWarnings(which.min(apply(temporal_matrix,
                                                   MARGIN = 2, min,
                                                   na.rm = TRUE)))
    plot_column <- min_result
    min_index <- which.min(temporal_matrix[, names(min_result)])
    row_index <- row.names(temporal_matrix)[min_index]
  }

  # The fourth return element is being created: rowMeans of optimal sequence:
  dataf <- data.frame(rowMeans(env_data[, as.numeric(plot_column):
                                         (as.numeric(plot_column) +
                                            as.numeric(row_index) - 1)],
                              na.rm = TRUE))

  colnames(dataf) <- "Optimized rowNames"

  final_list[[4]] <- dataf

  # Additional check:
  if (method == "lm" & measure == "r.squared"){
    temporal_df <- data.frame(cbind(dataf, response))
    temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
    temporal_summary <- summary(temporal_model)
    optimized_result <- temporal_summary$r.squared
  }

  if (method == "lm" & measure == "adj.r.squared"){
    temporal_df <- data.frame(cbind(dataf, response))
    temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
    temporal_summary <- summary(temporal_model)
    optimized_result <- temporal_summary$adj.r.squared
  }

  if (method == "brnn" & measure == "r.squared"){
    temporal_df <- data.frame(cbind(dataf, response))
    temporal_model <- brnn(Optimized.rowNames ~ ., data = temporal_df,
                           neurons = neurons, tol = 1e-6)
    temporal_predictions <- try(predict.brnn(temporal_model,
                                             temporal_df), silent = TRUE)
    optimized_result <- 1 - (sum((temporal_df[, 1] -
                                    temporal_predictions) ^ 2) /
                                 sum((temporal_df[, 1] -
                                        mean(temporal_df[, 1])) ^ 2))
  }

  if (method == "brnn" & measure == "adj.r.squared"){
    temporal_df <- data.frame(cbind(dataf, response))
    temporal_model <- brnn(Optimized.rowNames ~ .,
                           data = temporal_df, neurons = neurons, tol = 1e-6)
    temporal_predictions <- try(predict.brnn(temporal_model, temporal_df),
                                silent = TRUE)
    temporal_r_squared <- 1 - (sum((temporal_df[, 1] -
                                      temporal_predictions) ^ 2) /
                               sum((temporal_df[, 1] -
                                      mean(temporal_df[, 1])) ^ 2))
    optimized_result <- 1 - ((1 - temporal_r_squared) *
                                     ((nrow(temporal_df) - 1)) /
                               (nrow(temporal_df) -
                                  ncol(as.data.frame(response[, 1])) - 1))
  }

  if (method == "cor"){
    optimized_result <- cor(dataf, response)
  }

  return(final_list)
}
