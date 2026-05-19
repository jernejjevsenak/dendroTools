#' daily_response
#'
#' Function calculates all possible values of a selected statistical metric
#' between one or more response variables and daily sequences of environmental
#' data. Calculations are based on moving window which is defined with two
#' arguments: window width and a location in a matrix of daily sequences of
#' environmental data. Window width could be fixed (use fixed_width) or
#' variable width (use lower_limit and upper_limit arguments). In this case,
#' all window widths between lower and upper limit will be used. All calculated
#' metrics are stored in a matrix. The location of stored calculated metric
#' in the matrix is indicating a window width (row names) and a location in a
#' matrix of daily sequences of environmental data (column names).
#'
#' @param response a data frame with tree-ring proxy variables as columns and
#' (optional) years as row names. Row.names should be matched with those from a
#' env_data data frame. If not, set row_names_subset = TRUE.
#' @param env_data a data frame of daily sequences of environmental data as
#' columns and years as row names. Each row represents a year and
#' each column represents a day of a year. Row.names should be matched with
#' those from a response data frame. If not, set row_names_subset = TRUE.
#' Alternatively, env_data could be a tidy data with three columns,
#' i.e. Year, DOY and third column representing values of mean temperatures,
#' sum of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data to TRUE.
#' @param method a character string specifying which method to use. Current
#' possibilities are "cor" (default), "lm" and "brnn".
#' @param cor_method a character string indicating which correlation
#' coefficient is to be computed. One of "pearson" (default), "kendall", or
#' "spearman".
#' @param metric a character string specifying which metric to use. Current
#' possibilities are "r.squared" and "adj.r.squared". If method = "cor",
#' metric is not relevant.
#' @param lower_limit lower limit of window width
#' @param upper_limit upper limit of window width
#' @param fixed_width fixed width used for calculation. If fixed_width is
#' assigned a value, upper_limit and lower_limit will be ignored
#' @param previous_year logical. If TRUE, previous-year climate data are included
#' in the analysis. If FALSE, no previous-year climate data are included and
#' number_previous_years is ignored.
#' @param number_previous_years integer between 1 and 5 specifying how many
#' previous years should be included in the environmental matrix when
#' previous_year = TRUE. For example, number_previous_years = 2 uses climate
#' data from years t - 2, t - 1 and t for response year t. If NULL and
#' previous_year = TRUE, one previous year is included.
#' @param neurons positive integer that indicates the number of neurons used
#'  for brnn method
#' @param brnn_smooth if set to TRUE, a smoothing algorithm is applied that
#' removes unrealistic calculations which are a result of neural net failure.
#' @param remove_insignificant if set to TRUE, removes all correlations bellow
#' the significant threshold level, based on a selected alpha. For "lm" and
#' "brnn" method, squared correlation is used as a threshold
#' @param alpha significance level used to remove insignificant calculations.
#' @param row_names_subset if set to TRUE, row.names are used to subset
#' env_data and response data frames. Only years from both data frames are
#' kept.
#' @param aggregate_function character string specifying how the daily data
#' should be aggregated. The default is 'mean', the other options are 'median',
#' 'sum', 'min', 'max' and 'quantile'.
#' @param quantile_prob numeric value between 0 and 1 specifying the quantile
#' probability used when aggregate_function = 'quantile'. For example,
#' quantile_prob = 0.95 calculates the 95th percentile. The default is 0.5.
#' @param temporal_stability_check character string, specifying, how temporal stability
#' between the optimal selection and response variable(s) will be analysed. Current
#' possibilities are "sequential", "progressive" and "running_window". Sequential check
#' will split data into k splits and calculate selected metric for each split. Progressive
#' check will split data into k splits, calculate metric for the first split and then
#' progressively add 1 split at a time and calculate selected metric. For running window,
#' select the length of running window with the k_running_window argument.
#' @param k_running_window the length of running window for temporal stability check.
#' Applicable only if temporal_stability argument is set to running window.
#' @param k integer, number of breaks (splits) for temporal stability and cross validation
#' analysis.
#' @param cross_validation_type character string, specifying, how to perform cross validation
#' between the optimal selection and response variables. If the argument is set to "blocked",
#' years will not be shuffled. If the argument is set to "randomized", years will be shuffled.
#' @param subset_years a subset of years to be analyzed. Should be given in the form of
#' subset_years = c(1980, 2005)
#' @param ylimits limit of the y axes for plot_extreme. It should be given in
#' the form of: ylimits = c(0,1)
#' @param seed optional seed argument for reproducible results
#' @param tidy_env_data if set to TRUE, env_data should be inserted as a data frame with three
#' columns: "Year", "DOY", "Precipitation/Temperature/etc."
#' @param reference_window character string, the reference_window argument describes,
#' how each calculation is referred. There are three different options: 'start'
#' (default), 'end' and 'middle'. If the reference_window argument is set to 'start',
#' then each calculation is related to the starting day of window. If the
#' reference_window argument is set to 'middle', each calculation is related to the
#' middle day of window calculation. If the reference_window argument is set to
#' 'end', then each calculation is related to the ending day of window calculation.
#' For example, if we consider correlations with window from DOY 15 to DOY 35. If
#' reference window is set to  'start', then this calculation will be related to the
#' DOY 15. If the reference window is set to 'end', then this calculation will be
#' related to the DOY 35. If the reference_window is set to 'middle', then this
#' calculation is related to DOY 25.
#' The optimal selection, which describes the optimal consecutive days that returns
#' the highest calculated metric and is obtained by the $plot_extreme output, is the
#' same for all three reference windows.
#' @param boot logical, if TRUE, bootstrap procedure will be used to calculate
#' estimates correlation coefficients, R squared or adjusted R squared metrices
#' @param boot_n The number of bootstrap replicates
#' @param boot_ci_type A character string representing the type of bootstrap intervals
#' required. The value should be any subset of the values c("norm","basic", "stud",
#' "perc", "bca").
#' @param boot_conf_int A scalar or vector containing the confidence level(s) of
#' the required interval(s)
#' @param day_interval a vector of two values defining the interval of days used
#' for calculations. Positive values indicate days in the current year. Negative
#' values indicate days in the previous-year block and are used only when
#' previous_year = TRUE. If previous_year = FALSE and negative values are
#' supplied, day_interval is ignored, a warning is issued, and the analysis is
#' performed for the current year only using day_interval = c(1, 366). If
#' number_previous_years > 1, the previous-year block starts with the earliest
#' included previous year. For example, previous_year = TRUE,
#' number_previous_years = 2 and day_interval = c(-1, 366) analyses the full
#' sequence from DOY 1 of year t - 2 to DOY 366 of year t.
#' @param dc_method a character string to determine the method to detrend climate
#' data.  Possible values are "none" (default) and "SLD" which refers to Simple
#' Linear Detrending
#' @param cor_na_use an optional character string giving a method for computing
#' covariances in the presence of missing values for correlation coefficients.
#' This must be (an abbreviation of) one of the strings "everything" (default),
#' "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". See
#' also the documentation for the base cor() function.
#' @param skip_window_length an integer specifying the frequency of window
#' selection for the calculations of climate-growth relationships. The default
#' value is 1, indicating that every window is included in the calculations.
#' When set to a value greater than 1, the function selectively processes
#' windows at regular intervals defined by this parameter. For instance, if
#' skip_window_length = 2, the function processes every second window.
#' Similarly, if skip_window_length = 3, every third window is processed,
#' skipping two windows in between each selected one. This parameter allows for
#' controlling the granularity of the analysis and can help in reducing
#' computation time by focusing on a subset of the data.
#' @param skip_window_position an integer specifying the frequency of window
#' positions used in the calculations of climate-growth relationships. The
#' default value is 1, indicating that every window position is included in the
#' calculations. When set to a value greater than 1, the function selectively
#' processes window positions at regular intervals defined by this parameter.
#' For instance, if skip_window_position = 2, the function processes every
#' second window position. Similarly, if skip_window_position = 3, every third
#' window position is processed, skipping two positions in between each selected
#' one. This parameter allows for controlling the granularity of the analysis
#' and can help in reducing computation time by focusing on a subset of the data.
#'
#' @return a list with 19 elements:
#' \enumerate{
#'  \item $calculations - a matrix with calculated metrics
#'  \item $method - the character string of a method
#'  \item $metric - the character string indicating the metric used for calculations
#'  \item $analysed_period - the character string specifying the analysed period based on the information from row names. If there are no row names, this argument is given as NA
#'  \item $optimized_return - data frame with two columns, response variable and aggregated (averaged) daily data that return the optimal results. This data.frame could be directly used to calibrate a model for climate reconstruction
#'  \item $optimized_return_all - a data frame with aggregated daily data, that returned the optimal result for the entire env_data (and not only subset of analysed years)
#'  \item $transfer_function - a ggplot object: scatter plot of optimized return and a transfer line of the selected method
#'  \item $temporal_stability -  a data frame with calculations of selected metric for different temporal subsets
#'  \item $cross_validation - a data frame with cross validation results
#'  \item $plot_heatmap - ggplot2 object: a heatmap of calculated metrics
#'  \item $plot_extreme - ggplot2 object: line plot of a row with the highest value in a matrix of calculated metrics
#'  \item $type - the character string describing type of analysis: daily or monthly
#'  \item $reference_window - character string, which reference window was used for calculations
#'  \item $boot_lower - matrix with lower limit of confidence intervals of bootstrap calculations
#'  \item $boot_upper - matrix with upper limit of confidence intervals of bootstrap calculations
#'  \item $aggregated_climate - matrix with all aggregated climate series
#'  \item $previous_year - logical indicating whether previous-year climate data were used
#'  \item $number_previous_years - integer indicating how many previous years were used
#'}
#' @export
#'
#' @examples
#' \donttest{
#'
#' # The examples below are enclosed within donttest{} to minimize the execution
#' # time during R package checks. Additionally, all examples include the
#' # parameters `skip_window_length` and `skip_window_position`, which limit the
#' # number of combinations evaluated in climate-growth correlation calculations.
#' # To explore all possible combinations, users should set both parameters to 1.
#'
#' # Load the dendroTools R package
#' library(dendroTools)
#'
#' # Load data
#' data(data_MVA)
#' data(data_TRW)
#' data(data_TRW_1)
#' data(example_proxies_individual)
#' data(example_proxies_1)
#' data(LJ_daily_temperatures)
#'
#' example_basic <- daily_response(response = data_MVA,
#'                           env_data = LJ_daily_temperatures,
#'                           row_names_subset = TRUE,
#'                           fixed_width = 25,
#'                           lower_limit = 35, upper_limit = 45,
#'                           remove_insignificant = FALSE,
#'                           aggregate_function = 'median',
#'                           alpha = 0.05, cor_method = "spearman",
#'                           previous_year = FALSE, boot = TRUE,
#'                           boot_n = 10,
#'                           skip_window_length = 50,
#'                           skip_window_position = 50,
#'                           reference_window = "end", k = 5,
#'                           dc_method = "SLD",
#'                           day_interval = c(1, 250))
#'
#' # 1 Example with fixed width. Lower and upper limits are ignored.
#' example_daily_response <- daily_response(response = data_MVA,
#'     env_data = LJ_daily_temperatures,
#'     method = "cor", fixed_width = 40, cor_method = "spearman",
#'     row_names_subset = TRUE, previous_year = TRUE,
#'     remove_insignificant = TRUE, boot = TRUE,
#'     alpha = 0.005, aggregate_function = 'mean',
#'     day_interval = c(-100, 250), skip_window_length = 100,
#'     reference_window = "start", skip_window_position = 100)
#'
#' # summary(example_daily_response)
#' # plot(example_daily_response, type = 1)
#' # plot(example_daily_response, type = 2)
#'
#' # 2 Example for past and present. Use subset_years argument.
#' example_MVA_early <- daily_response(response = data_MVA,
#'     env_data = LJ_daily_temperatures, cor_method = "kendall",
#'     method = "lm", lower_limit = 21, upper_limit = 91,
#'     row_names_subset = TRUE, previous_year = TRUE,
#'     remove_insignificant = TRUE, alpha = 0.05,
#'     subset_years = c(1940, 1980),
#'     fixed_width = 45,
#'     aggregate_function = 'sum',
#'     skip_window_length = 50,
#'     skip_window_position = 50)
#'
#' example_MVA_late <- daily_response(response = data_MVA,
#'     env_data = LJ_daily_temperatures,
#'     method = "cor", lower_limit = 21, upper_limit = 60,
#'     row_names_subset = TRUE, previous_year = TRUE,
#'     remove_insignificant = TRUE, alpha = 0.05,
#'     subset_years = c(1981, 2010),
#'     skip_window_length = 50,
#'     skip_window_position = 50)
#'
#' # plot(example_MVA_early, type = 1)
#' # plot(example_MVA_late, type = 1)
#' # plot(example_MVA_early, type = 2)
#' # plot(example_MVA_late, type = 2)
#'
#' # 3 Example with negative correlations
#' example_neg_cor <- daily_response(response = data_TRW_1,
#'     env_data = LJ_daily_temperatures, previous_year = TRUE,
#'     method = "cor", lower_limit = 21, upper_limit = 90,
#'     row_names_subset = TRUE, remove_insignificant = TRUE,
#'     alpha = 0.05, skip_window_length = 50,
#'     skip_window_position = 50)
#'
#' # summary(example_neg_cor)
#' # plot(example_neg_cor, type = 1)
#' # plot(example_neg_cor, type = 2)
#'
#' # 4 Example of multiproxy analysis
#' # summary(example_proxies_1)
#' # cor(example_proxies_1)
#'
#' example_multiproxy <- daily_response(response = example_proxies_1,
#'    env_data = LJ_daily_temperatures,
#'    method = "lm", metric = "adj.r.squared",
#'    lower_limit = 21, upper_limit = 180,
#'    row_names_subset = TRUE, previous_year = FALSE,
#'    remove_insignificant = TRUE, alpha = 0.05,
#'    skip_window_length = 50,
#'    skip_window_position = 50)
#'
#' # plot(example_multiproxy, type = 1)
#'
#' # 5 Example to test the temporal stability
#' example_MVA_ts <- daily_response(response = data_MVA,
#'    env_data = LJ_daily_temperatures, method = "brnn",
#'    lower_limit = 100, metric = "adj.r.squared", upper_limit = 180,
#'    row_names_subset = TRUE, remove_insignificant = TRUE, alpha = 0.05,
#'    temporal_stability_check = "running_window", k_running_window = 10,
#'    skip_window_length = 50, skip_window_position = 50)
#'
#' # Check the results for temporal stability
#' # example_MVA_ts$temporal_stability
#'
#' # 6 Example with nonlinear brnn estimation
#' example_brnn <- daily_response(response = data_MVA,
#'    env_data = LJ_daily_temperatures, method = "brnn", boot = FALSE,
#'    lower_limit = 100, metric = "adj.r.squared", upper_limit = 101,
#'    row_names_subset = TRUE, remove_insignificant = TRUE, boot_n = 10,
#'    skip_window_length = 50, skip_window_position = 50)
#'
#' # summary(example_brnn)
#'
#' # Example using quantiles for the aggregation
#' example_q95 <- daily_response(
#' response = data_MVA,
#' env_data = LJ_daily_temperatures,
#' method = "cor",
#' fixed_width = 30,
#' row_names_subset = TRUE,
#' previous_year = TRUE,
#' aggregate_function = "quantile",
#' quantile_prob = 0.95
#' )
#'
#' # Example using two previous years plus the current year
#' example_two_previous_years <- daily_response(
#' response = data_MVA,
#' env_data = LJ_daily_temperatures,
#' method = "cor",
#' fixed_width = 60,
#' row_names_subset = TRUE,
#' previous_year = TRUE,
#' number_previous_years = 2,
#' day_interval = c(-1, 250),
#' skip_window_length = 50,
#' skip_window_position = 50
#' )
#'
#' # Example with previous-previous year effects
#' example_basic <- daily_response(response = data_MVA,
#'                           env_data = LJ_daily_temperatures,
#'                           row_names_subset = TRUE,
#'                           reference_window = "end",
#'                           lower_limit = 35, upper_limit = 75,
#'                           previous_year = TRUE,
#'                           remove_insignificant = TRUE,
#'                           number_previous_years = 3)
#' plot_heatmap(example_basic)
#' plot_extreme(example_basic)
#' summary(example_basic)
#' }

daily_response <- function(response, env_data, method = "cor",
                           metric = "r.squared", cor_method = "pearson",
                           lower_limit = 30, upper_limit = 90, fixed_width = 0,
                           previous_year = FALSE, number_previous_years = NULL,
                           neurons = 1,
                           brnn_smooth = TRUE, remove_insignificant = FALSE,
                           alpha = .05, row_names_subset = FALSE,
                           aggregate_function = 'mean',
                           quantile_prob = 0.5,
                           temporal_stability_check = "sequential", k = 2,
                           k_running_window = 30, cross_validation_type = "blocked",
                           subset_years = NULL,
                           ylimits = NULL, seed = NULL, tidy_env_data = FALSE,
                           reference_window = 'start',  boot = FALSE, boot_n = 1000,
                           boot_ci_type = "norm", boot_conf_int = 0.95,
                           day_interval = NULL,
                           dc_method = NULL,
                           cor_na_use = "everything",
                           skip_window_length = 1,
                           skip_window_position = 1
) {

  ##############################################################################
  # 1 day interval is organized

  days_per_year <- 366L

  # previous_year is the master switch for previous-year analyses.
  # If previous_year = FALSE, number_previous_years is ignored.
  if (!is.logical(previous_year) ||
      length(previous_year) != 1 ||
      is.na(previous_year)) {
    stop("previous_year must be either TRUE or FALSE.")
  }

  if (isFALSE(previous_year)) {

    if (!is.null(number_previous_years) &&
        is.numeric(number_previous_years) &&
        length(number_previous_years) == 1 &&
        !is.na(number_previous_years) &&
        number_previous_years > 0) {
      warning(paste0("number_previous_years is ignored because ",
                     "previous_year = FALSE."))
    }

    number_previous_years <- 0L

  } else {

    # Backward compatibility:
    # previous_year = TRUE still means one previous year unless
    # number_previous_years is explicitly supplied.
    if (is.null(number_previous_years)) {
      number_previous_years <- 1L
    }

    if (!is.numeric(number_previous_years) ||
        length(number_previous_years) != 1 ||
        is.na(number_previous_years) ||
        number_previous_years %% 1 != 0 ||
        number_previous_years < 1 ||
        number_previous_years > 5) {
      stop(paste0("When previous_year = TRUE, number_previous_years must be ",
                  "a single integer between 1 and 5."))
    }

    number_previous_years <- as.integer(number_previous_years)
  }

  year_block_width <- days_per_year * (number_previous_years + 1L)

  # Default interval:
  # current year only: c(1, 366)
  # previous-year analysis: c(-1, 366), i.e. from DOY 1 of the earliest
  # included previous year to DOY 366 of the current year.
  if (is.null(day_interval)) {
    day_interval <- if (isTRUE(previous_year)) c(-1, days_per_year) else c(1, days_per_year)
  }

  if (!is.numeric(day_interval) ||
      length(day_interval) != 2 ||
      any(is.na(day_interval))) {
    stop("day_interval must be a numeric vector of length 2.")
  }

  # Negative values in day_interval imply previous-year climate data.
  # If previous_year = FALSE, previous-year analyses are explicitly switched off,
  # so day_interval is ignored and replaced by the full current-year interval.
  if (isFALSE(previous_year) && any(day_interval < 0)) {

    day_interval <- c(1, days_per_year)

    warning(paste0("Negative values were supplied in day_interval, but ",
                   "previous_year = FALSE. The day_interval argument is ignored ",
                   "and the analysis is performed for the current year only ",
                   "using day_interval = c(1, 366)."))
  }

  offset_start <- day_interval[1]
  offset_end <- day_interval[2]

  if (offset_start == 0 || offset_end == 0) {
    stop("day_interval cannot contain 0. Use negative values for previous-year days and positive values for current-year days.")
  }

  # If both limits are positive while previous_year = TRUE, previous years are
  # not actually included in the selected interval. In that case, keep the
  # analysis current-year only.
  if (offset_start > 0 && offset_end > 0 && isTRUE(previous_year)) {

    number_previous_years <- 0L
    previous_year <- FALSE
    year_block_width <- days_per_year

    warning(paste0("Previous years are not included in selected day_interval. ",
                   "The analysis is performed for the current year only."))
  }

  # Convert the user-facing day_interval to column positions in the lagged matrix.
  # For number_previous_years = 2, the internal matrix is:
  # columns 1:366     = year t - 2
  # columns 367:732   = year t - 1
  # columns 733:1098  = year t
  if (offset_start < 0 && offset_end < 0) {

    offset_start <- abs(offset_start)
    offset_end <- abs(offset_end)

  } else if (offset_start < 0 && offset_end > 0) {

    offset_start <- abs(offset_start)
    offset_end <- offset_end + days_per_year * number_previous_years

  } else if (offset_start > 0 && offset_end > 0) {

    # Current-year-only interval; already in current-year coordinates.
    offset_start <- offset_start
    offset_end <- offset_end

  } else {

    stop("day_interval must not run from current-year days back to previous-year days.")

  }

  if (offset_start > offset_end) {
    stop("day_interval is invalid after conversion. The start of the interval is after the end.")
  }

  if (offset_start < 1 || offset_end > year_block_width) {
    stop(paste0("day_interval is outside the available climate sequence. ",
                "For number_previous_years = ", number_previous_years,
                ", the available width is ", year_block_width, " days."))
  }

  # Calculate the max_window allowed
  max_window <- offset_end - offset_start + 1

  # If max_window is smaller than upper_limit, upper_limit must be reduced
  if (upper_limit > max_window) {

    upper_limit <- max_window

    if (fixed_width == 0) {
      warning(paste0("The upper_limit is outside your day_interval and",
                     " therefore reduced to the maximum allowed: ", max_window, "."))
    }
  }

  if (lower_limit > max_window) {

    lower_limit <- max_window

    if (fixed_width == 0) {
      warning(paste0("The lower_limit is outside your day_interval and",
                     " therefore reduced to the minimum allowed: ", max_window, "."))
    }
  }

  # Also correction for fixed-window approach
  if (fixed_width > max_window) {

    stop(paste0("The selected fixed_width is outside your day_interval. ",
                "Decrease the fixed_width argument to at least: ", max_window, "."))
  }

  # offset_end is converted to the number of days that must be skipped at the
  # end of the full environmental matrix.
  offset_end <- year_block_width - offset_end

  ##############################################################################

  if (fixed_width != 0){
    lower_limit = 30
    upper_limit = 200
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Defining global variables
  median <- NULL
  proxy <- NULL
  optimized_return <- NULL
  transfer_f <- NULL
  journal_theme <- NULL
  CV <- NULL
  Period <- NULL
  Years <- NULL
  yearABC <- NULL
  RMSE <- NULL
  RE <- NULL
  CE <- NULL
  DE <- NULL
  d <- NULL

  temporal_matrix_lower <- NULL
  temporal_matrix_upper <- NULL

  # Check selected aggregation function
  allowed_aggregate_functions <- c("mean", "median", "sum", "min", "max", "quantile")

  if (!(aggregate_function %in% allowed_aggregate_functions)) {
    stop(paste0(
      "aggregate_function is '", aggregate_function,
      "'. Instead it should be one of: ",
      paste(allowed_aggregate_functions, collapse = ", "), "."
    ))
  }

  # Check quantile probability
  if (!is.numeric(quantile_prob) ||
      length(quantile_prob) != 1 ||
      is.na(quantile_prob) ||
      quantile_prob < 0 ||
      quantile_prob > 1) {
    stop("quantile_prob must be a single numeric value between 0 and 1.")
  }

  # Internal helper for aggregating daily data across rows
  aggregate_daily_window <- function(x) {

    x <- data.frame(x)

    if (aggregate_function == "mean") {
      return(rowMeans(x, na.rm = TRUE))
    }

    if (aggregate_function == "median") {
      return(apply(x, 1, median, na.rm = TRUE))
    }

    if (aggregate_function == "sum") {
      return(apply(x, 1, sum, na.rm = TRUE))
    }

    if (aggregate_function == "min") {
      return(apply(x, 1, min, na.rm = TRUE))
    }

    if (aggregate_function == "max") {
      return(apply(x, 1, max, na.rm = TRUE))
    }

    if (aggregate_function == "quantile") {
      return(apply(
        x, 1, quantile,
        probs = quantile_prob,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      ))
    }
  }

  # Internal helper for constructing a multi-year climate matrix.
  # Rows are current response years. Columns are ordered from the oldest previous
  # year to the current year.
  build_lagged_env_data <- function(x, n_previous_years) {

    x <- data.frame(x, check.names = FALSE)

    if (n_previous_years == 0L) {
      return(x)
    }

    env_years <- suppressWarnings(as.integer(row.names(x)))

    if (any(is.na(env_years))) {
      stop("For previous-year analyses, row.names of env_data must be calendar years.")
    }

    if (any(duplicated(env_years))) {
      stop("Duplicated years are present in row.names(env_data).")
    }

    row.names(x) <- as.character(env_years)
    x <- x[order(env_years), , drop = FALSE]
    env_years <- as.integer(row.names(x))

    current_years <- env_years[
      vapply(env_years, function(y) {
        all((y - n_previous_years):y %in% env_years)
      }, logical(1))
    ]

    if (length(current_years) == 0L) {
      stop(paste0("No years in env_data have the required ",
                  n_previous_years, " previous year(s)."))
    }

    lag_blocks <- lapply(seq(n_previous_years, 0L), function(lag_i) {

      block <- x[as.character(current_years - lag_i), , drop = FALSE]
      row.names(block) <- as.character(current_years)
      colnames(block) <- paste0("Y", -lag_i, "_", colnames(x))

      block
    })

    out <- do.call(cbind, lag_blocks)
    row.names(out) <- as.character(current_years)

    data.frame(out, check.names = FALSE)
  }

  # If there is a column name samp.depth in response data frame, warning is given
  if ("samp.depth" %in% colnames(response)){

    samp.depth_index <- grep("samp.depth", colnames(response))
    response <- response[, -samp.depth_index,F]

    warning("Removed the samp.depth from response data frame")
  }

  # If there is more than 2 columns in response data frame, give a warning
  if (ncol(response) > 1){
    warning(paste0("Your response data frame has more than 1 column! Are you doing a multiproxy research?",
                   " If so, OK. If not, check your response data frame!"))
  }


  # If env_data is given in tidy version, transformation is needed
  if (tidy_env_data == TRUE){

    n_col_tidy_DF <- ncol(env_data)
    colnames_tidy_DF <- colnames(env_data)

    if (ncol(env_data) != 3){
      stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                 "env_data should have 3 columns, but it has", n_col_tidy_DF, "instead!"))
    }

    if (colnames_tidy_DF[1] != "Year"){
      stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                 "The first column name of the env_data should be 'Year', but it is",
                 colnames_tidy_DF[1], "instead!"))
    }

    if (colnames_tidy_DF[2] != "DOY"){
      stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                 "The second column name of the env_data should be 'DOY', but it is",
                 colnames_tidy_DF[2], "instead!"))
    }

    value_variable = colnames(env_data)[3]
    env_data <- dcast(env_data, Year~DOY, value.var = value_variable)
    env_data <- years_to_rownames(env_data, "Year")



  }

  # PART 1 - general data arrangements, warnings and stops
  # Both bojects (response and env_data) are converted to data frames
  response <- data.frame(response)
  env_data <- data.frame(env_data)

  # Here we save the original response data that will be used later
  response_original <- response

  # Previous-year analyses require alignment by calendar year names. This avoids
  # accidental row-position mismatches after adding lagged climate years.
  if (number_previous_years > 0L && row_names_subset == FALSE) {
    row_names_subset <- TRUE
    warning(paste0("Previous-year analyses require alignment by year names. ",
                   "row_names_subset is set to TRUE."))
  }

  # Build current-year or multi-year environmental matrix.
  # For number_previous_years = 0 this returns env_data unchanged.
  env_data <- build_lagged_env_data(env_data, number_previous_years)

  # This object is used later for optimized_return_all. It should represent the
  # full lagged climate matrix before subset_years is applied.
  env_data_original <- env_data

  # For metric calculations, both objects need to have the same length,
  # with the exception, when row_names_subset is set to TRUE
  # Stop message in case both data frames do not have the same length
  if (nrow(response) != nrow(env_data) & row_names_subset == FALSE)
    stop("Length of env_data and response records differ")

  # Stop in case of method == "cor" and ncol(proxies) > 1
  # Correlations could be calculated only for one variable
  if (method == "cor" & ncol(response) > 1)
    stop(paste("More than 1 variable in response data frame not suitable ",
               "for 'cor' method. Use 'lm' or 'brnn'"))

  #######################################################
  # Rules for selected window limits

  if (fixed_width < 0 | fixed_width > year_block_width)
    stop(paste0("fixed_width should be between 0 and ", year_block_width))

  if (lower_limit > upper_limit)
    stop("lower_limit can not be higher than upper_limit!")

  if (lower_limit > year_block_width | lower_limit < 1)
    stop(paste0("lower_limit out of bounds! It should be between 1 and ",
                year_block_width))

  if (upper_limit > year_block_width | upper_limit < 1)
    stop(paste0("upper_limit out of bounds! It should be between 1 and ",
                year_block_width))

  # Make sure the selected method is appropriate
  if (!is.null(dc_method)){

    if (!(dc_method %in% c("SLD"))){

      stop(paste0('dc_method should be SLD but instead it is:',dc_method))

    }
  }



  # Warn users in case of missing values. The original threshold was 270
  # missing days for one 366-day climate year, so it is scaled when multiple
  # previous years are included.
  missing_day_threshold <- 270 * (number_previous_years + 1L)
  env_temp <- env_data[row.names(env_data) %in% row.names(response),]

  # Subset of years
  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    subset_seq <- seq(lower_subset, upper_subset)
    env_temp <- subset(env_temp, row.names(env_temp) %in% subset_seq)
  }

  na_problem <- data.frame(na_sum = rowSums(is.na(env_temp)))
  na_problem <- na_problem[na_problem$na_sum > missing_day_threshold, , F]
  problematic_years <- paste0(row.names(na_problem), sep = "", collapse=", ")

  if (nrow(na_problem) > 0){

    warning(paste0("Problematic years with missing values are present: ", problematic_years))

  }

  # Previous-year climate data have already been added by
  # build_lagged_env_data().

  # If row_names_subset == TRUE, data is subseted and ordered based on matching
  # row.names. Additionally, number of characters in row.names is checked.
  # There should be at least three characters (assuming years before 100 will
  # never be analysed, there is no such environmental data available)
  if (row_names_subset == TRUE & nchar(row.names(env_data)[1]) >= 3){

    ncol_response <- ncol(response)

    colnames_response <- colnames(response)

    env_data$yearABC <- row.names(env_data)
    response$yearABC <- row.names(response)

    temporal_data <- merge(response, env_data, by = "yearABC")

    response <- data.frame(temporal_data[, c(2:(1 + ncol_response))],
                           row.names = temporal_data$yearABC)
    colnames(response) <- colnames_response

    env_data <- data.frame(temporal_data[, c((1 + ncol_response + 1):
                                               ncol(temporal_data))],
                           row.names = temporal_data$yearABC)
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


  # Subset of years
  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    if (lower_subset > upper_subset){
      stop("Change the order of elements in the subset_years argument! First element should be lower than the second!")
    }

    subset_seq <- seq(lower_subset, upper_subset)

    if (any(!(subset_seq %in% row.names(response)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  " in the response data frame. Change the subset_years argument"))
    }

    if (any(!(subset_seq %in% row.names(env_data)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  " in the env_data data frame. Change the subset_years argument"))
    }

    response <- subset(response, row.names(response) %in% subset_seq)
    env_data <- subset(env_data, row.names(env_data) %in% subset_seq)
  }

  # NA values are not allowed and must be removed from response data.frame
  # exception if cor_na_us accounts for missing values
  if (sum(is.na(response)) > 0 & !(cor_na_use %in% c("complete.obs", "na.or.complete", "pairwise.complete.obs"))){

    prob_year <- row.names(response[is.na(response), , drop = F])

    stop(paste0("NA is not allowed in response data frame. ",
                "Problematic year is ", prob_year))

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

  # this is a list for climate and and holder for saving mm
  list_climate <- list()
  mm <- 1


  # A.1 method = "cor"
  if (fixed_width != 0 & method == "cor") {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be stored in this matrix
    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if(interactive()){

      if (fixed_width != max_window){
        pb <- txtProgressBar(min = 0, max = ceiling((ncol(env_data) - fixed_width - offset_end - offset_start + 1)/(skip_window_position)),
                             style = 3)
      }

    }

    b = 0

    # An iterating loop. In each itteration x is calculated and represents
    # response (dependent) variable. X is a moving average. Window width of
    # a moving window is fixed_width. Next, statistical metric is calculated
    # based on a selected method (cor, lm or brnn). Calculation is stored in
    # temporal matrix.

    for (j in (seq((0 + offset_start -1), (ncol(env_data) - max((fixed_width + offset_end), offset_end)), by = skip_window_position))) {

      b = b + 1

      x <- aggregate_daily_window(
        env_data[1:nrow(env_data), (1 + j):(j + fixed_width), drop = FALSE]
      )

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          tmp_model <- lm(x ~ seq(1:length(x)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x - tmp_pred

          x <- data.frame(x = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)

      }

      x_list <- x
      colnames(x_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x_list) <- row.names(env_data)
      list_climate[[mm]] <- x_list
      mm = mm + 1

      if (boot == FALSE){

        temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method, use = cor_na_use)
        temporal_lower <- NA
        temporal_upper <- NA

      } else if (boot == TRUE){

        temp_df_boot <- cbind(response[, 1], x[, 1])
        calc <- boot(temp_df_boot, boot_f_cor, cor.type = cor_method, R = boot_n)

        temporal_correlation <- colMeans(calc$t)[1]

        ci_int <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type), silent = TRUE)

        if (class(ci_int)[[1]] == "try-error"){

          temporal_lower <- NA
          temporal_upper <- NA

        } else {

          if (boot_ci_type == "norm"){

            temporal_lower <- ci_int$norm[2]
            temporal_upper <- ci_int$norm[3]

          } else if (boot_ci_type == "perc"){

            temporal_lower <- ci_int$perc[4]
            temporal_upper <- ci_int$perc[5]

          } else if (boot_ci_type == "stud") {

            temporal_lower <- ci_int$student[4]
            temporal_upper <- ci_int$student[5]

          } else if (boot_ci_type == "basic") {

            temporal_lower <- ci_int$basic[4]
            temporal_upper <- ci_int$basic[5]

          } else if (boot_ci_type == "bca") {

            temporal_lower <- ci_int$bca[4]
            temporal_upper <- ci_int$bca[5]

          } else {

            stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

          }
        }
      } else {
        print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }



      # Each calculation is printed. Reason: usually it takes several minutes
      # to go through all loops and therefore, users might think that R is
      # not responding. But if each calculation is printed, user could be
      # confident, that R is responding.
      if (reference_window == 'start'){
        temporal_matrix[1, j + 1] <- temporal_correlation
        temporal_matrix_lower[1, j + 1] <- temporal_lower
        temporal_matrix_upper[1, j + 1] <- temporal_upper
      } else if (reference_window == 'end'){
        temporal_matrix[1, j + fixed_width] <- temporal_correlation
        temporal_matrix_lower[1, j + fixed_width] <- temporal_lower
        temporal_matrix_upper[1, j + fixed_width] <- temporal_upper
      } else if (reference_window == 'middle'){
        temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_correlation
        temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_lower
        temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_upper
      }

      if(interactive()){
        if (fixed_width != max_window){setTxtProgressBar(pb, b)}
      }
    }


    if(interactive()){
      if (fixed_width != max_window){close(pb)}
    }

    # temporal_matrix is given rownames and colnames. Rownames represent a
    # window width used fot calculations. Colnames represent the position of
    # moving window in a original env_data data frame.
    row.names(temporal_matrix) <- fixed_width
    row.names(temporal_matrix_lower) <- fixed_width
    row.names(temporal_matrix_upper) <- fixed_width

    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # A.2 method == "lm"
  # For a description see A.1
  if (fixed_width != 0 & method == "lm") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (fixed_width != max_window){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = ceiling((ncol(env_data) - fixed_width - offset_end - offset_start + 1)/(skip_window_position)),
                             style = 3)}
    }

    b = 0

    for (j in (seq((0 + offset_start -1), (ncol(env_data) - max((fixed_width + offset_end), offset_end)), by = skip_window_position))) {



      b = b + 1

      x <- aggregate_daily_window(
        env_data[1:nrow(env_data), (1 + j):(j + fixed_width), drop = FALSE]
      )

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          tmp_model <- lm(x ~ seq(1:length(x)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x - tmp_pred

          x <- data.frame(x = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)

      }

      x_list <- x
      colnames(x_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x_list) <- row.names(env_data)
      list_climate[[mm]] <- x_list
      mm = mm + 1


      if (boot == FALSE){

        temporal_df <- data.frame(cbind(x, response))
        temporal_model <- lm(x ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        temporal_r_squared <- temporal_summary$r.squared
        temporal_adj_r_squared <- temporal_summary$adj.r.squared

        temporal_r_squared_lower<- NA
        temporal_r_squared_upper<- NA
        temporal_adj_r_squared_lower <- NA
        temporal_adj_r_squared_upper <- NA

      } else if (boot == TRUE){

        temporal_df <- data.frame(cbind(x, response))
        calc <- boot(data = temporal_df, statistic = boot_f_lm, R = boot_n, lm.formula = "x ~ .")

        temporal_r_squared <- colMeans(calc$t)[1]
        temporal_adj_r_squared <- colMeans(calc$t)[2]

        ci_int_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1), silent = TRUE)
        ci_int_adj_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2), silent = TRUE)


        if (class(ci_int_r_squared)[[1]] == "try-error"){

          temporal_r_squared_lower<- NA
          temporal_r_squared_upper<- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else {

          if (boot_ci_type == "norm"){

            temporal_r_squared_lower <- ci_int_r_squared$norm[2]
            temporal_r_squared_upper <- ci_int_r_squared$norm[3]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$norm[2]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$norm[3]

          } else if (boot_ci_type == "perc"){

            temporal_r_squared_lower <- ci_int_r_squared$perc[4]
            temporal_r_squared_upper <- ci_int_r_squared$perc[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$perc[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$perc[5]

          } else if (boot_ci_type == "stud") {

            temporal_r_squared_lower <- ci_int_r_squared$studen[4]
            temporal_r_squared_upper <- ci_int_r_squared$studen[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$studen[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$studen[5]

          } else if (boot_ci_type == "basic") {

            temporal_r_squared_lower <- ci_int_r_squared$basic[4]
            temporal_r_squared_upper <- ci_int_r_squared$basic[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$basic[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$basic[5]

          } else if (boot_ci_type == "bca") {

            temporal_r_squared_lower <- ci_int_r_squared$bca[4]
            temporal_r_squared_upper <- ci_int_r_squared$bca[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$bca[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$bca[5]

          } else {

            stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

          }

        }

      } else {
        stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }

      if (metric == "r.squared"){

        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_r_squared_upper

        } else if (reference_window == 'end') {

          temporal_matrix[1, j + fixed_width] <- temporal_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_upper

        }
      }

      if (metric == "adj.r.squared"){
        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_adj_r_squared_upper

        } else if (reference_window == 'end'){

          temporal_matrix[1, j + fixed_width] <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_adj_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_upper
        }
      }

      if(interactive()){
        if (fixed_width != max_window){setTxtProgressBar(pb, b)}
      }
    }

    if(interactive()){
      if (fixed_width != max_window){close(pb)}
    }


    row.names(temporal_matrix) <- fixed_width
    row.names(temporal_matrix_lower) <- fixed_width
    row.names(temporal_matrix_upper) <- fixed_width

    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # A.3 method == "brnn"
  # For a description see A.1
  if (fixed_width != 0 & method == "brnn") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (fixed_width != max_window){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = ceiling((ncol(env_data) - fixed_width - offset_end - offset_start + 1)/(skip_window_position)),
                             style = 3)}
    }
    b = 0

    for (j in (seq((0 + offset_start -1), (ncol(env_data) - max((fixed_width + offset_end), offset_end)), by = skip_window_position))) {

      b = b + 1

      x <- aggregate_daily_window(
        env_data[1:nrow(env_data), (1 + j):(j + fixed_width), drop = FALSE]
      )

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          tmp_model <- lm(x ~ seq(1:length(x)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x - tmp_pred

          x <- data.frame(x = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)

      }

      x_list <- x
      colnames(x_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x_list) <- row.names(env_data)
      list_climate[[mm]] <- x_list
      mm = mm + 1

      if (boot == FALSE){

        temporal_df <- data.frame(cbind(x, response))

        capture.output(temporal_model <- try(brnn(x ~ ., data = temporal_df,
                                                  neurons = neurons, tol = 1e-6),
                                             silent = TRUE))

        temporal_predictions <- try(predict.brnn(temporal_model, temporal_df),
                                    silent = TRUE)


        if (class(temporal_model)[[1]] != "try-error"){

          temporal_r_squared <- 1 - (sum((x[, 1] - temporal_predictions) ^ 2) /
                                       sum((x[, 1] - mean(x[, 1])) ^ 2))
          temporal_adj_r_squared <- 1 - ((1 - temporal_r_squared) *
                                           ((nrow(x) - 1)) /
                                           (nrow(x) -
                                              ncol(as.data.frame(temporal_df))
                                            -  1 + 1))


        }

        temporal_r_squared_lower <- NA
        temporal_r_squared_upper <- NA

        temporal_adj_r_squared_lower <- NA
        temporal_adj_r_squared_upper <- NA

      } else if (boot == TRUE){


        temporal_df <- data.frame(cbind(x, response))
        calc <- boot(data = temporal_df, statistic = boot_f_brnn, R = boot_n, brnn.formula = "x ~ .", neurons = neurons)

        temporal_r_squared <- colMeans(calc$t)[1]
        temporal_adj_r_squared <- colMeans(calc$t)[2]

        ci_int_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1), silent = TRUE)
        ci_int_adj_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2), silent = TRUE)

        if (class(ci_int_r_squared)[[1]] == "try-error"){

          temporal_r_squared_lower<- NA
          temporal_r_squared_upper<- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else {

          if (boot_ci_type == "norm"){

            temporal_r_squared_lower <- ci_int_r_squared$norm[2]
            temporal_r_squared_upper <- ci_int_r_squared$norm[3]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$norm[2]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$norm[3]

          } else if (boot_ci_type == "perc"){

            temporal_r_squared_lower <- ci_int_r_squared$perc[4]
            temporal_r_squared_upper <- ci_int_r_squared$perc[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$perc[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$perc[5]

          } else if (boot_ci_type == "stud") {

            temporal_r_squared_lower <- ci_int_r_squared$studen[4]
            temporal_r_squared_upper <- ci_int_r_squared$studen[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$studen[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$studen[5]

          } else if (boot_ci_type == "basic") {

            temporal_r_squared_lower <- ci_int_r_squared$basic[4]
            temporal_r_squared_upper <- ci_int_r_squared$basic[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$basic[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$basic[5]

          } else if (boot_ci_type == "bca") {

            temporal_r_squared_lower <- ci_int_r_squared$bca[4]
            temporal_r_squared_upper <- ci_int_r_squared$bca[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$bca[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$bca[5]

          } else {

            stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

          }

        }

      } else {

        stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))

      }



      if (metric == "r.squared"){

        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_r_squared_upper

        } else if (reference_window == 'end') {

          temporal_matrix[1, j + fixed_width] <- temporal_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_upper

        }
      }

      if (metric == "adj.r.squared"){
        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_adj_r_squared_upper

        } else if (reference_window == 'end'){

          temporal_matrix[1, j + fixed_width] <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_adj_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_upper
        }
      }
      if(interactive()){
        if (fixed_width != max_window){setTxtProgressBar(pb, b)}
      }
    }

    if(interactive()){
      if (fixed_width != max_window){close(pb)}
    }

    row.names(temporal_matrix) <- fixed_width
    row.names(temporal_matrix_lower) <- fixed_width
    row.names(temporal_matrix_upper) <- fixed_width

    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # B fixed_width == 0, in this case, lower_limit and upper_limit arguments
  # will be used to define window width of a moving window.
  # B.1 method == "cor"

  if (fixed_width == 0 & method == "cor") {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be stored in this matrix

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    # An iterating double loop: 1 outer loop) iterating from lower_limit :
    # upper_limit defines windo.width used for a moving window. 2) inner loop
    # defines the starting position of a moving window.
    # In each itteration, x is calculated and represents a response (dependent)
    # variable. x is a moving average, based on rowMeans/apply function.
    # Next, statistical metric is calculated based on a selected method (cor,
    # lm or brnn). Calculation is stored in temporal matrix in a proper place.
    # The position of stored calculation is informative later used for
    # indiciating optimal values.

    if (upper_limit != lower_limit){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = ceiling((upper_limit - lower_limit)/(skip_window_length*skip_window_position)), style = 3)

      }
    }

    b = 0



    for (K in seq(lower_limit, upper_limit, by = skip_window_length)) {

      b = b + 1

      for (j in seq((0 + offset_start -1), (ncol(env_data) - max((K + offset_end), offset_end)), by = skip_window_position)) {

        x <- aggregate_daily_window(
          env_data[1:nrow(env_data), (1 + j):(j + K), drop = FALSE]
        )

        if (!is.null(dc_method)){

          if (dc_method == "SLD"){

            tmp_model <- lm(x ~ seq(1:length(x)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x - tmp_pred

            x <- data.frame(x = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x <- matrix(x, nrow = nrow(env_data), ncol = 1)

        }

        # }
        x_list <- x
        colnames(x_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x_list) <- row.names(env_data)
        list_climate[[mm]] <- x_list
        mm = mm + 1

        if (boot == FALSE){

          temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method, use = cor_na_use)
          temporal_lower <- NA
          temporal_upper <- NA

        } else if (boot == TRUE){

          temp_df_boot <- cbind(response[, 1], x[, 1])
          calc <- boot(temp_df_boot, boot_f_cor, cor.type = cor_method, R = boot_n)

          temporal_correlation <- colMeans(calc$t)[1]

          ci_int <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type), silent = TRUE)

          if (class(ci_int)[[1]] == "try-error"){

            temporal_lower <- NA
            temporal_upper <- NA

          } else {

            if (boot_ci_type == "norm"){

              temporal_lower <- ci_int$norm[2]
              temporal_upper <- ci_int$norm[3]

            } else if (boot_ci_type == "perc"){

              temporal_lower <- ci_int$perc[4]
              temporal_upper <- ci_int$perc[5]

            } else if (boot_ci_type == "stud") {

              temporal_lower <- ci_int$student[4]
              temporal_upper <- ci_int$student[5]

            } else if (boot_ci_type == "basic") {

              temporal_lower <- ci_int$basic[4]
              temporal_upper <- ci_int$basic[5]

            } else if (boot_ci_type == "bca") {

              temporal_lower <- ci_int$bca[4]
              temporal_upper <- ci_int$bca[5]

            } else {

              stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

            }
          }
        } else {
          print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
        }

        if (reference_window == 'start'){
          temporal_matrix[(K - lower_limit) + 1, j + 1] <- temporal_correlation
          temporal_matrix_lower[(K - lower_limit) + 1, j + 1] <- temporal_lower
          temporal_matrix_upper[(K - lower_limit) + 1, j + 1] <- temporal_upper
        } else if (reference_window == 'end'){
          temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_correlation
          temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_lower
          temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_upper
        } else if (reference_window == 'middle'){
          temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_correlation
          temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_lower
          temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_upper
        }


      }

      if(interactive()){
        if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
      }
    }

    if(interactive()){
      if (upper_limit != lower_limit){close(pb)}
    }

    # temporal_matrix is given rownames and colnames. Rownames represent a
    # window width used fot calculations. Colnames represent the position of
    # moving window in a original env_data data frame.
    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit,
                                       by = 1))
    row.names(temporal_matrix) <- temporal_rownames
    row.names(temporal_matrix_lower) <- temporal_rownames
    row.names(temporal_matrix_upper) <- temporal_rownames


    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # B.2 method == "lm"
  # For a description see B.1
  if (fixed_width == 0 & method == "lm") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (upper_limit != lower_limit){
      if(interactive()){
        pb <- txtProgressBar(min = 0, max = ceiling((upper_limit - lower_limit)/(skip_window_length*skip_window_position)), style = 3)
      }
    }
    b = 0

    for (K in seq(lower_limit, upper_limit, by = skip_window_length)) {

      b = b + 1

      for (j in seq((0 + offset_start -1), (ncol(env_data) - max((K + offset_end), offset_end)), by = skip_window_position)) {


        x <- aggregate_daily_window(
          env_data[1:nrow(env_data), (1 + j):(j + K), drop = FALSE]
        )

        if (!is.null(dc_method)){

          if (dc_method == "SLD"){

            tmp_model <- lm(x ~ seq(1:length(x)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x - tmp_pred

            x <- data.frame(x = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x <- matrix(x, nrow = nrow(env_data), ncol = 1)

        }

        x_list <- x
        colnames(x_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x_list) <- row.names(env_data)
        list_climate[[mm]] <- x_list
        mm = mm + 1

        if (boot == FALSE){

          temporal_df <- data.frame(cbind(x, response))
          temporal_model <- lm(x ~ ., data = temporal_df)
          temporal_summary <- summary(temporal_model)
          temporal_r_squared <- temporal_summary$r.squared
          temporal_adj_r_squared <- temporal_summary$adj.r.squared

          temporal_r_squared_lower<- NA
          temporal_r_squared_upper<- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else if (boot == TRUE){

          temporal_df <- data.frame(cbind(x, response))
          calc <- boot(data = temporal_df, statistic = boot_f_lm, R = boot_n, lm.formula = "x ~ .")

          temporal_r_squared <- colMeans(calc$t)[1]
          temporal_adj_r_squared <- colMeans(calc$t)[2]

          ci_int_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1), silent = TRUE)
          ci_int_adj_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2), silent = TRUE)


          if (class(ci_int_r_squared)[[1]] == "try-error"){

            temporal_r_squared_lower<- NA
            temporal_r_squared_upper<- NA
            temporal_adj_r_squared_lower <- NA
            temporal_adj_r_squared_upper <- NA

          } else {

            if (boot_ci_type == "norm"){

              temporal_r_squared_lower <- ci_int_r_squared$norm[2]
              temporal_r_squared_upper <- ci_int_r_squared$norm[3]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$norm[2]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$norm[3]

            } else if (boot_ci_type == "perc"){

              temporal_r_squared_lower <- ci_int_r_squared$perc[4]
              temporal_r_squared_upper <- ci_int_r_squared$perc[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$perc[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$perc[5]

            } else if (boot_ci_type == "stud") {

              temporal_r_squared_lower <- ci_int_r_squared$studen[4]
              temporal_r_squared_upper <- ci_int_r_squared$studen[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$studen[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$studen[5]

            } else if (boot_ci_type == "basic") {

              temporal_r_squared_lower <- ci_int_r_squared$basic[4]
              temporal_r_squared_upper <- ci_int_r_squared$basic[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$basic[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$basic[5]

            } else if (boot_ci_type == "bca") {

              temporal_r_squared_lower <- ci_int_r_squared$bca[4]
              temporal_r_squared_upper <- ci_int_r_squared$bca[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$bca[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$bca[5]

            } else {

              stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

            }

          }

        } else {
          stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))
        }

        if (metric == "r.squared"){

          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_upper
          }

        }

        if (metric == "adj.r.squared"){


          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_upper
          }

        }
      }
      if(interactive()){
        if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
      }
    }

    if(interactive()){
      if (upper_limit != lower_limit){close(pb)}
    }

    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit, by = 1))
    row.names(temporal_matrix) <- temporal_rownames
    row.names(temporal_matrix_lower) <- temporal_rownames
    row.names(temporal_matrix_upper) <- temporal_rownames

    temporal_colnames <- as.vector(seq(from = 1, to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames

  }

  # B.3 method == "brnn"
  # For a description see B.1
  if (fixed_width == 0 & method == "brnn") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }
    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (upper_limit != lower_limit){
      if(interactive()){
        pb <- txtProgressBar(min = 0, max = ceiling((upper_limit - lower_limit)/(skip_window_length*skip_window_position)), style = 3)
      }
    }

    b = 0

    for (K in seq(lower_limit, upper_limit, by = skip_window_length)) {

      b = b + 1

      for (j in seq((0 + offset_start -1), (ncol(env_data) - max((K + offset_end), offset_end)), by = skip_window_position)) {

        x <- aggregate_daily_window(
          env_data[1:nrow(env_data), (1 + j):(j + K), drop = FALSE]
        )

        if (!is.null(dc_method)){

          if (dc_method == "SLD"){

            tmp_model <- lm(x ~ seq(1:length(x)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x - tmp_pred

            x <- data.frame(x = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x <- matrix(x, nrow = nrow(env_data), ncol = 1)

        }

        x_list <- x
        colnames(x_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x_list) <- row.names(env_data)
        list_climate[[mm]] <- x_list
        mm = mm + 1

        if (boot == FALSE){

          temporal_df <- data.frame(cbind(x, response))

          capture.output(temporal_model <- try(brnn(x ~ ., data = temporal_df,
                                                    neurons = neurons, tol = 1e-6),
                                               silent = TRUE))

          temporal_predictions <- try(predict.brnn(temporal_model, temporal_df),
                                      silent = TRUE)


          if (class(temporal_model)[[1]] != "try-error"){

            temporal_r_squared <- 1 - (sum((x[, 1] - temporal_predictions) ^ 2) /
                                         sum((x[, 1] - mean(x[, 1])) ^ 2))
            temporal_adj_r_squared <- 1 - ((1 - temporal_r_squared) *
                                             ((nrow(x) - 1)) /
                                             (nrow(x) -
                                                ncol(as.data.frame(temporal_df))
                                              -  1 + 1))


          } else {
            temporal_r_squared <- NA
            temporal_adj_r_squared <- NA
          }

          temporal_r_squared_lower <- NA
          temporal_r_squared_upper <- NA

          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA


        } else if (boot == TRUE){

          temporal_df <- data.frame(cbind(x, response))
          calc <- boot(data = temporal_df, statistic = boot_f_brnn, R = boot_n, brnn.formula = "x ~ .", neurons = neurons)

          temporal_r_squared <- colMeans(calc$t)[1]
          temporal_adj_r_squared <- colMeans(calc$t)[2]

          ci_int_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1), silent = TRUE)
          ci_int_adj_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2), silent = TRUE)

          if (class(ci_int_r_squared)[[1]] == "try-error"){

            temporal_r_squared_lower<- NA
            temporal_r_squared_upper<- NA
            temporal_adj_r_squared_lower <- NA
            temporal_adj_r_squared_upper <- NA

          } else {

            if (boot_ci_type == "norm"){

              temporal_r_squared_lower <- ci_int_r_squared$norm[2]
              temporal_r_squared_upper <- ci_int_r_squared$norm[3]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$norm[2]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$norm[3]

            } else if (boot_ci_type == "perc"){

              temporal_r_squared_lower <- ci_int_r_squared$perc[4]
              temporal_r_squared_upper <- ci_int_r_squared$perc[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$perc[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$perc[5]

            } else if (boot_ci_type == "stud") {

              temporal_r_squared_lower <- ci_int_r_squared$studen[4]
              temporal_r_squared_upper <- ci_int_r_squared$studen[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$studen[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$studen[5]

            } else if (boot_ci_type == "basic") {

              temporal_r_squared_lower <- ci_int_r_squared$basic[4]
              temporal_r_squared_upper <- ci_int_r_squared$basic[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$basic[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$basic[5]

            } else if (boot_ci_type == "bca") {

              temporal_r_squared_lower <- ci_int_r_squared$bca[4]
              temporal_r_squared_upper <- ci_int_r_squared$bca[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$bca[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$bca[5]

            } else {

              stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

            }

          }

        } else {

          stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))

        }

        if (metric == "r.squared"){

          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_upper
          }

        }

        if (metric == "adj.r.squared"){


          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_upper
          }

        }

      }
      if(interactive()){
        if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
      }
    }

    if(interactive()){
      if (upper_limit != lower_limit){close(pb)}
    }

    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit, by = 1))
    row.names(temporal_matrix) <- temporal_rownames
    row.names(temporal_matrix_lower) <- temporal_rownames
    row.names(temporal_matrix_upper) <- temporal_rownames

    temporal_colnames <- as.vector(seq(from = 1, to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames

  }

  # PART 3: smoothing function, if brnn method is used. It turnes out, that
  # brnn sometimes (1 - 3 % of calculations) fails to construct a realistic
  # result. In those cases, much lower r.squared (or adj.r.squared) are
  # calculated. smooth_matrix function removes unrealistic calculations and
  # replace them with mean of values in a window 3 x 3. Maximum value
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

    # 1 Method is correlation
    if (method == "cor") {

      temporal_matrix[abs(temporal_matrix) < abs(critical_threshold_cor)] <- NA

      # 2 lm and brnn method
    } else if (method == "lm" | method == "brnn") {
      temporal_matrix[abs(temporal_matrix) < abs(critical_threshold_cor2)] <- NA
    }
  }




  if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){

    warning("All calculations are insignificant! Please change the alpha argument.")

    final_list <- list(calculations = temporal_matrix,
                       method = method,
                       metric = cor_method,
                       analysed_period = NA,
                       optimized_return = NA,
                       optimized_return_all = NA,
                       transfer_function = NA, temporal_stability = NA,
                       cross_validation = NA,
                       plot_heatmap = NA,
                       plot_extreme = NA,
                       type = "daily",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate = NA,
                       previous_year = previous_year,
                       number_previous_years = number_previous_years)

    class(final_list) <- 'dmrs'

    return(final_list)


  } else {

    ########################################################################
    # PART 4: Final list is being created and returned as a function output#
    ########################################################################

    # The first three elements of the final list are already created: calculated
    # values, method and metric used.
    # Here we create the fourth element: the optimal sequence of days that
    # returns the best selected statistical metric. We name it optimal_return.

    # In case of negative correlations, different strategy is applied.
    # For more detailed description see plot_extreme()

    if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){
      stop("All calculations are insignificant! Change the alpha argument!")
    }

    overall_max <- max(temporal_matrix, na.rm = TRUE)
    overall_min <- min(temporal_matrix, na.rm = TRUE)

    # absolute vales of overall_maximum and overall_minimum are compared and
    # one of the following two if functions is used
    # There are unimportant warnings produced:
    # no non-missing arguments to max; returning -Inf

    if ((abs(overall_max) >= abs(overall_min)) == TRUE) {

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

    # The fourth return element is being created: rowMeans/ apply of optimal sequence:
    # So, here we consider more options, based on the reference_winow
    # option 1: reference window = "start"
    if (reference_window == 'start'){

      dataf <- data.frame(aggregate_daily_window(
        env_data[, as.numeric(plot_column):
                   (as.numeric(plot_column) + as.numeric(row_index) - 1),
                 drop = FALSE]
      ))

      # if detrending was applied, should also be applied here
      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          dataf <- as.numeric(dataf[, 1])
          tmp_model <- lm(dataf ~ seq(1:length(dataf)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- dataf - tmp_pred

          dataf <- data.frame(tmp_res / sd(tmp_res, na.rm = TRUE))

        }

      }

      dataf_full <- cbind(response, dataf)
      colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      dataf_original <- data.frame(aggregate_daily_window(
        env_data_original[, as.numeric(plot_column):
                            (as.numeric(plot_column) + as.numeric(row_index) - 1),
                          drop = FALSE]
      ))

      dataf_full_original <- dataf_original

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          dataf_full_original <- as.numeric(dataf_full_original[, 1])
          tmp_model <- lm(dataf_full_original ~ seq(1:length(dataf_full_original)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- dataf_full_original - tmp_pred

          dataf_full_original <- data.frame(tmp_res / sd(tmp_res, na.rm = TRUE))

        }

      }

      colnames(dataf_full_original) <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      # Additional check: (we should get the same metric as before in the loop)
      if (method == "lm" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$r.squared
      }

      if (method == "lm" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$adj.r.squared
      }

      if (method == "brnn" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        capture.output(temporal_model <- brnn(Optimized.rowNames ~ ., data = temporal_df,
                                              neurons = neurons, tol = 1e-6))
        temporal_predictions <- try(predict.brnn(temporal_model,
                                                 temporal_df), silent = TRUE)
        optimized_result <- 1 - (sum((temporal_df[, 1] -
                                        temporal_predictions) ^ 2) /
                                   sum((temporal_df[, 1] -
                                          mean(temporal_df[, 1])) ^ 2))
      }

      if (method == "brnn" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        capture.output(temporal_model <- brnn(Optimized.rowNames ~ .,
                                              data = temporal_df, neurons = neurons, tol = 1e-6))
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
        optimized_result <- cor(dataf, response, method = cor_method, use = cor_na_use)
      }

      # Just give a nicer colname
      colnames(dataf) <- "Optimized return"

    }


    # Option 2: reference window = "end"
    if (reference_window == 'end'){

      dataf <- data.frame(aggregate_daily_window(
        env_data[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                   as.numeric(plot_column),
                 drop = FALSE]
      ))

      # if detrending was applied, should also be applied here
      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          dataf <- as.numeric(dataf[, 1])
          tmp_model <- lm(dataf ~ seq(1:length(dataf)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- dataf - tmp_pred

          dataf <- data.frame(tmp_res / sd(tmp_res, na.rm = TRUE))

        }

      }

      dataf_full <- cbind(response, dataf)
      colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      dataf_original <- data.frame(aggregate_daily_window(
        env_data_original[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                            as.numeric(plot_column),
                          drop = FALSE]
      ))

      dataf_full_original <- dataf_original

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          dataf_full_original <- as.numeric(dataf_full_original[, 1])
          tmp_model <- lm(dataf_full_original ~ seq(1:length(dataf_full_original)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- dataf_full_original - tmp_pred

          dataf_full_original <- data.frame(tmp_res / sd(tmp_res, na.rm = TRUE))

        }

      }

      colnames(dataf_full_original) <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      # Additional check: (we should get the same metric as before in the loop)
      if (method == "lm" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$r.squared
      }

      if (method == "lm" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$adj.r.squared
      }

      if (method == "brnn" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        capture.output(temporal_model <- brnn(Optimized.rowNames ~ ., data = temporal_df,
                                              neurons = neurons, tol = 1e-6))
        temporal_predictions <- try(predict.brnn(temporal_model,
                                                 temporal_df), silent = TRUE)
        optimized_result <- 1 - (sum((temporal_df[, 1] -
                                        temporal_predictions) ^ 2) /
                                   sum((temporal_df[, 1] -
                                          mean(temporal_df[, 1])) ^ 2))
      }

      if (method == "brnn" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        capture.output(temporal_model <- brnn(Optimized.rowNames ~ .,
                                              data = temporal_df, neurons = neurons, tol = 1e-6))
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
        optimized_result <- cor(dataf, response, method = cor_method, use = cor_na_use)
      }

      # Just give a nicer colname
      colnames(dataf) <- "Optimized return"

    }


    # Option 3: reference window = "middle"
    if (reference_window == 'middle'){

      if (as.numeric(row_index) %% 2 == 0){
        adjustment_1 <- 0
        adjustment_2 <- 1
      } else {
        adjustment_1 <- 1
        adjustment_2 <- 2
      }

      dataf <- data.frame(aggregate_daily_window(
        env_data[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                   (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                 drop = FALSE]
      ))

      # if detrending was applied, should also be applied here
      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          dataf <- as.numeric(dataf[, 1])
          tmp_model <- lm(dataf ~ seq(1:length(dataf)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- dataf - tmp_pred

          dataf <- data.frame(tmp_res / sd(tmp_res, na.rm = TRUE))

        }

      }

      dataf_full <- cbind(response, dataf)
      colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      dataf_original <- data.frame(aggregate_daily_window(
        env_data_original[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                            (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                          drop = FALSE]
      ))

      dataf_full_original <- dataf_original

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          dataf_full_original <- as.numeric(dataf_full_original[, 1])
          tmp_model <- lm(dataf_full_original ~ seq(1:length(dataf_full_original)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- dataf_full_original - tmp_pred

          dataf_full_original <- data.frame(tmp_res / sd(tmp_res, na.rm = TRUE))

        }

      }

      colnames(dataf_full_original) <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      # Additional check: (we should get the same metric as before in the loop)
      if (method == "lm" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$r.squared
      }

      if (method == "lm" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$adj.r.squared
      }

      if (method == "brnn" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        capture.output(temporal_model <- brnn(Optimized.rowNames ~ ., data = temporal_df,
                                              neurons = neurons, tol = 1e-6))
        temporal_predictions <- try(predict.brnn(temporal_model,
                                                 temporal_df), silent = TRUE)
        optimized_result <- 1 - (sum((temporal_df[, 1] -
                                        temporal_predictions) ^ 2) /
                                   sum((temporal_df[, 1] -
                                          mean(temporal_df[, 1])) ^ 2))
      }

      if (method == "brnn" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        capture.output(temporal_model <- brnn(Optimized.rowNames ~ .,
                                              data = temporal_df, neurons = neurons, tol = 1e-6))
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
        optimized_result <- cor(dataf, response, method = cor_method, use = cor_na_use)
      }

      # Just give a nicer colname
      colnames(dataf) <- "Optimized return"

    }
    ##############################################################################
    # If detrending was used, it also needs to be applied on optimized return

    # Element 5
    # Here we create the fifth element of the final list: Analysed period in the
    # form of min(year) - max(year), e.g. 1950 - 2015
    min_env_data <- min(as.numeric(row.names(env_data)))
    min_response <- min(as.numeric(row.names(response)))

    max_env_data <- max(as.numeric(row.names(env_data)))
    max_response <- max(as.numeric(row.names(response)))

    min_together <- min(min_env_data, min_response)
    max_together <- min(max_env_data, max_response)


    analysed_period <- paste(as.character(min_together),
                             as.character(max_together),
                             sep = " - ")
    if (nchar(analysed_period) < 9) {
      analysed_period <- NA
    }

    # Here, the transfer function is being created
    transfer_data = data.frame(proxy = response[,1], optimized_return =dataf[,1])
    lm_model = lm(optimized_return ~ proxy, data = transfer_data)
    capture.output(brnn_model <- try(brnn(optimized_return ~ proxy, data = transfer_data, neurons = neurons), silent = TRUE))
    full_range = data.frame(proxy = seq(from = min(response[,1], na.rm = TRUE), to = max(response[,1], na.rm = TRUE), length.out = 100))

    if (method == "lm" | method == "cor"){
      full_range$transfer_f = predict(lm_model, full_range)
    }

    if (method == "brnn"){
      full_range$transfer_f = predict.brnn(brnn_model, full_range)
    }

    # String for titles
    if (method == "cor"){
      title_string <- "Correlation Coefficients"
    } else if (method == "lm"){
      title_string <- "Linear Regression"
    } else if (method == "brnn"){
      title_string <- "ANN With Bayesian Regularization"
    } else (print("The selection of method is not correct"))

    # The definition of theme
    journal_theme <- theme_bw() +
      theme(axis.text = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 18), text = element_text(size = 18),
            plot.title = element_text(size = 16,  face = "bold"))

    p1 <- ggplot(transfer_data, aes(proxy, optimized_return)) +
      geom_point() +
      geom_line(aes(proxy, transfer_f), full_range) +
      journal_theme +
      ggtitle(paste("Analysed Period:", analysed_period, "\nMethod:", title_string))

    # If there is more than one independent variable in the model,
    # transfer function is not given, since we should return a 3d model
    if (ncol(response) > 1){
      p1 <- "No transfer function is created for two or more response variables"
    }


    #######################################################################
    ############## The temporal stability of optimized_return #############
    #######################################################################

    dataset = data.frame(optimized_return = dataf[,1], proxy = response)

    empty_list = list()
    empty_list_period = list()
    empty_list_significance = list()

    temporal_stability <- data.frame()

    # 1. Progressive stability check
    if (temporal_stability_check == "progressive"){
      foldi <- seq(1:k)
      #foldi <- paste("fold_", foldi)
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      for (m in 1:k){
        #Segement your data by fold using the which() function
        trainIndexes <- which(folds <= m, arr.ind = TRUE)
        dataset_temp <- dataset[trainIndexes, ]
        MAKS <- max(as.numeric(row.names(dataset_temp)))
        MIN <- min(as.numeric(row.names(dataset_temp)))
        empty_list_period[[m]] <- paste(MIN, "-", MAKS)

        if (method == "cor"){
          calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method, use = cor_na_use)
          sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exact=F)$p.value
          empty_list_significance[[m]] <- sig
          empty_list[[m]] <- calculation
          colname = "correlation"
        } else if (method == "lm" & metric == "r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          colname = "r.squared"
          empty_list[[m]] <- summary(MLR)$r.squared
          empty_list_significance[[m]] <- NA
        } else if (method == "lm" & metric == "adj.r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          empty_list[[m]] <- summary(MLR)$adj.r.squared
          empty_list_significance[[m]] <- NA
          colname = "adj.r.squared"
        } else if (method == "brnn" & metric == "r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
            predictions <- predict(BRNN, dataset_temp, neurons = neurons)
            r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                                sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))
            empty_list[[m]] <- r_squared
            empty_list_significance[[m]] <- NA
            colname = "r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "r.squared"
            empty_list_significance[[m]] <- NA
          }
        } else if (method == "brnn" & metric == "adj.r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
            predictions <- predict(BRNN, dataset_temp, neurons = neurons)
            r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                                sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))

            adj_r_squared <- 1 - ((1 - r_squared) * ((nrow(dataset_temp) - 1)) /
                                    (nrow(dataset_temp) - ncol(as.data.frame(response[, 1])) -  1))
            empty_list[[m]] <- adj_r_squared
            empty_list_significance[[m]] <- NA
            colname = "adj.r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "adj.r.squared"
          }
        }
      }
      m1 <- do.call(rbind, empty_list)
      m2 <- do.call(rbind, empty_list_period)
      m3 <- do.call(rbind, empty_list_significance)

      temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(m3, 4), nsmall = 3)))
      colnames(temporal_stability) <-c("Period", colname, "p value")
      temporal_stability
    }

    # 2. Sequential stability check
    if (temporal_stability_check == "sequential"){
      foldi <- seq(1:k)
      #foldi <- paste("fold_", foldi)
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      for (m in 1:k){
        #Segement your data by fold using the which() function
        trainIndexes <- which(folds == m, arr.ind = TRUE)
        dataset_temp <- dataset[trainIndexes, ]

        MAKS <- max(as.numeric(row.names(dataset_temp)))
        MIN <- min(as.numeric(row.names(dataset_temp)))
        empty_list_period[[m]] <- paste(MIN, "-", MAKS)

        if (method == "cor"){
          calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method, use = cor_na_use)
          sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exact=F)$p.value
          empty_list[[m]] <- calculation
          empty_list_significance[[m]] <- sig
          colname = "correlation"

        } else if (method == "lm" & metric == "r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          colname = "r.squared"
          empty_list[[m]] <- summary(MLR)$r.squared
          empty_list_significance[[m]] <- NA
        } else if (method == "lm" & metric == "adj.r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          empty_list[[m]] <- summary(MLR)$adj.r.squared
          empty_list_significance[[m]] <- NA
          colname = "adj.r.squared"
        } else if (method == "brnn" & metric == "r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
            predictions <- predict(BRNN, dataset_temp, neurons = neurons)
            r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                                sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))
            empty_list[[m]] <- r_squared
            empty_list_significance[[m]] <- NA
            colname = "r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "r.squared"
            empty_list_significance[[m]] <- NA
          }
        } else if (method == "brnn" & metric == "adj.r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
            predictions <- predict(BRNN, dataset_temp, neurons = neurons)
            r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                                sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))

            adj_r_squared <- 1 - ((1 - r_squared) * ((nrow(dataset_temp) - 1)) /
                                    (nrow(dataset_temp) - ncol(as.data.frame(response[, 1])) -  1))
            empty_list[[m]] <- adj_r_squared
            empty_list_significance[[m]] <- NA
            colname = "adj.r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "adj.r.squared"
            empty_list_significance[[m]] <- NA
          }
        }
      }
      m1 <- do.call(rbind, empty_list)
      m2 <- do.call(rbind, empty_list_period)
      m3 <- do.call(rbind, empty_list_significance)

      temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(m3, 4), nsmall = 3)))
      colnames(temporal_stability) <-c("Period", colname, "p value")
      temporal_stability

    }


    # 3. running_window
    if (temporal_stability_check == "running_window"){

      k_end <- nrow(dataset)
      place_list <- 1
      empty_list_datasets <- list()
      empty_list_period <- list()

      empty_list <- list()
      empty_list_significance <- list()

      place_list = 1

      if (k_end < k_running_window){
        stop("k_running_window is less than the number of analysed years. Reduce the argument k_running_window")
      }

      for (w in 0:(k_end - k_running_window)){

        dataset_temp <- dataset[(1+w):(k_running_window + w),]
        empty_list_datasets[[place_list]] <- dataset_temp

        MAKS <- max(as.numeric(row.names(dataset_temp)))
        MIN <- min(as.numeric(row.names(dataset_temp)))
        empty_list_period[[place_list]] <- paste(MIN, "-", MAKS)

        place_list <- place_list + 1


      }

      for (m in 1:length(empty_list_datasets)){

        dataset_temp <- empty_list_datasets[[m]]

        if (method == "cor"){
          calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method, use = cor_na_use)
          sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exact=F)$p.value
          empty_list[[m]] <- calculation
          empty_list_significance[[m]] <- sig
          colname = "correlation"

        } else if (method == "lm" & metric == "r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          colname = "r.squared"
          empty_list[[m]] <- summary(MLR)$r.squared
          empty_list_significance[[m]] <- NA
        } else if (method == "lm" & metric == "adj.r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          empty_list[[m]] <- summary(MLR)$adj.r.squared
          empty_list_significance[[m]] <- NA
          colname = "adj.r.squared"
        } else if (method == "brnn" & metric == "r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
            predictions <- predict(BRNN, dataset_temp, neurons = neurons)
            r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                                sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))
            empty_list[[m]] <- r_squared
            empty_list_significance[[m]] <- NA
            colname = "r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "r.squared"
            empty_list_significance[[m]] <- NA
          }
        } else if (method == "brnn" & metric == "adj.r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
            predictions <- predict(BRNN, dataset_temp, neurons = neurons)
            r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                                sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))

            adj_r_squared <- 1 - ((1 - r_squared) * ((nrow(dataset_temp) - 1)) /
                                    (nrow(dataset_temp) - ncol(as.data.frame(response[, 1])) -  1))
            empty_list[[m]] <- adj_r_squared
            empty_list_significance[[m]] <- NA
            colname = "adj.r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "adj.r.squared"
            empty_list_significance[[m]] <- NA
          }
        }
      }
      m1 <- do.call(rbind, empty_list)
      m2 <- do.call(rbind, empty_list_period)
      m3 <- do.call(rbind, empty_list_significance)

      temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(as.numeric(m3), digits = 3), nsmall = 3)))
      colnames(temporal_stability) <-c("Period", colname, "p value")
      temporal_stability
    }

    #########################################################################
    ################## Out of sample estimates ##############################
    #########################################################################

    dataset = data.frame(optimized_return =dataf[,1], proxy = response)

    empty_list = list()
    empty_list_period = list()

    if (cross_validation_type == "blocked"){
      dataset <- dataset
    } else if (cross_validation_type == "randomized"){
      dataset <- dataset[sample(nrow(dataset)), ]
    } else (stop(paste("The cross_validation_type is not selected correctly! It is ", cross_validation_type,
                       ". It should be 'blocked' or 'randomized'!", sep = "")))

    foldi <- seq(1:k)
    #foldi <- paste("fold_", foldi)
    folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)
    bl = 1

    for (m in 1:k){
      #Segement your data by fold using the which() function
      testIndexes <- which(folds == m, arr.ind = TRUE)
      test <- dataset[testIndexes, ]
      train <- dataset[-testIndexes, ]

      MAKS <- max(as.numeric(row.names(train)))
      MIN <- min(as.numeric(row.names(train)))
      empty_list_period[[bl]] <- paste(MIN, "-", MAKS)
      bl <- bl + 1

      MAKS <- max(as.numeric(row.names(test)))
      MIN <- min(as.numeric(row.names(test)))
      empty_list_period[[bl]] <- paste(MIN, "-", MAKS)
      bl <- bl + 1

      if (method == "lm" | method == "cor"){
        MLR <- lm(optimized_return ~ ., data = train)
        train_predicted <- predict(MLR, train)
        test_predicted <- predict(MLR, test)
        train_observed <- train[, 1]
        test_observed <- test[, 1]
        calculations <- calculate_metrics(train_predicted, test_predicted,
                                          train_observed, test_observed, test = test,
                                          formula = optimized_return ~ ., digits = 15)

        empty_list[[m]] <- calculations
      }

      if (method == "brnn"){
        capture.output(BRNN <- try(brnn(optimized_return ~ ., data = train, neurons = neurons), silent = TRUE))
        if (class(BRNN)[[1]] != "try-error"){
          train_predicted <- predict(BRNN, train)
          test_predicted <- predict(BRNN, test)
          train_observed <- train[, 1]
          test_observed <- test[, 1]
          calculations <- calculate_metrics(train_predicted, test_predicted,
                                            train_observed, test_observed, digits = 15,
                                            test = test,
                                            formula = optimized_return ~ .)

          empty_list[[m]] <- calculations

        } else {
          empty_list[[m]] <- NA
        }
      }

    }
    m1 <- do.call(rbind, empty_list)
    # m1 <- m1[, -c(3, 4, 7)]
    m2 <- do.call(rbind, empty_list_period)

    cross_validation <- cbind(Years = m2, m1)
    cross_validation$Period <- c("Calibration", "Validation")
    cross_validation$CV <- rep(1:k, each = 2)
    row.names(cross_validation) <- NULL

    if (cross_validation_type == "blocked"){
      cross_validation <- dplyr::select(cross_validation, CV, Period, Years, cor, RMSE, RRSE, d, RE, CE, DE)
    }
    if (cross_validation_type == "randomized"){
      cross_validation <- dplyr::select(cross_validation, CV, Period, cor, RMSE, RRSE, d, RE, CE, DE)
    }



    ################################################################
    #### Here the final list is being filled with six elements #####
    ################################################################




    # When metohod == "cor", different final_list is created
    if (method == "lm" | method == "brnn") {
      final_list <- list(calculations = temporal_matrix, method = method,
                         metric = metric, analysed_period = analysed_period,
                         optimized_return = dataf_full,
                         optimized_return_all = dataf_full_original,
                         transfer_function = p1, temporal_stability = temporal_stability,
                         cross_validation = cross_validation,
                         aggregated_climate = do.call(cbind, list_climate),
                         previous_year = previous_year,
                         number_previous_years = number_previous_years)
    }

    if (method == "cor"){
      final_list <- list(calculations = temporal_matrix, method = method,
                         metric = cor_method, analysed_period = analysed_period,
                         optimized_return = dataf_full,
                         optimized_return_all = dataf_full_original,
                         transfer_function = p1, temporal_stability = temporal_stability,
                         cross_validation = cross_validation,
                         aggregated_climate = do.call(cbind, list_climate),
                         previous_year = previous_year,
                         number_previous_years = number_previous_years)
    }


    plot_heatmapA <- plot_heatmap(final_list, reference_window = reference_window, type = "daily")
    plot_extremeA <- plot_extreme(final_list, ylimits = ylimits, reference_window = reference_window, type = "daily")

    # Here, for the sake of simplicity, we create final list again
    if (method == "lm" | method == "brnn") {
      final_list <- list(calculations = temporal_matrix, method = method,
                         metric = metric, analysed_period = analysed_period,
                         optimized_return = dataf_full,
                         optimized_return_all = dataf_full_original,
                         transfer_function = p1, temporal_stability = temporal_stability,
                         cross_validation = cross_validation,
                         plot_heatmap = plot_heatmapA,
                         plot_extreme = plot_extremeA,
                         type = "daily",
                         reference_window = reference_window,
                         boot_lower = temporal_matrix_lower,
                         boot_upper = temporal_matrix_upper,
                         aggregated_climate = do.call(cbind, list_climate),
                         previous_year = previous_year,
                         number_previous_years = number_previous_years)
    }

    if (method == "cor"){

      final_list <- list(calculations = temporal_matrix, method = method,
                         metric = cor_method, analysed_period = analysed_period,
                         optimized_return = dataf_full,
                         optimized_return_all = dataf_full_original,
                         transfer_function = p1, temporal_stability = temporal_stability,
                         cross_validation = cross_validation,
                         plot_heatmap = plot_heatmapA,
                         plot_extreme = plot_extremeA,
                         type = "daily",
                         reference_window = reference_window,
                         boot_lower = temporal_matrix_lower,
                         boot_upper = temporal_matrix_upper,
                         aggregated_climate = do.call(cbind, list_climate),
                         previous_year = previous_year,
                         number_previous_years = number_previous_years)
    }

    class(final_list) <- 'dmrs'

    return(final_list)
  }

}
