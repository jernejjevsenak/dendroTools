#' monthly_response_seascorr
#'
#' Function calculates all possible partial correlation coefficients between
#' tree-ring chronology and monthly environmental (usually climate) data.
#' All calculated (partial) correlation coefficients are stored in a matrix.
#' The location of stored correlation in the matrix is indicating a window
#' width (row names) and a location in a matrix of monthly sequences of
#' environmental data (column names).
#'
#' @param response a data frame with tree-ring proxy variable and (optional)
#' years as row names. Row.names should be matched with those from env_data_primary
#' and env_data_control data frame. If not, set the row_names_subset argument to
#' TRUE.
#' @param env_data_primary primary data frame of monthly sequences of environmental
#' data as columns and years as row names. Each row represents a year and
#' each column represents a day of a year. Row.names should be matched with
#' those from the response data frame. If not, set the argument row_names_subset
#' to TRUE. Alternatively, env_data_primary could be a tidy data with three columns,
#' i.e. Year, Month and third column representing values of mean temperatures,
#' sum of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data_primary to TRUE.
#' @param env_data_control a data frame of monthly sequences of environmental data as
#' columns and years as row names. This data is used as control for calculations of
#' partial correlation coefficients. Each row represents a year and each column
#' represents a day of a year. Row.names should be matched with those from the
#' response data frame. If not, set the row_names_subset argument to TRUE.
#' Alternatively, env_data_control could be a tidy data with three columns,
#' i.e. Year, Month and third column representing values of mean temperatures, sum
#' of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data_control to TRUE.
#' @param pcor_method a character string indicating which partial correlation
#' coefficient is to be computed. One of "pearson" (default), "kendall", or
#' "spearman", can be abbreviated.
#' @param lower_limit lower limit of window width (i.e. number of consecutive months
#' to be used for calculations)
#' @param upper_limit upper limit of window width (i.e. number of consecutive months
#' to be used for calculations)
#' @param fixed_width fixed width used for calculations (i.e. number of consecutive
#' months to be used for calculations)
#' @param previous_year logical. If TRUE, previous-year climate data are included
#' in the analysis. If FALSE, no previous-year climate data are included and
#' number_previous_years is ignored.
#' @param number_previous_years integer between 1 and 5 specifying how many
#' previous years should be included in the environmental matrices when
#' previous_year = TRUE. For example, number_previous_years = 2 uses climate
#' data from years t - 2, t - 1 and t for response year t. If NULL and
#' previous_year = TRUE, one previous year is included.
#' @param reference_window character string, the reference_window argument describes,
#' how each calculation is referred. There are two different options: 'start'
#' (default) and 'end'. If the reference_window argument is set to 'start',
#' then each calculation is related to the starting month of window. If the
#' reference_window argument is set to 'end', then each calculation is related
#' to the ending day of window calculation.
#' @param remove_insignificant if set to TRUE, removes all correlations bellow
#' the significant threshold level, based on a selected alpha. If
#' remove_insignificant_boot = TRUE and boot = TRUE, the bootstrap-based
#' filtering is used instead.
#' @param remove_insignificant_boot if set to TRUE and boot = TRUE, removes all
#' partial correlations for which the bootstrap confidence interval includes
#' zero. If both remove_insignificant and remove_insignificant_boot are TRUE,
#' the bootstrap confidence-interval filter takes priority. If boot = FALSE,
#' remove_insignificant_boot is ignored with a warning and remove_insignificant
#' is also ignored.
#' @param alpha significance level used to remove insignificant calculations.
#' @param row_names_subset if set to TRUE, row.names are used to subset
#' env_data_primary, env_data_control and response data frames. Only years from
#' all three data frames are kept.
#' @param aggregate_function_env_data_primary character string specifying how the
#' monthly data from env_data_primary should be aggregated. The default is 'mean',
#' the other options are 'median', 'sum' and 'quantile'.
#' @param quantile_prob_env_data_primary numeric value between 0 and 1 specifying
#' the quantile probability used when aggregate_function_env_data_primary =
#' 'quantile'. For example, 0.95 calculates the 95th percentile. The default is 0.5.
#' @param aggregate_function_env_data_control character string specifying how the
#' monthly data from env_data_control should be aggregated. The default is 'mean',
#' the other options are 'median', 'sum' and 'quantile'.
#' @param quantile_prob_env_data_control numeric value between 0 and 1 specifying
#' the quantile probability used when aggregate_function_env_data_control =
#' 'quantile'. For example, 0.95 calculates the 95th percentile. The default is 0.5.
#' @param temporal_stability_check character string, specifying, how temporal stability
#' between the optimal selection and response variable(s) will be analysed. Current
#' possibilities are "sequential", "progressive" and "running_window". Sequential check
#' will split data into k splits and calculate selected metric for each split. Progressive
#' check will split data into k splits, calculate metric for the first split and then
#' progressively add 1 split at a time and calculate selected metric. For running window,
#' select the length of running window with the k_running_window argument.
#' @param k_running_window the length of running window for temporal stability check.
#' Applicable only if temporal_stability argument is set to running window.
#' @param k integer, number of breaks (splits) for temporal stability
#' @param subset_years a subset of years to be analyzed. Should be given in the form of
#' subset_years = c(1980, 2005)
#' @param ylimits limit of the y axes for plot_extreme. It should be given in
#' the form of: ylimits = c(0,1)
#' @param seed optional seed argument for reproducible results
#' @param tidy_env_data_primary if set to TRUE, env_data_primary should be inserted as a
#' data frame with three columns: "Year", "Month", "Precipitation/Temperature/etc."
#' @param tidy_env_data_control if set to TRUE, env_data_control should be inserted as a
#' data frame with three columns: "Year", "Month", "Precipitation/Temperature/etc."
#' @param boot logical, if TRUE, bootstrap procedure will be used to calculate
#' partial correlation coefficients
#' @param boot_n The number of bootstrap replicates
#' @param boot_ci_type A character string representing the type of bootstrap intervals
#' required. The value should be any subset of the values c("norm","basic", "stud",
#' "perc", "bca").
#' @param boot_conf_int A scalar or vector containing the confidence level(s) of
#' the required interval(s)
#' @param month_interval a vector of two values defining the interval of months
#' used for calculations. Positive values indicate months in the current year.
#' Negative values indicate months in the previous-year block and are used only
#' when previous_year = TRUE. If previous_year = FALSE and negative values are
#' supplied, month_interval is ignored and the analysis is performed for the
#' current year only using month_interval = c(1, 12). If number_previous_years > 1,
#' the previous-year block starts with the earliest included previous year. For
#' example, previous_year = TRUE, number_previous_years = 2 and
#' month_interval = c(-1, 12) analyses the full sequence from January of year
#' t - 2 to December of year t.
#' @param dc_method a character string to determine the method to detrend climate
#' data.  Possible values are "none" (default) and "SLD" which refers to Simple
#' Linear Detrending
#' @param pcor_na_use an optional character string giving a method for computing
#' covariances in the presence of missing values for partial correlation
#' coefficients. This must be (an abbreviation of) one of the strings "all.obs",
#' "everything", "complete.obs", "na.or.complete", or "pairwise.complete.obs"
#' (default). See also the documentation for the base partial.r in psych R package
#'
#' @return a list with elements including:
#' \enumerate{
#'  \item $calculations - a matrix with calculated metrics
#'  \item $method - the character string of a method
#'  \item $metric - the character string indicating the metric used for calculations
#'  \item $analysed_period - the character string specifying the analysed period based on the information from row names. If there are no row names, this argument is given as NA
#'  \item $optimized_return - data frame with two columns, response variable and aggregated (averaged) monthly data that return the optimal results. This data.frame could be directly used to calibrate a model for climate reconstruction
#'  \item $optimized_return_all - a data frame with aggregated monthly data, that returned the optimal result for the entire env_data_primary (and not only subset of analysed years)
#'  \item $transfer_function - a ggplot object: scatter plot of optimized return and a transfer line of the selected method
#'  \item $temporal_stability - a data frame with calculations of selected metric for different temporal subsets
#'  \item $cross_validation - not available for partial correlation method
#'  \item $plot_heatmap - ggplot2 object: a heatmap of calculated metrics
#'  \item $plot_extreme - ggplot2 object: line plot of a row with the highest value in a matrix of calculated metrics
#'  \item $type - the character string describing type of analysis: monthly or monthly
#'  \item $reference_window - character string, which reference window was used for calculations
#'  \item $boot_lower - matrix with lower limits of bootstrap confidence intervals
#'  \item $boot_upper - matrix with upper limits of bootstrap confidence intervals
#'  \item $aggregated_climate_primary - matrix with all aggregated climate series of primary data
#'  \item $aggregated_climate_control - matrix with all aggregated climate series of control data
#'  \item $boot_ci_filter - list with information on bootstrap confidence-interval filtering
#'  \item $previous_year - logical indicating whether previous-year climate data were used
#'  \item $number_previous_years - integer indicating how many previous years were used
#'}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load the dendroTools R package
#' library(dendroTools)
#'
#' # Load data
#' data(data_MVA)
#' data(data_TRW)
#' data(data_TRW_1)
#' data(example_proxies_individual)
#' data(example_proxies_1)
#' data(LJ_monthly_temperatures)
#' data(LJ_monthly_precipitation)
#'
#' # 1 Basic example
#' example_basic <- monthly_response_seascorr(response = data_MVA,
#'    fixed_width = 11,
#'    env_data_primary = LJ_monthly_temperatures,
#'    env_data_control = LJ_monthly_precipitation,
#'    row_names_subset = TRUE,
#'    remove_insignificant = TRUE,
#'    reference_window = "start",
#'    aggregate_function_env_data_primary = 'median',
#'    aggregate_function_env_data_control = 'median',
#'    alpha = 0.05, pcor_method = "spearman",
#'    tidy_env_data_primary = FALSE,
#'    tidy_env_data_control = TRUE,
#'    previous_year = TRUE)
#'
#' # summary(example_basic)
#' # plot(example_basic, type = 1)
#' # plot(example_basic, type = 2)
#' # example_basic$optimized_return
#' # example_basic$optimized_return_all
#' # example_basic$temporal_stability
#'
#' # 2 Extended example
#' example_extended <- monthly_response_seascorr(response = data_MVA,
#'    env_data_primary = LJ_monthly_temperatures,
#'    env_data_control = LJ_monthly_precipitation,
#'    row_names_subset = TRUE, dc_method = "SLD",
#'    remove_insignificant = FALSE,
#'    aggregate_function_env_data_primary = 'mean',
#'    aggregate_function_env_data_control = 'mean',
#'    alpha = 0.05, pcor_na_use = "pairwise.complete",
#'    reference_window = "end",
#'    tidy_env_data_primary = FALSE,
#'    tidy_env_data_control = TRUE)
#'
#' # summary(example_extended)
#' # plot(example_extended, type = 1)
#' # plot(example_extended, type = 2)
#' # example_extended$optimized_return
#' # example_extended$optimized_return_all
#'
#' # 3 Example using quantiles for primary and control variables
#' example_quantile <- monthly_response_seascorr(
#'    response = data_MVA,
#'    env_data_primary = LJ_monthly_temperatures,
#'    env_data_control = LJ_monthly_precipitation,
#'    row_names_subset = TRUE,
#'    fixed_width = 3,
#'    previous_year = TRUE,
#'    aggregate_function_env_data_primary = "quantile",
#'    quantile_prob_env_data_primary = 0.50,
#'    aggregate_function_env_data_control = "quantile",
#'    quantile_prob_env_data_control = 0.50,
#'    tidy_env_data_primary = FALSE,
#'    tidy_env_data_control = TRUE)
#'
#'  # summary(example_quantile)
#'
#' }

monthly_response_seascorr <- function(response, env_data_primary, env_data_control,
                                      previous_year = FALSE, number_previous_years = NULL,
                                      pcor_method = "pearson",
                                      remove_insignificant = TRUE,
                                      remove_insignificant_boot = FALSE,
                                      lower_limit = 1,
                                      upper_limit = 12, fixed_width = 0,
                                      alpha = .05, row_names_subset = FALSE,
                                      reference_window = "start",
                                      aggregate_function_env_data_primary = 'mean',
                                      quantile_prob_env_data_primary = 0.5,
                                      aggregate_function_env_data_control = 'mean',
                                      quantile_prob_env_data_control = 0.5,
                                      temporal_stability_check = "sequential", k = 2,
                                      k_running_window = 30,
                                      subset_years = NULL,
                                      ylimits = NULL, seed = NULL, tidy_env_data_primary = FALSE,
                                      tidy_env_data_control = FALSE,  boot = FALSE, boot_n = 1000,
                                      boot_ci_type = "norm", boot_conf_int = 0.95,
                                      month_interval = NULL,
                                      dc_method = NULL,
                                      pcor_na_use = "pairwise.complete") {

  ##############################################################################
  # 1 month interval is organized

  months_per_year <- 12L

  if (reference_window == "middle"){

    stop(paste0("reference_window should be 'start' or 'end'. 'middle' reference_window is not implemented for the monthly_response_seascorr"))
  }

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

  year_block_width <- months_per_year * (number_previous_years + 1L)

  # Default interval:
  # current year only: c(1, 12)
  # previous-year analysis: c(-1, 12), i.e. from January of the earliest
  # included previous year to December of the current year.
  if (is.null(month_interval)) {
    month_interval <- if (isTRUE(previous_year)) c(-1, months_per_year) else c(1, months_per_year)
  }

  if (!is.numeric(month_interval) ||
      length(month_interval) != 2 ||
      any(is.na(month_interval))) {
    stop("month_interval must be a numeric vector of length 2.")
  }

  offset_start <- month_interval[1]
  offset_end <- month_interval[2]

  if (offset_start == 0 || offset_end == 0) {
    stop("month_interval cannot contain 0. Use negative values for previous-year months and positive values for current-year months.")
  }

  # Negative values in month_interval would imply previous-year climate data.
  # If previous_year = FALSE, previous-year data are not used. In this case,
  # ignore month_interval and use the full current year.
  if (isFALSE(previous_year) && (offset_start < 0 || offset_end < 0)) {

    warning(paste0("Negative values were supplied in month_interval, but ",
                   "previous_year = FALSE. The month_interval argument is ignored ",
                   "and the analysis is performed for the current year only ",
                   "using month_interval = c(1, 12)."))

    month_interval <- c(1, months_per_year)
    offset_start <- month_interval[1]
    offset_end <- month_interval[2]
  }

  # If both limits are positive while previous_year = TRUE, previous years are
  # not actually included in the selected interval. In that case, keep the
  # analysis current-year only.
  if (offset_start > 0 && offset_end > 0 && isTRUE(previous_year)) {

    number_previous_years <- 0L
    previous_year <- FALSE
    year_block_width <- months_per_year

    warning(paste0("Previous years are not included in selected month_interval. ",
                   "The analysis is performed for the current year only."))
  }

  # Convert the user-facing month_interval to column positions in the lagged matrix.
  # For number_previous_years = 2, the internal matrix is:
  # columns 1:12     = year t - 2
  # columns 13:24    = year t - 1
  # columns 25:36    = year t
  if (offset_start < 0 && offset_end < 0) {

    offset_start <- abs(offset_start)
    offset_end <- abs(offset_end)

  } else if (offset_start < 0 && offset_end > 0) {

    offset_start <- abs(offset_start)
    offset_end <- offset_end + months_per_year * number_previous_years

  } else if (offset_start > 0 && offset_end > 0) {

    # Current-year-only interval; already in current-year coordinates.
    offset_start <- offset_start
    offset_end <- offset_end

  } else {

    stop("month_interval must not run from current-year months back to previous-year months.")
  }

  if (offset_start > offset_end) {
    stop("month_interval is invalid after conversion. The start of the interval is after the end.")
  }

  if (offset_start < 1 || offset_end > year_block_width) {
    stop(paste0("month_interval is outside the available climate sequence. ",
                "For number_previous_years = ", number_previous_years,
                ", the available width is ", year_block_width, " months."))
  }

  # Calculate the max_window allowed
  max_window <- offset_end - offset_start + 1

  # If max_window is smaller than upper_limit, upper_limit must be reduced
  if (upper_limit > max_window){

    upper_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The upper_limit is outside your month_interval and",
                     " therefore reduced to the maximum allowed: ", max_window, "."))
    }
  }

  if (lower_limit > max_window){
    lower_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The lower_limit is outside your month_interval and",
                     " therefore reduced to the minimum allowed: ", max_window, "."))
    }
  }

  # Also correction for fixed_window approach - can only be ERROR
  if (fixed_width > max_window){

    stop(paste0("The selected fixed_width is outside your month_interval. ",
                "Decrease the fixed_width argument to at least: ", max_window, "."))
  }

  # offset_end is converted to the number of months that must be skipped at the
  # end of the full environmental matrix.
  offset_end <- year_block_width - offset_end

  ################################################################################

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.logical(remove_insignificant) ||
      length(remove_insignificant) != 1 ||
      is.na(remove_insignificant)) {
    stop("remove_insignificant must be either TRUE or FALSE.")
  }

  if (!is.logical(remove_insignificant_boot) ||
      length(remove_insignificant_boot) != 1 ||
      is.na(remove_insignificant_boot)) {
    stop("remove_insignificant_boot must be either TRUE or FALSE.")
  }

  boot_ci_filter <- list(
    applied = FALSE,
    criterion = NA_character_,
    removed_cells = NA_integer_,
    retained_cells = NA_integer_,
    total_cells = NA_integer_
  )

  if (isTRUE(remove_insignificant_boot) && isFALSE(boot)) {

    warning(paste0("remove_insignificant_boot = TRUE, but boot = FALSE. ",
                   "The bootstrap confidence-interval filter cannot be ",
                   "applied. Both remove_insignificant_boot and ",
                   "remove_insignificant are ignored."))

    remove_insignificant_boot <- FALSE
    remove_insignificant <- FALSE
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
  env_data_primary_control <- NULL
  metric <- NULL
  temporal_matrix_lower <- NULL
  temporal_matrix_upper <- NULL

  # Check selected aggregation functions
  allowed_aggregate_functions <- c("mean", "median", "sum", "quantile")

  check_aggregate_function <- function(aggregate_function, argument_name) {
    if (!(aggregate_function %in% allowed_aggregate_functions)) {
      stop(paste0(
        argument_name, " is '", aggregate_function,
        "'. Instead it should be one of: ",
        paste(allowed_aggregate_functions, collapse = ", "), "."
      ))
    }
  }

  check_quantile_prob <- function(quantile_prob, argument_name) {
    if (!is.numeric(quantile_prob) ||
        length(quantile_prob) != 1 ||
        is.na(quantile_prob) ||
        quantile_prob < 0 ||
        quantile_prob > 1) {
      stop(paste0(argument_name, " must be a single numeric value between 0 and 1."))
    }
  }

  check_aggregate_function(aggregate_function_env_data_primary,
                           "aggregate_function_env_data_primary")
  check_aggregate_function(aggregate_function_env_data_control,
                           "aggregate_function_env_data_control")
  check_quantile_prob(quantile_prob_env_data_primary,
                      "quantile_prob_env_data_primary")
  check_quantile_prob(quantile_prob_env_data_control,
                      "quantile_prob_env_data_control")

  # Internal helper for aggregating monthly data across rows
  aggregate_monthly_window <- function(x, aggregate_function, quantile_prob) {

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

  aggregate_monthly_window_primary <- function(x) {
    aggregate_monthly_window(x,
                             aggregate_function = aggregate_function_env_data_primary,
                             quantile_prob = quantile_prob_env_data_primary)
  }

  aggregate_monthly_window_control <- function(x) {
    aggregate_monthly_window(x,
                             aggregate_function = aggregate_function_env_data_control,
                             quantile_prob = quantile_prob_env_data_control)
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
      stop("For previous-year analyses, row.names of environmental data must be calendar years.")
    }

    if (any(duplicated(env_years))) {
      stop("Duplicated years are present in row.names of environmental data.")
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
      stop(paste0("No years in environmental data have the required ",
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

  # lower_limit = 1
  # upper_limit = 12
  # fixed_width = 0

  if (reference_window == "middle"){

    stop(paste0("reference_window should be 'start' or 'end'. 'middle' reference_window is not implemented for the monthly_response"))
  }

  # if fixed width is not 0, then change lower and upper limits, just to avoid
  # error messages
  if (fixed_width != 0){
    lower_limit = 1
    upper_limit = 12
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


  # If env_data_control is given in tidy version, transformation is needed
  if (tidy_env_data_control == TRUE){

    n_col_tidy_DF <- ncol(env_data_control)
    colnames_tidy_DF <- colnames(env_data_control)

    if (ncol(env_data_control) != 3){
      stop(paste("env_data_control was inserted in tidy version (tidy_env_data_control is set to TRUE).",
                 "env_data_control should have 3 columns, but it has", n_col_tidy_DF, "instead!"))
    }

    if (colnames_tidy_DF[1] != "Year"){
      stop(paste("env_data_control was inserted in tidy version (tidy_env_data_control is set to TRUE).",
                 "The first column name of the env_data_control should be 'Year', but it is",
                 colnames_tidy_DF[1], "instead!"))
    }

    if (colnames_tidy_DF[2] != "Month"){
      stop(paste("env_data_control was inserted in tidy version (tidy_env_data_control is set to TRUE).",
                 "The second column name of the env_data_control should be 'Month', but it is",
                 colnames_tidy_DF[2], "instead!"))
    }

    value_variable = colnames(env_data_control)[3]
    env_data_control <- dcast(env_data_control, Year~Month, value.var = value_variable)
    env_data_control <- years_to_rownames(env_data_control, "Year")
  }


  # If env_data_primary is given in tidy version, transformation is needed
  if (tidy_env_data_primary == TRUE){

    n_col_tidy_DF <- ncol(env_data_primary)
    colnames_tidy_DF <- colnames(env_data_primary)

    if (ncol(env_data_primary) != 3){
      stop(paste("env_data_primary was inserted in tidy version (tidy_env_data_primary is set to TRUE).",
                 "env_data_primary should have 3 columns, but it has", n_col_tidy_DF, "instead!"))
    }

    if (colnames_tidy_DF[1] != "Year"){
      stop(paste("env_data_primary was inserted in tidy version (tidy_env_data_primary is set to TRUE).",
                 "The first column name of the env_data_primary should be 'Year', but it is",
                 colnames_tidy_DF[1], "instead!"))
    }

    if (colnames_tidy_DF[2] != "Month"){
      stop(paste("env_data_primary was inserted in tidy version (tidy_env_data_primary is set to TRUE).",
                 "The second column name of the env_data_primary should be 'Month', but it is",
                 colnames_tidy_DF[2], "instead!"))
    }

    value_variable = colnames(env_data_primary)[3]
    env_data_primary <- dcast(env_data_primary, Year~Month, value.var = value_variable)
    env_data_primary <- years_to_rownames(env_data_primary, "Year")
  }



  # PART 1 - general data arrangements, warnings and stops
  # Both objects (response, env_data_primary and env_data_control) are converted
  # to data frames.
  response <- data.frame(response)
  env_data_primary <- data.frame(env_data_primary)
  env_data_control <- data.frame(env_data_control)

  # Here we save the original response data that will be used later
  response_original <- response

  # Previous-year analyses require alignment by calendar year names. This avoids
  # accidental row-position mismatches after adding lagged climate years.
  if (number_previous_years > 0L && row_names_subset == FALSE) {
    row_names_subset <- TRUE
    warning(paste0("Previous-year analyses require alignment by year names. ",
                   "row_names_subset is set to TRUE."))
  }

  # Build current-year or multi-year environmental matrices.
  # For number_previous_years = 0 these return the data unchanged.
  env_data_primary <- build_lagged_env_data(env_data_primary, number_previous_years)
  env_data_control <- build_lagged_env_data(env_data_control, number_previous_years)

  # These objects are used later for optimized_return_all. They should represent
  # the full lagged climate matrices before subset_years is applied.
  env_data_primary_original <- env_data_primary
  env_data_control_original <- env_data_control

  # For metric calculations, all objects need to have the same length,
  # with the exception when row_names_subset is set to TRUE.
  if (nrow(response) != nrow(env_data_primary) & row_names_subset == FALSE)
    stop(paste0("Length of env_data_primary and response records differ",
                " You can use row_names_subset = TRUE"))

  if (nrow(response) != nrow(env_data_control) & row_names_subset == FALSE)
    stop(paste0("Length of env_data_control and response records differ",
                " You can use row_names_subset = TRUE"))

  #######################################################
  # Rules for selected window limits

  if (fixed_width < 0 | fixed_width > year_block_width)
    stop(paste0("fixed_width should be between 0 and ", year_block_width))

  # Correlations could be calculated only for one variable
  if (ncol(response) > 1)
    stop(paste("More than 1 variable in response data frame not suitable ",
               "for 'pcor' method'"))

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

  # Warn users in case of missing values. The original threshold was 8 missing
  # months for one 12-month climate year, so it is scaled when multiple
  # previous years are included.
  missing_month_threshold <- 8 * (number_previous_years + 1L)

  # A) env_data_primary
  env_temp_primary <- env_data_primary[row.names(env_data_primary) %in% row.names(response),]

  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    subset_seq <- seq(lower_subset, upper_subset)
    env_temp_primary <- subset(env_temp_primary, row.names(env_temp_primary) %in% subset_seq)
  }

  na_problem <- data.frame(na_sum = rowSums(is.na(env_temp_primary)))
  na_problem <- na_problem[na_problem$na_sum > missing_month_threshold, , F]
  problematic_years <- paste0(row.names(na_problem), sep = "", collapse=", ")

  if (nrow(na_problem) > 0){

    warning(paste0("Problematic years with missing values are present in env_data_primary: ", problematic_years))

  }

  # B) env_data_control
  env_temp_conotrol <- env_data_control[row.names(env_data_control) %in% row.names(response),]

  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    subset_seq <- seq(lower_subset, upper_subset)
    env_temp_conotrol <- subset(env_temp_conotrol, row.names(env_temp_conotrol) %in% subset_seq)
  }

  na_problem <- data.frame(na_sum = rowSums(is.na(env_temp_conotrol)))
  na_problem <- na_problem[na_problem$na_sum > missing_month_threshold, , F]
  problematic_years <- paste0(row.names(na_problem), sep = "", collapse=", ")

  if (nrow(na_problem) > 0){

    warning(paste0("Problematic years with missing values are present in env_data_control: ", problematic_years))

  }

  ##############################################################################

  # Previous-year climate data have already been added by
  # build_lagged_env_data().

  # If row_names_subset == TRUE, data is subset and ordered based on matching
  # row.names. Additionally, number of characters in row.names is checked.
  # There should be at least three characters (assuming years before 100 will
  # never be analysed, there is no such environmental data available)
  if (row_names_subset == TRUE & nchar(row.names(env_data_primary)[1]) >= 3){

    env_data_primary$yearABC <- row.names(env_data_primary)
    env_data_control$yearABC <- row.names(env_data_control)
    response$yearABC <- row.names(response)

    intersects <- intersect(intersect(response$yearABC, env_data_control$yearABC), env_data_primary$yearABC)

    response <- dplyr::filter(response, yearABC %in% intersects)
    response <- dplyr::arrange(response, -as.numeric(yearABC))
    response <- years_to_rownames(response, "yearABC")

    env_data_primary <- dplyr::filter(env_data_primary, yearABC %in% intersects)
    env_data_primary <- dplyr::arrange(env_data_primary, -as.numeric(yearABC))
    env_data_primary <- years_to_rownames(env_data_primary, "yearABC")

    env_data_control <- dplyr::filter(env_data_control, yearABC %in% intersects)
    env_data_control <- dplyr::arrange(env_data_control, -as.numeric(yearABC))
    env_data_control <- years_to_rownames(env_data_control, "yearABC")

  }

  # if row.names of env_data_primary and the response data frames are not equal,
  # warning is given.
  if (sum(row.names(env_data_primary) == row.names(response)) != nrow(env_data_primary)) {
    warning("row.names between env_data_primary and response do not match!")
  }

  # If row_names_subset == TRUE, but row.names does not appear to be years,
  # error is given.
  if (row_names_subset == TRUE & nchar(row.names(env_data_primary)[1]) < 3){
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

    if (any(!(subset_seq %in% row.names(env_data_primary)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  " in the env_data_primary data frame. Change the subset_years argument"))
    }

    if (any(!(subset_seq %in% row.names(env_data_control)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  "in the env_data_control data frame. Change the subset_years argument"))
    }

    response <- subset(response, row.names(response) %in% subset_seq)
    env_data_primary <- subset(env_data_primary, row.names(env_data_primary) %in% subset_seq)
    env_data_control <- subset(env_data_control, row.names(env_data_control) %in% subset_seq)
  }

  # NA values are not allowed and must be removed from response data.frame
  # exception if cor_na_us accounts for missing values
  if (sum(is.na(response)) > 0 & !(pcor_na_use %in% c("complete.obs", "na.or.complete", "pairwise.complete.obs"))){

    prob_year <- row.names(response[is.na(response), , drop = F])

    stop(paste0("NA is not allowed in response data frame. ",
                "Problematic year is ", prob_year))

  }


  # PART 2 - the calculation of day-wise correlations
  # A) fixed window approach
  # B) Variable window approach

  # these are lists for climate and and holder for saving mm1 and mm2
  list_climate_primary <- list()
  list_climate_control <- list()

  mm1 <- 1
  mm2 <- 1

  # A) The fixed window approach
  if (fixed_width != 0) {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be saved in this matrix
    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data_primary) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data_primary)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data_primary) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }

    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (fixed_width != max_window){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (ncol(env_data_primary) - fixed_width - offset_end - offset_start + 1),
                             style = 3)}

    }

    b = 0

    # An iterating loop. In each iteration x is calculated and represents
    # response (dependent) variable. X is a moving average. Window width of
    # a moving window is fixed_width. Next, partial correlation is stored in
    # temporal matrix.

    for (j in (0 + offset_start -1): (ncol(env_data_primary) - max((fixed_width + offset_end), offset_end))) {

      b = b + 1

      x1 <- aggregate_monthly_window_primary(
        env_data_primary[1:nrow(env_data_primary),
                         (1 + j):(j + fixed_width),
                         drop = FALSE]
      )

      x2 <- aggregate_monthly_window_control(
        env_data_control[1:nrow(env_data_control),
                         (1 + j):(j + fixed_width),
                         drop = FALSE]
      )

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          tmp_model <- lm(x1 ~ seq(1:length(x1)))
          tmp_pred <- predict(tmp_model)

          if (length(x1) != length(tmp_pred)) {
            warning("Note missing values in your env_data")
          }

          tmp_res <- suppressWarnings(x1 - tmp_pred)

          x1 <- data.frame(x1 = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x1 <- matrix(x1, nrow = nrow(env_data_primary), ncol = 1)

      }

      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          tmp_model <- lm(x2 ~ seq(1:length(x2)))
          tmp_pred <- predict(tmp_model)

          if (length(x2) != length(tmp_pred)) {
            warning("Note missing values in your env_data")
          }

          tmp_res <- suppressWarnings(x2 - tmp_pred)

          x2 <- data.frame(x2 = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x2 <- matrix(x2, nrow = nrow(env_data_primary), ncol = 1)

      }

      my_temporal_data <- cbind(response[, 1], x1, x2)
      colnames(my_temporal_data) <- c("x", "y", "z")

      x1_list <- x1
      colnames(x1_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x1_list) <- row.names(env_data_primary)
      list_climate_primary[[mm1]] <- x1_list
      mm1 = mm1 + 1

      x2_list <- x2
      colnames(x2_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x2_list) <- row.names(env_data_control)
      list_climate_control[[mm2]] <- x2_list
      mm2 = mm2 + 1

      if (boot == FALSE){

        temporal_correlation <- partial.r(data=my_temporal_data, x=c("x","y"), y="z",
                                          use=pcor_na_use,method = pcor_method)[2]

        if (reference_window == 'start'){
          temporal_matrix[1, j + 1] <- as.numeric(temporal_correlation)[1]
        } else if (reference_window == 'end'){
          temporal_matrix[1, j + fixed_width] <- as.numeric(temporal_correlation)[1]
        } else if (reference_window == 'middle'){
          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- as.numeric(temporal_correlation)[1]
        }

      } else if (boot == TRUE){

        calc <- boot(data = my_temporal_data, statistic = boot_f_pcor, R = boot_n, cor.type = pcor_method)

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

      } else {
        print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }

      if(interactive()){

        if (fixed_width != max_window){setTxtProgressBar(pb, b)}

      }
    }

    if(interactive()){

      if (fixed_width != max_window){close(pb)}

    }

    # temporal_matrix is given rownames and colnames. Rownames represent a
    # window width used for calculations. Colnames represent the position of
    # moving window in a original env_data_primary data frame.
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

  if (fixed_width == 0) {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be saved in this matrix

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data_primary) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data_primary)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data_primary) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }

    # An iterating double loop: 1 outer loop) iterating from lower_limit to
    # upper_limit defines window width used for a moving window. 2) inner loop
    # defines the starting position of a moving window.
    # In each iteration, x is calculated and represents a response (dependent)
    # variable. x is a moving average, based on rowMeans/apply function.
    # Next, statistical metric is calculated based on a selected method (pcor).
    # Calculation is stored in temporal matrix in a proper place. The position of
    # stored calculation is informative later used for indicating optimal values.

    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (upper_limit != lower_limit){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)

      }
    }

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1

      for (j in (0 + offset_start -1): (ncol(env_data_primary) - max((K + offset_end), offset_end))) {

        x1 <- aggregate_monthly_window_primary(
          env_data_primary[1:nrow(env_data_primary),
                           (1 + j):(j + K),
                           drop = FALSE]
        )

        x2 <- aggregate_monthly_window_control(
          env_data_control[1:nrow(env_data_control),
                           (1 + j):(j + K),
                           drop = FALSE]
        )

        if (!is.null(dc_method)){

          if (dc_method == "SLD"){

            tmp_model <- lm(x1 ~ seq(1:length(x1)))
            tmp_pred <- predict(tmp_model)

            if (length(x1) != length(tmp_pred)) {
              warning("Note missing values in your env_data")
            }

            tmp_res <- suppressWarnings(x1 - tmp_pred)

            x1 <- data.frame(x1 = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x1 <- matrix(x1, nrow = nrow(env_data_primary), ncol = 1)

        }

        if (!is.null(dc_method)){

          if (dc_method == "SLD"){

            tmp_model <- lm(x2 ~ seq(1:length(x2)))
            tmp_pred <- predict(tmp_model)

            if (length(x2) != length(tmp_pred)) {
              warning("Note missing values in your env_data")
            }

            tmp_res <- suppressWarnings(x2 - tmp_pred)

            x2 <- data.frame(x2 = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x2 <- matrix(x2, nrow = nrow(env_data_primary), ncol = 1)

        }

        my_temporal_data <- cbind(response[, 1], x1, x2)
        colnames(my_temporal_data) <- c("x", "y", "z")

        x1_list <- x1
        colnames(x1_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x1_list) <- row.names(env_data_primary)
        list_climate_primary[[mm1]] <- x1_list
        mm1 = mm1 + 1

        x2_list <- x2
        colnames(x2_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x2_list) <- row.names(env_data_control)
        list_climate_control[[mm2]] <- x2_list
        mm2 = mm2 + 1

        if (boot == FALSE){

          temporal_correlation <- partial.r(data=my_temporal_data, x=c("x","y"), y="z",
                                            use=pcor_na_use,method = pcor_method)[2]

          if (reference_window == 'start'){

            temporal_matrix[(K - lower_limit) + 1, j + 1] <- as.numeric(temporal_correlation)[1]

          } else if (reference_window == 'end'){

            temporal_matrix[(K - lower_limit) + 1, j + K] <- as.numeric(temporal_correlation)[1]

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- as.numeric(temporal_correlation)[1]

          }

        } else if (boot == TRUE){

          calc <- boot(data = my_temporal_data, statistic = boot_f_pcor, R = boot_n, cor.type = pcor_method)

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

        } else {
          print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
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
    # moving window in a original env_data_primary data frame.
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

  # To enhance the visualization, insignificant values can be removed
  # either using the standard alpha-based threshold or, when requested, using
  # bootstrap confidence intervals. The bootstrap filter takes priority when
  # both remove_insignificant and remove_insignificant_boot are TRUE.
  if (isTRUE(remove_insignificant_boot) && isTRUE(boot)) {

    if (!all(dim(temporal_matrix) == dim(temporal_matrix_lower)) ||
        !all(dim(temporal_matrix) == dim(temporal_matrix_upper))) {
      stop("temporal_matrix, boot_lower and boot_upper must have the same dimensions.")
    }

    boot_keep <- (temporal_matrix_lower > 0 & temporal_matrix_upper > 0) |
      (temporal_matrix_lower < 0 & temporal_matrix_upper < 0)

    boot_keep[is.na(boot_keep)] <- FALSE

    temporal_matrix[!boot_keep] <- NA

    boot_ci_filter <- list(
      applied = TRUE,
      criterion = "Bootstrap confidence interval excludes zero",
      removed_cells = sum(!boot_keep, na.rm = TRUE),
      retained_cells = sum(boot_keep, na.rm = TRUE),
      total_cells = length(boot_keep)
    )

  } else if (isTRUE(remove_insignificant)) {

    critical_threshold_cor <- critical_r(nrow(response), alpha = alpha)
    critical_threshold_cor2 <- critical_threshold_cor ^ 2

    temporal_matrix[abs(temporal_matrix) < abs(critical_threshold_cor)] <- NA
  }

  if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){

    warning("All calculations are insignificant! Please change the alpha argument.")

    final_list <- list(calculations = temporal_matrix,
                       method = "pcor",
                       metric = pcor_method, analysed_period = NA,
                       optimized_return = NA,
                       optimized_return_all = NA,
                       transfer_function = NA, temporal_stability = NA,
                       cross_validation = NA,
                       plot_heatmap = NA,
                       plot_extreme = NA,
                       type = "monthly",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate_primary = NA,
                       aggregated_climate_control = NA,
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       boot_ci_filter = boot_ci_filter)

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
      max_calculation <- max(temporal_matrix, na.rm = TRUE)

      max_result <- suppressWarnings(which.max(apply(temporal_matrix,
                                                     MARGIN = 2, max,
                                                     na.rm = TRUE)))
      plot_column <- max_result
      max_index <- which.max(temporal_matrix[, names(max_result)])
      row_index <- row.names(temporal_matrix)[max_index]
    }

    if ((abs(overall_max) < abs(overall_min)) == TRUE) {

      max_calculation <- min(temporal_matrix, na.rm = TRUE)

      min_result <- suppressWarnings(which.min(apply(temporal_matrix,
                                                     MARGIN = 2, min,
                                                     na.rm = TRUE)))
      plot_column <- min_result
      min_index <- which.min(temporal_matrix[, names(min_result)])
      row_index <- row.names(temporal_matrix)[min_index]
    }


    ###############################################################

    # The fourth return element is being created: rowMeans/apply of optimal sequence.
    # Here we consider options based on reference_window.

    # 1. reference window = "start"
    if (reference_window == 'start'){

      x1 <- data.frame(aggregate_monthly_window_primary(
        env_data_primary[, as.numeric(plot_column):
                           (as.numeric(plot_column) + as.numeric(row_index) - 1),
                         drop = FALSE]
      ))

      x2 <- data.frame(aggregate_monthly_window_control(
        env_data_control[, as.numeric(plot_column):
                           (as.numeric(plot_column) + as.numeric(row_index) - 1),
                         drop = FALSE]
      ))

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      x1_original <- data.frame(aggregate_monthly_window_primary(
        env_data_primary_original[, as.numeric(plot_column):
                                    (as.numeric(plot_column) + as.numeric(row_index) - 1),
                                  drop = FALSE]
      ))

      x2_original <- data.frame(aggregate_monthly_window_control(
        env_data_control_original[, as.numeric(plot_column):
                                    (as.numeric(plot_column) + as.numeric(row_index) - 1),
                                  drop = FALSE]
      ))

    }

    # Option 2, reference window = "end"
    if (reference_window == 'end'){

      x1 <- data.frame(aggregate_monthly_window_primary(
        env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                           as.numeric(plot_column),
                         drop = FALSE]
      ))

      x2 <- data.frame(aggregate_monthly_window_control(
        env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                           as.numeric(plot_column),
                         drop = FALSE]
      ))

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      x1_original <- data.frame(aggregate_monthly_window_primary(
        env_data_primary_original[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                    as.numeric(plot_column),
                                  drop = FALSE]
      ))

      x2_original <- data.frame(aggregate_monthly_window_control(
        env_data_control_original[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                    as.numeric(plot_column),
                                  drop = FALSE]
      ))

    }

    # Third option: reference window = "middle"
    if (reference_window == 'middle'){

      if (as.numeric(row_index) %% 2 == 0){
        adjustment_1 <- 0
        adjustment_2 <- 1
      } else {
        adjustment_1 <- 1
        adjustment_2 <- 2
      }

      x1 <- data.frame(aggregate_monthly_window_primary(
        env_data_primary[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                           (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                         drop = FALSE]
      ))

      x2 <- data.frame(aggregate_monthly_window_control(
        env_data_control[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                           (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                         drop = FALSE]
      ))

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      x1_original <- data.frame(aggregate_monthly_window_primary(
        env_data_primary_original[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                                    (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                                  drop = FALSE]
      ))

      x2_original <- data.frame(aggregate_monthly_window_control(
        env_data_control_original[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                                    (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                                  drop = FALSE]
      ))

    }

    # Final output list

    if (!is.null(dc_method)){

      if (dc_method == "SLD"){

        x1 <- x1[,1]
        tmp_model <- lm(x1 ~ seq(1:length(x1)))
        tmp_pred <- predict(tmp_model)

        if (length(x1) != length(tmp_pred)) {
          warning("Note missing values in your env_data")
        }

        tmp_res <- suppressWarnings(x1 - tmp_pred)

        x1 <- data.frame(x1 = tmp_res/sd(tmp_res, na.rm = TRUE))

        x2 <- x2[,1]
        tmp_model <- lm(x2 ~ seq(1:length(x2)))
        tmp_pred <- predict(tmp_model)

        if (length(x2) != length(tmp_pred)) {
          warning("Note missing values in your env_data")
        }

        tmp_res <- x2 - tmp_pred

        x2 <- data.frame(x2 = tmp_res/sd(tmp_res, na.rm = TRUE))

      }
    }


    x1_full <- cbind(response, x1, x2)
    colnames(x1_full)[ncol(x1_full)-1] <- "Optimized_return_primary"
    colnames(x1_full)[ncol(x1_full)] <- "Optimized_return_control"

    colnames(x1) <- "Optimized.rowNames.primary"
    colnames(x2) <- "Optimized.rowNames.control"


    # original

    if (!is.null(dc_method)){

      if (dc_method == "SLD"){

        x1_original <- x1_original[,1]
        tmp_model <- lm(x1_original ~ seq(1:length(x1_original)))
        tmp_pred <- predict(tmp_model)
        tmp_res <- suppressWarnings(x1_original - tmp_pred)
        x1_original <- data.frame(x1_original = tmp_res/sd(tmp_res, na.rm = TRUE))

        x2_original <- x2_original[,1]
        tmp_model <- lm(x2_original ~ seq(1:length(x2_original)))
        tmp_pred <- predict(tmp_model)
        tmp_res <- suppressWarnings(x2_original - tmp_pred)
        x2_original <- data.frame(x2_original = tmp_res/sd(tmp_res, na.rm = TRUE))

      }
    }

    x1_full_original <- merge(x1_original, x2_original, by = 0, all = TRUE)
    colnames(x1_full_original)[1] <- "year"
    x1_full_original <- years_to_rownames(x1_full_original, "year")
    colnames(x1_full_original)[ncol(x1_full_original)-1] <- "Optimized_return_primary"
    colnames(x1_full_original)[ncol(x1_full_original)] <- "Optimized_return_control"

    colnames(x1) <- "Optimized.rowNames"

    my_temporal_data <- cbind(x1_full[,1], x1_full[,2], x1_full[,3])
    colnames(my_temporal_data) <- c("x", "y", "z")
    test_calculation <- partial.r(data=my_temporal_data, x=c("x","y"), y="z", use=pcor_na_use,method = pcor_method)[2]

    test_logical <- as.numeric(max_calculation) == as.numeric(test_calculation)

    # Element 5
    # Here we create the fifth element of the final list: Analysed period in the
    # form of min(year) - max(year), e.g. 1950 - 2015
    min_env_data_primary <- min(as.numeric(row.names(env_data_primary)))
    min_response <- min(as.numeric(row.names(response)))

    max_env_data_primary <- max(as.numeric(row.names(env_data_primary)))
    max_response <- max(as.numeric(row.names(response)))

    min_together <- min(min_env_data_primary, min_response)
    max_together <- min(max_env_data_primary, max_response)


    analysed_period <- paste(as.character(min_together),
                             as.character(max_together),
                             sep = " - ")
    if (nchar(analysed_period) < 9) {
      analysed_period <- NA
    }

    # Here, the transfer function is being created
    transfer_data = data.frame(proxy = response[,1], optimized_return =x1[,1])
    lm_model = lm(optimized_return ~ proxy, data = transfer_data)
    full_range = data.frame(proxy = seq(from = min(response[,1], na.rm = TRUE), to = max(response[,1], na.rm = TRUE), length.out = 100))
    full_range$transfer_f = predict(lm_model, full_range)

    # String for titles

    title_string <- "Partial Correlation Coefficients"

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

    dataset = x1_full

    empty_list = list()
    empty_list_period = list()
    empty_list_significance = list()

    temporal_stability <- data.frame()

    # 1. Progressive stability check
    if (temporal_stability_check == "progressive"){

      foldi <- seq(1:k)
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      for (m in 1:k){

        #Segement your data by fold using the which() function
        trainIndexes <- which(folds <= m, arr.ind = TRUE)
        dataset_temp <- dataset[trainIndexes, ]
        MAKS <- max(as.numeric(row.names(dataset_temp)))
        MIN <- min(as.numeric(row.names(dataset_temp)))
        empty_list_period[[m]] <- paste(MIN, "-", MAKS)

        par.r <- partial.r(data=dataset_temp, x=c(1,2), y=3, use=pcor_na_use,method = pcor_method)
        calculation <- par.r[2]
        sig <- corr.p(par.r,n = nrow(dataset_temp - 2))$p[2]

        empty_list_significance[[m]] <- sig
        empty_list[[m]] <- calculation
        colname = "partial correlation"

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
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      for (m in 1:k){

        #Segement your data by fold using the which() function
        trainIndexes <- which(folds == m, arr.ind = TRUE)
        dataset_temp <- dataset[trainIndexes, ]

        MAKS <- max(as.numeric(row.names(dataset_temp)))
        MIN <- min(as.numeric(row.names(dataset_temp)))
        empty_list_period[[m]] <- paste(MIN, "-", MAKS)

        par.r <- partial.r(data=dataset_temp, x=c(1,2), y=3, use=pcor_na_use,method = pcor_method)
        calculation <- par.r[2]
        sig <- corr.p(par.r,n = nrow(dataset_temp - 2))$p[2]

        empty_list_significance[[m]] <- sig
        empty_list[[m]] <- calculation
        colname = "partial correlation"

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

        par.r <- partial.r(data=dataset_temp, x=c(1,2), y=3, use=pcor_na_use,method = pcor_method)
        calculation <- par.r[2]
        sig <- corr.p(par.r,n = nrow(dataset_temp - 2))$p[2]

        empty_list_significance[[m]] <- sig
        empty_list[[m]] <- calculation
        colname = "partial correlation"

      }
      m1 <- do.call(rbind, empty_list)
      m2 <- do.call(rbind, empty_list_period)
      m3 <- do.call(rbind, empty_list_significance)

      temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(as.numeric(m3), digits = 3), nsmall = 3)))
      colnames(temporal_stability) <-c("Period", colname, "p value")
      temporal_stability
    }

    ################################################################
    #### Here the final list is being filled with six elements #####
    ################################################################

    final_list <- list(calculations = temporal_matrix, method = "pcor",
                       metric = pcor_method, analysed_period = analysed_period,
                       optimized_return = x1_full,
                       optimized_return_all = x1_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = NA,
                       aggregated_climate_primary = do.call(cbind, list_climate_primary),
                       aggregated_climate_control = do.call(cbind, list_climate_control),
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       boot_ci_filter = boot_ci_filter)

    plot_heatmapA <- plot_heatmap(final_list, reference_window = reference_window, type = "monthly")
    plot_extremeA <- plot_extreme(final_list, ylimits = ylimits, reference_window = reference_window, type = "monthly")

    final_list <- list(calculations = temporal_matrix, method = "pcor",
                       metric = pcor_method, analysed_period = analysed_period,
                       optimized_return = x1_full,
                       optimized_return_all = x1_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = NA,
                       plot_heatmap = plot_heatmapA,
                       plot_extreme = plot_extremeA,
                       type = "monthly",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate_primary = do.call(cbind, list_climate_primary),
                       aggregated_climate_control = do.call(cbind, list_climate_control),
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       boot_ci_filter = boot_ci_filter)

    class(final_list) <- "dmrs"

    return(final_list)

  }
}
