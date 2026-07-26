#' daily_response_seascorr
#'
#' Function calculates all possible partial correlation coefficients between
#' tree-ring chronology and daily environmental (usually climate) data.
#' Calculations are based on moving window which is defined with two
#' arguments: lower_limit and upper_limit. All calculated (partial) correlation
#' coefficients are stored in a matrix. The location of stored correlation
#' in the matrix is indicating a window width (row names) and a location in a
#' matrix of daily sequences of environmental data (column names).
#' @param response a data frame with tree-ring proxy variable and (optional)
#' years as row names. Row.names should be matched with those from env_data_primary
#' and env_data_control data frame. If not, set the row_names_subset argument to
#' TRUE.
#' @param env_data_primary primary data frame of daily sequences of environmental
#' data as columns and years as row names. Each row represents a year and
#' each column represents a day of a year. Row.names should be matched with
#' those from the response data frame. If not, set the argument row_names_subset
#' to TRUE. Alternatively, env_data_primary could be a tidy data with three columns,
#' i.e. Year, DOY and third column representing values of mean temperatures,
#' sum of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data_primary to TRUE.
#' @param env_data_control a data frame of daily sequences of environmental data as
#' columns and years as row names. This data is used as control for calculations of
#' partial correlation coefficients. Each row represents a year and each column
#' represents a day of a year. Row.names should be matched with those from the
#' response data frame. If not, set the row_names_subset argument to TRUE.
#' Alternatively, env_data_control could be a tidy data with three columns,
#' i.e. Year, DOY and third column representing values of mean temperatures, sum
#' of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data_control to TRUE.
#' @param lower_limit lower limit of window width
#' @param upper_limit upper limit of window width
#' @param fixed_width fixed width used for calculation. If fixed_width is
#' assigned a value, upper_limit and lower_limit will be ignored
#' @param pcor_method a character string indicating which partial correlation
#' coefficient is to be computed. One of "pearson" (default), "kendall", or
#' "spearman", can be abbreviated.
#' @param previous_year logical. Master switch for previous-year analyses. If
#' FALSE, only the current year is analysed and number_previous_years is ignored.
#' If TRUE, previous-year climate data are included.
#' @param number_previous_years integer between 1 and 5 specifying how many
#' previous years should be included when previous_year = TRUE. For example,
#' number_previous_years = 2 uses climate data from years t - 2, t - 1 and t
#' for response year t. If previous_year = TRUE and number_previous_years is
#' NULL, one previous year is used for backward compatibility.
#' @param remove_insignificant if set to TRUE, removes all correlations bellow
#' the significant threshold level, based on a selected alpha.
#' @param remove_insignificant_method the criterion used when
#' remove_insignificant = TRUE. One of "alpha" (per-window significance
#' threshold; the default), "boot" (bootstrap confidence interval excludes
#' zero; requires boot = TRUE), "fdr_BY" (false discovery rate,
#' Benjamini-Yekutieli, valid under arbitrary dependence between windows),
#' "fdr_BH" (Benjamini-Hochberg, valid under positive dependence and less
#' conservative), "maxstat" (max-statistic permutation test), "sidak"
#' (Sidak correction on the effective number of tests) or "cluster"
#' (cluster-extent inference). The "alpha" and "boot" criteria are applied per
#' window and therefore do not account for the number of windows screened; the
#' remaining criteria do.
#' Degrees of freedom account for the control variable throughout.
#'
#' "maxstat" builds a null distribution of the largest absolute partial
#' correlation obtained anywhere in the search by repeatedly shuffling the
#' response, and retains only windows exceeding the resulting threshold. It
#' makes no distributional assumptions but requires perm_n re-evaluations of
#' the whole search.
#'
#' "sidak" targets the same family-wise error rate but estimates it
#' analytically and requires no permutations. The effective number of
#' independent tests among the screened primary windows is obtained from the
#' eigenvalues of their correlation matrix following Li and Ji (2005), and a
#' Sidak correction is applied to that number rather than to the nominal count
#' of windows.
#'
#' "sidak" is an approximation to "maxstat", and a rougher one for partial
#' correlations than for ordinary ones, because the null distribution also
#' depends on the control series. Across 54 combinations of site and
#' primary/control pairing the two thresholds agreed to within 0.020 on
#' average, with a maximum discrepancy of 0.046. The approximation is close to
#' exact when the primary variable is strongly autocorrelated (mean error
#' -0.005 for temperature controlling for precipitation) and permissive when
#' it is not (mean error -0.031 for precipitation controlling for
#' temperature), because the eigenvalue approach does not represent temporal
#' autocorrelation. "maxstat" is therefore preferable for final results, while
#' "sidak" suits exploratory work.
#'
#' "cluster" is cluster-extent inference, adapted from neuroimaging. The
#' surface is thresholded at cluster_threshold, connected regions of
#' supra-threshold windows are identified, and each region's mass, the sum of
#' absolute partial correlations inside it, is compared against the largest
#' mass arising in permuted surfaces. It controls the family-wise error rate
#' over CLUSTERS rather than over windows, so it establishes that a region
#' contains signal without identifying which windows within it carry it.
#'
#' In tree-ring response surfaces this method has little discriminating power.
#' A real climate signal forms one large connected region whose mass greatly
#' exceeds anything the permuted surfaces produce, so every cluster that forms
#' is usually retained and the permutation step rarely rejects. The number of
#' windows kept is then determined by cluster_threshold rather than by the
#' data. It is provided for completeness and comparison; "maxstat" or "sidak"
#' are preferable when the question is which windows carry signal.
#' @param cluster_threshold the absolute partial correlation at which clusters
#' are formed when remove_insignificant_method = "cluster". Either "nominal",
#' the default, which uses the critical correlation at the chosen alpha, or a
#' single number between 0 and 1. This is a free parameter and it largely
#' determines how many windows are retained, so the value used is reported in
#' $significance_filter. Ignored by all other criteria.
#' @param perm_n number of permutations used when
#' remove_insignificant_method = "maxstat". Defaults to 999. Ignored by all
#' other criteria.
#' @param alpha significance level used to remove insignificant calculations.
#' @param row_names_subset if set to TRUE, row.names are used to subset
#' env_data_primary, env_data_control and response data frames. Only years from
#' all three data frames are kept.
#' @param aggregate_function_env_data_primary character string specifying how the
#' daily data from env_data_primary should be aggregated. The default is 'mean',
#' the other options are 'median', 'sum', 'min', 'max' and 'quantile'.
#' @param quantile_prob_env_data_primary numeric value between 0 and 1 specifying
#' the quantile probability used when aggregate_function_env_data_primary =
#' 'quantile'. For example, 0.95 calculates the 95th percentile. The default is 0.5.
#' @param aggregate_function_env_data_control character string specifying how the
#' daily data from env_data_control should be aggregated. The default is 'mean',
#' the other options are 'median', 'sum', 'min', 'max' and 'quantile'.
#' @param quantile_prob_env_data_control numeric value between 0 and 1 specifying
#' the quantile probability used when aggregate_function_env_data_control =
#' 'quantile'. For example, 0.95 calculates the 95th percentile. The default is 0.5.
#' @param temporal_stability_check character string, specifying how temporal stability
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
#' @param ylimits limit of the y axes for plot_extreme It should be given in the
#' form of: ylimits = c(0,1)
#' @param seed optional seed argument for reproducible results
#' @param tidy_env_data_primary if set to TRUE, env_data_primary should be inserted as a
#' data frame with three columns: "Year", "DOY", "Precipitation/Temperature/etc."
#' @param tidy_env_data_control if set to TRUE, env_data_control should be inserted as a
#' data frame with three columns: "Year", "DOY", "Precipitation/Temperature/etc."
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
#' partial correlation coefficients
#' @param boot_n The number of bootstrap replicates
#' @param boot_ci_type A character string representing the type of bootstrap intervals
#' required. The value should be any subset of the values c("norm","basic", "stud",
#' "perc", "bca").
#' @param boot_conf_int A scalar or vector containing the confidence level(s) of
#' the required interval(s)
#' @param day_interval a vector of two values defining the interval of days used
#' for calculations. Positive values indicate days in the current year. Negative
#' values indicate days in the previous-year block. If previous_year = FALSE,
#' negative values in day_interval are ignored with a warning and the analysis
#' is performed for the current year only using day_interval = c(1, 366). If
#' previous_year = TRUE and number_previous_years > 1, the previous-year block
#' starts with the earliest included previous year. For example, previous_year
#' = TRUE, number_previous_years = 2 and day_interval = c(-1, 366) analyses
#' the full sequence from DOY 1 of year t - 2 to DOY 366 of year t.
#' @param dc_method a character string to determine the method to detrend climate
#' data.  Possible values are "none" (default) and "SLD" which refers to Simple
#' Linear Detrending
#' @param pcor_na_use an optional character string giving a method for computing
#' covariances in the presence of missing values for partial correlation
#' coefficients. This must be (an abbreviation of) one of the strings "all.obs",
#' "everything", "complete.obs", "na.or.complete", or "pairwise.complete.obs"
#' (default). See also the documentation for the base partial.r in psych R package
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
#'  \item $optimized_return_all - a data frame with aggregated daily data, that returned the optimal result for the entire env_data_primary (and not only subset of analysed years)
#'  \item $transfer_function - a ggplot object: scatter plot of optimized return and a transfer line of the selected method
#'  \item $temporal_stability - a data frame with calculations of selected metric for different temporal subsets
#'  \item $cross_validation - not available for partial correlations
#'  \item $plot_heatmap - ggplot2 object: a heatmap of calculated metrics
#'  \item $plot_extreme - ggplot2 object: line plot of a row with the highest value in a matrix of calculated metrics
#'  \item $type - the character string describing type of analysis: daily or monthly
#'  \item $reference_window - character string, which reference window was used for calculations
#'  \item $aggregated_climate_primary - matrix with all aggregated climate series of primary data
#'  \item $aggregated_climate_control - matrix with all aggregated climate series of control data
#'  \item $previous_year - logical indicating whether previous-year climate data were used
#'  \item $number_previous_years - integer indicating how many previous years were used
#'  \item $boot_lower - matrix with lower limit of bootstrap confidence intervals
#'  \item $boot_upper - matrix with upper limit of bootstrap confidence intervals
#'  \item $significance_filter - information about the significance filter
#'  applied when remove_insignificant = TRUE, including the method used,
#'  the criterion, and how many cells were removed and retained
#'}
#'
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
#' data(LJ_daily_precipitation)
#'
#' # 1 Basic example using the partial correlation coefficient
#' example_basic <- daily_response_seascorr(response = data_MVA,
#'                           env_data_primary = LJ_daily_temperatures,
#'                           env_data_control = LJ_daily_precipitation,
#'                           row_names_subset = TRUE,
#'                           fixed_width = 25,
#'                           lower_limit = 35, upper_limit = 45,
#'                           remove_insignificant = FALSE,
#'                           aggregate_function_env_data_primary = 'mean',
#'                           aggregate_function_env_data_control = 'mean',
#'                           tidy_env_data_primary = FALSE,
#'                           tidy_env_data_control = TRUE,
#'                           alpha = 0.05, pcor_method = "spearman",
#'                           previous_year = FALSE, boot = TRUE,
#'                           boot_n = 10,
#'                           reference_window = "end", k = 5,
#'                           dc_method = "SLD",
#'                           day_interval = c(-100, 250),
#'                           skip_window_position = 50,
#'                           skip_window_length= 50
#'                           )
#' # summary(example_basic)
#' # plot(example_basic, type = 1)
#' # plot(example_basic, type = 2)
#' # example_basic$optimized_return
#' # example_basic$optimized_return_all
#' # example_basic$temporal_stability
#'
#' # 2 Example with fixed temporal time window
#' example_fixed_width <- daily_response_seascorr(response = data_MVA,
#'                           env_data_primary = LJ_daily_temperatures,
#'                           env_data_control = LJ_daily_precipitation,
#'                           row_names_subset = TRUE,
#'                           remove_insignificant = TRUE,
#'                           aggregate_function_env_data_primary = 'mean',
#'                           aggregate_function_env_data_control = 'mean',
#'                           alpha = 0.05,
#'                           dc_method = "SLD",
#'                           fixed_width = 45,
#'                           tidy_env_data_primary = FALSE,
#'                           tidy_env_data_control = TRUE,
#'                           reference_window = "end",
#'                           skip_window_position = 50,
#'                           skip_window_length= 50)
#'
#' # summary(example_fixed_width)
#' # plot(example_fixed_width, type = 1)
#' # plot(example_fixed_width, type = 2)
#' # example_fixed_width$optimized_return
#' # example_fixed_width$optimized_return_all
#'
#'
#' # 3 Example using quantiles for primary and control variables
#' example_quantile <- daily_response_seascorr(
#'   response = data_MVA,
#'   env_data_primary = LJ_daily_temperatures,
#'   env_data_control = LJ_daily_precipitation,
#'   row_names_subset = TRUE,
#'   fixed_width = 30,
#'   previous_year = TRUE,
#'   aggregate_function_env_data_primary = "quantile",
#'   quantile_prob_env_data_primary = 0.95,
#'   aggregate_function_env_data_control = "quantile",
#'   quantile_prob_env_data_control = 0.95,
#'   tidy_env_data_primary = FALSE,
#'   tidy_env_data_control = TRUE,
#'   skip_window_position = 50
#' )
#'
#' # 4 Example using two previous years plus the current year
#' example_two_previous_years <- daily_response_seascorr(
#'   response = data_MVA,
#'   env_data_primary = LJ_daily_temperatures,
#'   env_data_control = LJ_daily_precipitation,
#'   row_names_subset = TRUE,
#'   fixed_width = 60,
#'   previous_year = TRUE,
#'   number_previous_years = 2,
#'   day_interval = c(-1, 250),
#'   tidy_env_data_primary = FALSE,
#'   tidy_env_data_control = TRUE,
#'   skip_window_position = 50
#' )
#' }

daily_response_seascorr <- function(response, env_data_primary, env_data_control,
                                    lower_limit = 30,
                                    upper_limit = 90, fixed_width = 0,
                                    previous_year = FALSE, number_previous_years = NULL,
                                    pcor_method = "pearson",
                                    remove_insignificant = TRUE,
                                    remove_insignificant_method = "alpha",
                                    perm_n = 999,
                                    cluster_threshold = "nominal",
                                    alpha = .05, row_names_subset = FALSE,
                                    aggregate_function_env_data_primary = 'mean',
                                    quantile_prob_env_data_primary = 0.5,
                                    aggregate_function_env_data_control = 'mean',
                                    quantile_prob_env_data_control = 0.5,
                                    temporal_stability_check = "sequential", k = 2,
                                    k_running_window = 30,
                                    subset_years = NULL,
                                    ylimits = NULL, seed = NULL, tidy_env_data_primary = FALSE,
                                    tidy_env_data_control = FALSE,
                                    reference_window = 'start',  boot = FALSE, boot_n = 1000,
                                    boot_ci_type = "norm", boot_conf_int = 0.95,
                                    day_interval = NULL,
                                    dc_method = NULL,
                                    pcor_na_use = "everything",
                                    skip_window_length = 1,
                                    skip_window_position = 1){

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

  offset_start <- day_interval[1]
  offset_end <- day_interval[2]

  if (offset_start == 0 || offset_end == 0) {
    stop("day_interval cannot contain 0. Use negative values for previous-year days and positive values for current-year days.")
  }

  # Negative values in day_interval would imply previous-year climate data.
  # However, previous_year is the master switch. If previous_year = FALSE,
  # we ignore the supplied negative day_interval and perform a current-year
  # analysis using the default current-year interval.
  if (isFALSE(previous_year) && (offset_start < 0 || offset_end < 0)) {

    warning(paste0("Negative values were supplied in day_interval, but ",
                   "previous_year = FALSE. The day_interval argument is ",
                   "ignored and the analysis is performed for the current ",
                   "year only using day_interval = c(1, 366)."))

    day_interval <- c(1, days_per_year)
    offset_start <- day_interval[1]
    offset_end <- day_interval[2]
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
  if (upper_limit > max_window){

    upper_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The upper_limit is outside your day_interval and",
                     " therefore reduced to the maximum allowed: ",max_window,"."))
    }
  }

  if (lower_limit > max_window){

    lower_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The lower_limit is outside your day_interval and",
                     " therefore reduced to the minimum allowed: ",max_window,"."))
    }
  }

  # Also correction for fixed-window approach
  if (fixed_width > max_window){

    stop(paste0("The selected fixed_width is outside your day_interval. ",
                "Decrease the fixed_width argument to at least: ",max_window,"."))
  }

  # offset_end is converted to the number of days that must be skipped at the
  # end of the full environmental matrix.
  offset_end <- year_block_width - offset_end

  ##############################################################################


  # Studentized intervals need a variance estimate for every bootstrap
  # replicate. The statistics computed here return a single value per
  # replicate, so boot.ci() cannot build a studentized interval and previously
  # failed in a way that produced an uninformative error.

  # Missing values are not handled automatically. Tell the user explicitly,
  # so that an all-NA surface is never a silent surprise.
  if ((anyNA(env_data_primary) || anyNA(env_data_control)) &&
      !(pcor_na_use %in% c("complete.obs", "na.or.complete",
                           "pairwise.complete.obs"))) {

    warning(paste0("env_data_primary or env_data_control contains missing ",
                   "values and pcor_na_use = '", pcor_na_use, "'. Any window ",
                   "overlapping a missing value will return NA. Set ",
                   "pcor_na_use = 'pairwise.complete.obs' to compute over the ",
                   "available years instead."))
  }

  allowed_boot_ci_types <- c("norm", "perc", "basic", "bca")

  if (!(boot_ci_type %in% allowed_boot_ci_types)) {
    stop(paste0("boot_ci_type is '", boot_ci_type, "'. Instead it should be ",
                "one of: ", paste(allowed_boot_ci_types, collapse = ", "),
                ". Studentized intervals are not available because they ",
                "require a variance estimate for each bootstrap replicate."))
  }
  allowed_significance_methods <- c("alpha", "boot", "fdr_BY", "fdr_BH",
                                    "maxstat", "sidak", "cluster")

  if (!(remove_insignificant_method %in% allowed_significance_methods)) {
    stop(paste0("remove_insignificant_method is '", remove_insignificant_method,
                "'. Instead it should be one of: ",
                paste(allowed_significance_methods, collapse = ", "), "."))
  }

  if (!is.logical(remove_insignificant) ||
      length(remove_insignificant) != 1 ||
      is.na(remove_insignificant)) {
    stop("remove_insignificant must be either TRUE or FALSE.")
  }

  # The bootstrap filter needs bootstrap intervals to exist.
  if (isTRUE(remove_insignificant) &&
      remove_insignificant_method == "boot" && isFALSE(boot)) {

    warning(paste0("remove_insignificant_method = 'boot' requires boot = TRUE. ",
                   "No significance filter is applied."))

    remove_insignificant <- FALSE
  }

  # The max-statistic permutation test rebuilds the partial correlation of a
  # shuffled response against every window. Because the correlation between the
  # primary and control windows does not change when only the response is
  # shuffled, it is computed once and reused, so each permutation costs two
  # vectorised correlations rather than a partial correlation per window.

  # Cluster-extent inference needs the same permutation machinery as the
  # max-statistic test, plus a threshold at which clusters are formed.
  if (isTRUE(remove_insignificant) &&
      remove_insignificant_method == "cluster") {

    if (!is.numeric(perm_n) || length(perm_n) != 1 || is.na(perm_n) ||
        perm_n < 99) {
      stop("perm_n must be a single number of at least 99.")
    }

    if (!(identical(cluster_threshold, "nominal") ||
          (is.numeric(cluster_threshold) && length(cluster_threshold) == 1 &&
           !is.na(cluster_threshold) && cluster_threshold > 0 &&
           cluster_threshold < 1))) {
      stop(paste0("cluster_threshold must be \"nominal\" or a single number ",
                  "between 0 and 1 giving the absolute partial correlation at ",
                  "which clusters are formed."))
    }
  }

  if (isTRUE(remove_insignificant) &&
      remove_insignificant_method == "maxstat") {

    if (!is.numeric(perm_n) || length(perm_n) != 1 || is.na(perm_n) ||
        perm_n < 99) {
      stop("perm_n must be a single number of at least 99.")
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (fixed_width != 0){
    lower_limit = 30
    upper_limit = 200
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
  allowed_aggregate_functions <- c("mean", "median", "sum", "min", "max", "quantile")

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

  # Internal helper for aggregating daily data across rows
  aggregate_daily_window <- function(x, aggregate_function, quantile_prob) {

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

  aggregate_daily_window_primary <- function(x) {
    aggregate_daily_window(x,
                           aggregate_function = aggregate_function_env_data_primary,
                           quantile_prob = quantile_prob_env_data_primary)
  }

  aggregate_daily_window_control <- function(x) {
    aggregate_daily_window(x,
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

  # If there is a column name samp.depth in response data frame, warning is given
  if ("samp.depth" %in% colnames(response)){

    samp.depth_index <- grep("samp.depth", colnames(response))
    response <- response[, -samp.depth_index,F]

    warning("Removed the samp.depth from response data frame")
  }

  # If there are more than 2 columns in response data frame, give a warning
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

    if (colnames_tidy_DF[2] != "DOY"){
      stop(paste("env_data_control was inserted in tidy version (tidy_env_data_control is set to TRUE).",
                 "The second column name of the env_data_control should be 'DOY', but it is",
                 colnames_tidy_DF[2], "instead!"))
    }

    value_variable = colnames(env_data_control)[3]
    env_data_control <- dcast(env_data_control, Year~DOY, value.var = value_variable)
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

    if (colnames_tidy_DF[2] != "DOY"){
      stop(paste("env_data_primary was inserted in tidy version (tidy_env_data_primary is set to TRUE).",
                 "The second column name of the env_data_primary should be 'DOY', but it is",
                 colnames_tidy_DF[2], "instead!"))
    }

    value_variable = colnames(env_data_primary)[3]
    env_data_primary <- dcast(env_data_primary, Year~DOY, value.var = value_variable)
    env_data_primary <- years_to_rownames(env_data_primary, "Year")
  }

  # PART 1 - general data arrangements, warnings and stops
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

  if (lower_limit > upper_limit)
    stop("lower_limit can not be higher than upper_limit!")

  if (lower_limit > year_block_width | lower_limit < 1)
    stop(paste0("lower_limit out of bounds! It should be between 1 and ",
                year_block_width))

  if (upper_limit > year_block_width | upper_limit < 1)
    stop(paste0("upper_limit out of bounds! It should be between 1 and ",
                year_block_width))

  # Make sure the selected method is appropriate
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

  # A) env_data_primary
  env_temp_primary <- env_data_primary[row.names(env_data_primary) %in% row.names(response),]

  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    subset_seq <- seq(lower_subset, upper_subset)
    env_temp_primary <- subset(env_temp_primary, row.names(env_temp_primary) %in% subset_seq)
  }

  na_problem <- data.frame(na_sum = rowSums(is.na(env_temp_primary)))
  na_problem <- na_problem[na_problem$na_sum > missing_day_threshold, , F]
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
  na_problem <- na_problem[na_problem$na_sum > missing_day_threshold, , F]
  problematic_years <- paste0(row.names(na_problem), sep = "", collapse=", ")

  if (nrow(na_problem) > 0){

    warning(paste0("Problematic years with missing values are present in env_data_control: ", problematic_years))

  }

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

  # if row.names of env_data_primary and the response data frames are not equal - warning is given.
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
                  " in the env_data_control data frame. Change the subset_years argument"))
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

  # ---- shared window-coverage gate (method-agnostic) -----------------
  # A window is usable only if it has at least `min_complete` complete
  # response/predictor pairs. Used by cor, lm and brnn so all methods
  # return NA (not a silent partial fit) on fully/mostly missing windows.
  window_is_usable_pcor <- function(df, min_n = 4L) {
    sum(stats::complete.cases(df[, c("x", "y", "z")])) >= min_n
  }

  # PART 2 - the calculation of day-wise correlations
  # A) fixed window approach
  # B) Variable window approach

  # these are lists for climate and holder for saving mm1 and mm2
  list_climate_primary <- list()
  list_climate_control <- list()

  # Vectorised bootstrap of the partial correlation. The resample indices are
  # drawn once and shared by every window; each replicate is rebuilt from the
  # three pairwise correlations rather than by resampling. "bca" needs
  # jackknife estimates and keeps the boot() route.
  boot_fast_obj <- NULL
  if (boot == TRUE && boot_ci_type %in% c("norm", "perc", "basic")) {

    boot_fast_obj <- boot_fast_prep(y = response[, 1],
                                    boot_n = boot_n,
                                    cor_method = pcor_method,
                                    cor_na_use = pcor_na_use)
  }

  mm1 <- 1
  mm2 <- 1

  # A) The fixed_window approach
  if (fixed_width != 0) {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be stored in this matrix
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

        pb <- txtProgressBar(min = 0, max = ceiling(((ncol(env_data_primary) - fixed_width - offset_end - offset_start + 1))/(skip_window_position)), style = 3)

      }

    }

    b = 0

    # An iterating loop. In each iteration x is calculated and represents
    # response (dependent) variable. X is a moving average. Window width of
    # a moving window is fixed_width. Next, partial correlation is stored in
    # temporal matrix.
    for (j in (seq((0 + offset_start -1), (ncol(env_data_primary) - max((fixed_width + offset_end), offset_end)), by = skip_window_position))) {

      b = b + 1

      x1 <- aggregate_daily_window_primary(
        env_data_primary[1:nrow(env_data_primary),
                         (1 + j):(j + fixed_width),
                         drop = FALSE]
      )

      x2 <- aggregate_daily_window_control(
        env_data_control[1:nrow(env_data_control),
                         (1 + j):(j + fixed_width),
                         drop = FALSE]
      )

      if (!is.null(dc_method)){


        if (dc_method == "SLD"){

          tmp_model <- lm(x1 ~ seq(1:length(x1)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x1 - tmp_pred

          x1 <- data.frame(x1 = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x1 <- matrix(x1, nrow = nrow(env_data_primary), ncol = 1)

      }


      if (!is.null(dc_method)){

        if (dc_method == "SLD"){

          tmp_model <- lm(x2 ~ seq(1:length(x2)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x2 - tmp_pred

          x2 <- data.frame(x2 = tmp_res/sd(tmp_res, na.rm = TRUE))

        }

      } else {

        x2 <- matrix(x2, nrow = nrow(env_data_primary), ncol = 1)

      }

      my_temporal_data <- cbind(response[, 1], x1, x2)
      colnames(my_temporal_data) <- c("x", "y", "z")

      win_ok <- window_is_usable_pcor(my_temporal_data)

      # A window that fails the usability gate is NA in temporal_matrix, so it
      # must not enter the aggregated matrices either. The permutation nulls
      # used by "maxstat" and "cluster" are built from these columns and would
      # otherwise be taken over a larger window set than the observed surface.
      if (win_ok) {
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
      }

      if (!win_ok) {
        temporal_correlation <- NA
      } else if (boot == FALSE){

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

        if (!is.null(boot_fast_obj)) {

          boot_t <- boot_fast_pcor_replicates(boot_fast_obj,
                                              my_temporal_data[, "y"],
                                              my_temporal_data[, "z"])

          boot_t0 <- as.numeric(partial.r(data = my_temporal_data,
                                          x = c("x", "y"), y = "z",
                                          use = pcor_na_use,
                                          method = pcor_method)[2])

          temporal_correlation <- mean(boot_t)

          boot_ci_vals <- boot_fast_ci(boot_t, boot_t0,
                                       conf = boot_conf_int,
                                       type = boot_ci_type)
          temporal_lower <- boot_ci_vals[1]
          temporal_upper <- boot_ci_vals[2]
          ci_int <- NULL

        } else {

          calc <- boot(data = my_temporal_data, statistic = boot_f_pcor, R = boot_n,
                       cor.type = pcor_method, pcor.use = pcor_na_use)

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

    # An iterating double loop: 1 outer loop) iterating from lower_limit :
    # upper_limit defines windo.width used for a moving window. 2) inner loop
    # defines the starting position of a moving window.
    # In each iteration, x is calculated and represents a response (dependent)
    # variable. x is a moving average, based on rowMeans/apply function.
    # Next, statistical metric is calculated based on a selected method (partial).
    # Calculation is stored in temporal matrix in a proper place.
    # The position of stored calculation is informative later used for
    # indicating optimal values.

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

      for (j in seq((0 + offset_start -1), (ncol(env_data_primary) - max((K + offset_end), offset_end)), by = skip_window_position)) {

        x1 <- aggregate_daily_window_primary(
          env_data_primary[1:nrow(env_data_primary),
                           (1 + j):(j + K),
                           drop = FALSE]
        )

        x2 <- aggregate_daily_window_control(
          env_data_control[1:nrow(env_data_control),
                           (1 + j):(j + K),
                           drop = FALSE]
        )

        if (!is.null(dc_method)){


          if (dc_method == "SLD"){

            tmp_model <- lm(x1 ~ seq(1:length(x1)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x1 - tmp_pred

            x1 <- data.frame(x1 = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x1 <- matrix(x1, nrow = nrow(env_data_primary), ncol = 1)

        }


        if (!is.null(dc_method)){


          if (dc_method == "SLD"){

            tmp_model <- lm(x2 ~ seq(1:length(x2)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x2 - tmp_pred

            x2 <- data.frame(x2 = tmp_res/sd(tmp_res, na.rm = TRUE))

          }

        } else {

          x2 <- matrix(x2, nrow = nrow(env_data_primary), ncol = 1)

        }


        my_temporal_data <- cbind(response[, 1], x1, x2)
        colnames(my_temporal_data) <- c("x", "y", "z")

        win_ok <- window_is_usable_pcor(my_temporal_data)

        # A window that fails the usability gate is NA in temporal_matrix, so it
        # must not enter the aggregated matrices either. The permutation nulls
        # used by "maxstat" and "cluster" are built from these columns and would
        # otherwise be taken over a larger window set than the observed surface.
        if (win_ok) {
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
        }

        if (!win_ok) {

          temporal_correlation <- NA

        } else if (boot == FALSE){

          temporal_correlation <- partial.r(data = my_temporal_data, x=c("x","y"), y="z",
                                            use=pcor_na_use,method = pcor_method)[2]

          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1] <- as.numeric(temporal_correlation)[1]
          } else if (reference_window == 'end'){
            temporal_matrix[(K - lower_limit) + 1, j + K] <- as.numeric(temporal_correlation)[1]
          } else if (reference_window == 'middle'){
            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- as.numeric(temporal_correlation)[1]
          }

        } else if (boot == TRUE){

          if (!is.null(boot_fast_obj)) {

            boot_t <- boot_fast_pcor_replicates(boot_fast_obj,
                                                my_temporal_data[, "y"],
                                                my_temporal_data[, "z"])

            boot_t0 <- as.numeric(partial.r(data = my_temporal_data,
                                            x = c("x", "y"), y = "z",
                                            use = pcor_na_use,
                                            method = pcor_method)[2])

            temporal_correlation <- mean(boot_t)

            boot_ci_vals <- boot_fast_ci(boot_t, boot_t0,
                                         conf = boot_conf_int,
                                         type = boot_ci_type)
            temporal_lower <- boot_ci_vals[1]
            temporal_upper <- boot_ci_vals[2]
            ci_int <- NULL

          } else {

            calc <- boot(data = my_temporal_data, statistic = boot_f_pcor, R = boot_n,
                         cor.type = pcor_method, pcor.use = pcor_na_use)

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

  # To enhance the visualisation, insignificant values are removed if
  # remove_insignificant == TRUE. The criterion is set by
  # remove_insignificant_method:
  #   "alpha"  per-window significance threshold; does not account for how
  #            many windows were screened
  #   "boot"   bootstrap confidence interval excludes zero (per-window)
  #   "fdr_BY" false discovery rate, Benjamini-Yekutieli, valid under
  #            arbitrary dependence between windows
  #   "fdr_BH" false discovery rate, Benjamini-Hochberg, valid under positive
  #            dependence and less conservative than BY
  significance_filter <- list(
    applied = FALSE,
    method = NA_character_,
    criterion = NA_character_,
    removed_cells = 0L,
    retained_cells = sum(!is.na(temporal_matrix)),
    total_cells = length(temporal_matrix)
  )

  if (isTRUE(remove_insignificant)) {

    # Partial correlation controls for one variable, therefore the degrees of
    # freedom are n - 3 rather than n - 2.
    n_obs <- nrow(response)
    df_pcor <- n_obs - 3
    original_non_missing <- !is.na(temporal_matrix)
    keep <- NULL

    # 1 Plain alpha threshold
    if (remove_insignificant_method == "alpha") {

      critical_threshold_cor <- critical_r(n_obs, alpha = alpha)

      keep <- abs(temporal_matrix) >= abs(critical_threshold_cor)

      significance_filter$applied <- TRUE
      significance_filter$method <- "alpha"
      significance_filter$criterion <- paste0("|partial r| >= ",
                                              round(critical_threshold_cor, 3),
                                              " (alpha = ", alpha, ")")

      # 2 Bootstrap confidence interval
    } else if (remove_insignificant_method == "boot") {

      if (is.null(temporal_matrix_lower) ||
          is.null(temporal_matrix_upper) ||
          !all(dim(temporal_matrix) == dim(temporal_matrix_lower)) ||
          !all(dim(temporal_matrix) == dim(temporal_matrix_upper))) {

        warning(paste0("Bootstrap confidence-interval filtering was requested, ",
                       "but boot_lower and boot_upper are not available with the ",
                       "same dimensions as the calculation matrix. The bootstrap ",
                       "filter is ignored."))

      } else {

        keep <- (temporal_matrix_lower > 0 & temporal_matrix_upper > 0) |
          (temporal_matrix_lower < 0 & temporal_matrix_upper < 0)

        significance_filter$applied <- TRUE
        significance_filter$method <- "boot"
        significance_filter$criterion <- "Bootstrap confidence interval excludes zero"
      }

      # 3 False discovery rate
    } else if (remove_insignificant_method == "fdr_BY" |
               remove_insignificant_method == "fdr_BH") {

      fdr_type <- sub("^fdr_", "", remove_insignificant_method)

      t_surface <- temporal_matrix * sqrt(df_pcor / (1 - temporal_matrix ^ 2))
      p_surface <- 2 * stats::pt(-abs(t_surface), df = df_pcor)

      p_adjusted <- stats::p.adjust(as.vector(p_surface), method = fdr_type)
      p_adjusted <- matrix(p_adjusted, nrow = nrow(p_surface))

      keep <- p_adjusted < alpha

      significance_filter$applied <- TRUE
      significance_filter$method <- remove_insignificant_method
      significance_filter$criterion <- paste0("Adjusted p < ", alpha,
                                              " (", fdr_type, ")")

      # 4 Max-statistic permutation
    } else if (remove_insignificant_method == "maxstat") {

      climate_primary <- do.call(cbind, list_climate_primary)
      climate_control <- do.call(cbind, list_climate_control)

      if (is.null(climate_primary) || is.null(climate_control) ||
          ncol(climate_primary) < 1 ||
          ncol(climate_primary) != ncol(climate_control)) {

        warning(paste0("remove_insignificant_method = 'maxstat' was requested, ",
                       "but the matrices of aggregated primary and control ",
                       "windows are not available or do not match. The ",
                       "significance filter is ignored."))

      } else {

        shared_years <- intersect(row.names(response), row.names(climate_primary))

        if (length(shared_years) < 10) {

          warning(paste0("remove_insignificant_method = 'maxstat' was requested, ",
                         "but fewer than 10 years are shared between the response ",
                         "and the aggregated climate windows. The significance ",
                         "filter is ignored."))

        } else {

          y_perm  <- response[shared_years, 1]
          x1_perm <- climate_primary[shared_years, , drop = FALSE]
          x2_perm <- climate_control[shared_years, , drop = FALSE]

          # Correlation between each primary window and its paired control
          # window. This does not involve the response, so it is unchanged by
          # permutation and is computed only once.
          r_yz <- vapply(seq_len(ncol(x1_perm)), function(win_i) {
            suppressWarnings(stats::cor(x1_perm[, win_i], x2_perm[, win_i],
                                        method = pcor_method, use = pcor_na_use))
          }, numeric(1))

          denom_z <- sqrt(1 - r_yz ^ 2)

          # Null distribution of the largest absolute partial correlation
          # obtained anywhere in the search when the response is unrelated to
          # climate. The same shuffle is applied to both correlations so that
          # the response stays internally consistent within a permutation.
          null_max <- rep(NA_real_, perm_n)

          for (perm_i in seq_len(perm_n)) {

            y_shuffled <- sample(y_perm)

            r_xy <- suppressWarnings(
              as.numeric(stats::cor(y_shuffled, x1_perm,
                                    method = pcor_method, use = pcor_na_use)))

            r_xz <- suppressWarnings(
              as.numeric(stats::cor(y_shuffled, x2_perm,
                                    method = pcor_method, use = pcor_na_use)))

            perm_pcor <- (r_xy - r_xz * r_yz) / (sqrt(1 - r_xz ^ 2) * denom_z)

            null_max[perm_i] <- max(abs(perm_pcor), na.rm = TRUE)
          }

          # Exact permutation threshold: the k-th largest null value with
          # k = floor(alpha * (B + 1)), B being the number of valid
          # permutations. This is the same construction already used for
          # p_global below. stats::quantile() interpolates and lands just
          # under this value, which makes the test marginally liberal.
          null_sorted <- sort(null_max, decreasing = TRUE)
          critical_threshold_perm <- as.numeric(
            null_sorted[max(1L, floor(alpha * (length(null_sorted) + 1)))])

          observed_max <- max(abs(temporal_matrix), na.rm = TRUE)

          keep <- abs(temporal_matrix) >= critical_threshold_perm

          p_global <- (sum(null_max >= observed_max, na.rm = TRUE) + 1) /
            (perm_n + 1)

          significance_filter$applied <- TRUE
          significance_filter$method <- "maxstat"
          significance_filter$criterion <- paste0("|partial r| >= ",
                                                  round(critical_threshold_perm, 3),
                                                  " (max-statistic permutation, ",
                                                  perm_n, " permutations, alpha = ",
                                                  alpha, ")")
          significance_filter$critical_r <- critical_threshold_perm
          significance_filter$observed_max_r <- observed_max
          significance_filter$p_global <- p_global
          significance_filter$perm_n <- perm_n
          significance_filter$n_windows_searched <- ncol(x1_perm)
        }
      }

      # 5 Sidak correction on the effective number of tests
    } else if (remove_insignificant_method == "sidak") {

      climate_primary <- do.call(cbind, list_climate_primary)
      climate_control <- do.call(cbind, list_climate_control)

      if (is.null(climate_primary) || is.null(climate_control) ||
          ncol(climate_primary) < 2 ||
          ncol(climate_primary) != ncol(climate_control)) {

        warning(paste0("remove_insignificant_method = 'sidak' was requested, ",
                       "but the matrices of aggregated primary and control ",
                       "windows are not available or do not match. The ",
                       "significance filter is ignored."))

      } else {

        shared_years <- intersect(row.names(response), row.names(climate_primary))

        if (length(shared_years) < 10) {

          warning(paste0("remove_insignificant_method = 'sidak' was requested, ",
                         "but fewer than 10 years are shared between the ",
                         "response and the aggregated climate windows. The ",
                         "significance filter is ignored."))

        } else {

          n_obs_sidak <- length(shared_years)

          x1_sidak <- climate_primary[shared_years, , drop = FALSE]
          x2_sidak <- climate_control[shared_years, , drop = FALSE]

          complete_win <- colSums(is.na(x1_sidak)) == 0 &
            colSums(is.na(x2_sidak)) == 0
          x1_sidak <- x1_sidak[, complete_win, drop = FALSE]
          x2_sidak <- x2_sidak[, complete_win, drop = FALSE]

          if (ncol(x1_sidak) < 2) {

            warning(paste0("remove_insignificant_method = 'sidak' was ",
                           "requested, but fewer than two complete window ",
                           "pairs remain. The significance filter is ",
                           "ignored."))

          } else {

            x1_sidak <- x1_sidak[, apply(x1_sidak, 2, stats::sd) > 0,
                                 drop = FALSE]

            if (ncol(x1_sidak) < 2) {

              warning(paste0("remove_insignificant_method = 'sidak' was ",
                             "requested, but fewer than two primary windows ",
                             "vary. The significance filter is ignored."))

            } else {

              # Effective number of independent tests among the screened
              # primary windows, following Li and Ji (2005). The correlation
              # matrix among windows has rank at most the number of years, so
              # its non-zero eigenvalues come from the years-by-years matrix
              # rather than the windows-by-windows one.
              z_sidak <- scale(x1_sidak)
              gram <- tcrossprod(z_sidak) / (n_obs_sidak - 1)

              eig <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
              eig <- eig[eig > 1e-10]

              m_eff <- sum((eig >= 1) + (eig - floor(eig)))

              # One degree of freedom is spent on the control variable
              alpha_sidak <- 1 - (1 - alpha) ^ (1 / m_eff)
              t_sidak <- stats::qt(1 - alpha_sidak / 2, df = n_obs_sidak - 3)
              critical_threshold_sidak <- t_sidak /
                sqrt(t_sidak ^ 2 + (n_obs_sidak - 3))

              keep <- abs(temporal_matrix) >= critical_threshold_sidak

              significance_filter$applied <- TRUE
              significance_filter$method <- "sidak"
              significance_filter$criterion <- paste0("|partial r| >= ",
                                                      round(critical_threshold_sidak, 3),
                                                      " (Sidak on ", round(m_eff, 1),
                                                      " effective tests, alpha = ",
                                                      alpha, ")")
              significance_filter$critical_r <- critical_threshold_sidak
              significance_filter$effective_tests <- m_eff
              significance_filter$alpha_per_test <- alpha_sidak
              significance_filter$n_windows_searched <- ncol(x1_sidak)
            }
          }
        }
      }

      # 6 Cluster-extent inference
    } else if (remove_insignificant_method == "cluster") {

      climate_primary <- do.call(cbind, list_climate_primary)
      climate_control <- do.call(cbind, list_climate_control)

      if (is.null(climate_primary) || is.null(climate_control) ||
          ncol(climate_primary) < 1 ||
          ncol(climate_primary) != ncol(climate_control)) {

        warning(paste0("remove_insignificant_method = 'cluster' was requested, ",
                       "but the matrices of aggregated primary and control ",
                       "windows are not available or do not match. The ",
                       "significance filter is ignored."))

      } else {

        shared_years <- intersect(row.names(response), row.names(climate_primary))

        if (length(shared_years) < 10) {

          warning(paste0("remove_insignificant_method = 'cluster' was ",
                         "requested, but fewer than 10 years are shared ",
                         "between the response and the aggregated climate ",
                         "windows. The significance filter is ignored."))

        } else {

          y_perm <- response[shared_years, 1]
          x1_perm <- climate_primary[shared_years, , drop = FALSE]
          x2_perm <- climate_control[shared_years, , drop = FALSE]
          n_obs_cl <- length(shared_years)

          # Cluster inference needs the permuted statistics laid out on the
          # same surface as the observed ones, unlike the max-statistic test
          # which needs only their maximum. Each aggregated window is named
          # "start_end", which gives its width and position; the mapping is
          # checked against the observed surface before it is used.
          win_nm <- colnames(x1_perm)
          parts <- strsplit(win_nm, "_", fixed = TRUE)
          w_start <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 1)))
          w_end <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
          w_width <- w_end - w_start + 1

          surf_widths <- suppressWarnings(as.numeric(row.names(temporal_matrix)))

          if (nrow(temporal_matrix) == 1) {
            row_idx <- rep(1L, length(w_width))
          } else {
            row_idx <- match(w_width, surf_widths)
          }

          col_idx <- if (reference_window == "start") {
            w_start
          } else if (reference_window == "end") {
            w_end
          } else {
            round2(w_start + w_width / 2, 0)
          }

          usable <- !is.na(row_idx) & !is.na(col_idx) &
            row_idx >= 1 & row_idx <= nrow(temporal_matrix) &
            col_idx >= 1 & col_idx <= ncol(temporal_matrix)

          map_ok <- FALSE
          if (any(usable)) {

            r_xy_o <- suppressWarnings(
              as.numeric(stats::cor(y_perm, x1_perm[, usable, drop = FALSE],
                                    method = pcor_method, use = pcor_na_use)))
            r_xz_o <- suppressWarnings(
              as.numeric(stats::cor(y_perm, x2_perm[, usable, drop = FALSE],
                                    method = pcor_method, use = pcor_na_use)))
            r_yz_o <- vapply(which(usable), function(win_i) {
              suppressWarnings(stats::cor(x1_perm[, win_i], x2_perm[, win_i],
                                          method = pcor_method,
                                          use = pcor_na_use))
            }, numeric(1))

            obs_pcor <- (r_xy_o - r_xz_o * r_yz_o) /
              (sqrt(1 - r_xz_o ^ 2) * sqrt(1 - r_yz_o ^ 2))

            on_surface <- temporal_matrix[cbind(row_idx[usable],
                                                col_idx[usable])]
            both <- is.finite(obs_pcor) & is.finite(on_surface)
            map_ok <- sum(both) > 0 &&
              max(abs(obs_pcor[both] - on_surface[both])) < 1e-8
          }

          if (!map_ok) {

            warning(paste0("remove_insignificant_method = 'cluster' was ",
                           "requested, but the aggregated windows could not be ",
                           "matched to their positions on the response ",
                           "surface. The significance filter is ignored."))

          } else {

            idx_cells <- cbind(row_idx[usable], col_idx[usable])
            x1_use <- x1_perm[, usable, drop = FALSE]
            x2_use <- x2_perm[, usable, drop = FALSE]

            form_thr <- if (identical(cluster_threshold, "nominal")) {
              critical_r(n_obs_cl, alpha = alpha)
            } else {
              cluster_threshold
            }

            # the primary-control correlation does not involve the response,
            # so it is unchanged by permutation and computed only once
            r_yz <- vapply(seq_len(ncol(x1_use)), function(win_i) {
              suppressWarnings(stats::cor(x1_use[, win_i], x2_use[, win_i],
                                          method = pcor_method,
                                          use = pcor_na_use))
            }, numeric(1))
            denom_z <- sqrt(1 - r_yz ^ 2)

            null_clust <- rep(NA_real_, perm_n)
            null_max <- rep(NA_real_, perm_n)

            perm_surface <- matrix(NA_real_, nrow(temporal_matrix),
                                   ncol(temporal_matrix))

            for (perm_i in seq_len(perm_n)) {

              y_shuffled <- sample(y_perm)

              r_xy <- suppressWarnings(
                as.numeric(stats::cor(y_shuffled, x1_use,
                                      method = pcor_method, use = pcor_na_use)))
              r_xz <- suppressWarnings(
                as.numeric(stats::cor(y_shuffled, x2_use,
                                      method = pcor_method, use = pcor_na_use)))

              perm_pcor <- (r_xy - r_xz * r_yz) /
                (sqrt(1 - r_xz ^ 2) * denom_z)

              perm_surface[] <- NA_real_
              perm_surface[idx_cells] <- perm_pcor

              null_max[perm_i] <- max(abs(perm_pcor), na.rm = TRUE)
              null_clust[perm_i] <- max(cluster_fast_masses(perm_surface,
                                                            form_thr))
            }

            # Exact permutation thresholds, as in the "maxstat" branch above.
            null_clust_sorted <- sort(null_clust, decreasing = TRUE)
            crit_mass <- as.numeric(null_clust_sorted[
              max(1L, floor(alpha * (length(null_clust_sorted) + 1)))])

            null_max_sorted <- sort(null_max, decreasing = TRUE)
            crit_max_cl <- as.numeric(null_max_sorted[
              max(1L, floor(alpha * (length(null_max_sorted) + 1)))])

            obs_mask <- !is.na(temporal_matrix) &
              abs(temporal_matrix) > form_thr
            cc <- cluster_fast_label(obs_mask)

            keep <- matrix(FALSE, nrow(temporal_matrix),
                           ncol(temporal_matrix))
            n_retained <- 0L
            obs_masses <- numeric(0)

            if (cc$n > 0L) {
              for (cl in seq_len(cc$n)) {
                in_cl <- cc$labels == cl
                mass <- sum(abs(temporal_matrix[in_cl]), na.rm = TRUE)
                obs_masses <- c(obs_masses, mass)
                if (mass > crit_mass) {
                  keep[in_cl] <- TRUE
                  n_retained <- n_retained + 1L
                }
              }
            }

            significance_filter$applied <- TRUE
            significance_filter$method <- "cluster"
            significance_filter$criterion <- paste0("cluster mass > ",
                                                    round(crit_mass, 1),
                                                    ", clusters formed at ",
                                                    "|partial r| > ",
                                                    round(form_thr, 3), " (",
                                                    perm_n,
                                                    " permutations, alpha = ",
                                                    alpha, ")")
            significance_filter$forming_threshold <- form_thr
            significance_filter$cluster_mass_critical <- crit_mass
            significance_filter$n_clusters_formed <- cc$n
            significance_filter$n_clusters_retained <- n_retained
            significance_filter$largest_cluster_mass <-
              if (length(obs_masses)) max(obs_masses) else 0
            significance_filter$maxstat_critical_r <- crit_max_cl
            significance_filter$perm_n <- perm_n
            significance_filter$n_windows_searched <- ncol(x1_use)
          }
        }
      }
    }

    # Apply whichever mask was produced
    if (!is.null(keep)) {

      keep[is.na(keep)] <- FALSE

      temporal_matrix[!keep] <- NA

      # The same mask is applied to the confidence-interval matrices, so that
      # a filtered surface is not accompanied by intervals for windows the
      # filter removed.
      if (exists("temporal_matrix_lower") && !is.null(temporal_matrix_lower) &&
          all(dim(temporal_matrix_lower) == dim(temporal_matrix))) {
        temporal_matrix_lower[!keep] <- NA
        temporal_matrix_upper[!keep] <- NA
      }

      significance_filter$removed_cells <- sum(original_non_missing & !keep,
                                               na.rm = TRUE)
      significance_filter$retained_cells <- sum(!is.na(temporal_matrix))
      significance_filter$total_cells <- length(temporal_matrix)
    }
  }



  if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){

    if (isTRUE(significance_filter$applied) &&
        identical(significance_filter$method, "maxstat")) {

      warning(paste0("No window exceeded the max-statistic permutation ",
                     "threshold (global p = ",
                     signif(significance_filter$p_global, 3),
                     "). This indicates that the strongest partial correlation ",
                     "found anywhere in the search is no larger than expected ",
                     "when the response is unrelated to climate."))

    } else if (isTRUE(significance_filter$applied)) {
      warning(paste0("All calculations were removed by the significance filter (",
                     significance_filter$method, "). Change alpha or ",
                     "remove_insignificant_method."))
    } else {
      warning("All calculations are insignificant! Please change the alpha argument.")
    }

    final_list <- list(calculations = temporal_matrix,
                       method = "pcor",
                       metric = pcor_method, analysed_period = NA,
                       optimized_return = NA,
                       optimized_return_all = NA,
                       transfer_function = NA, temporal_stability = NA,
                       cross_validation = NA,
                       plot_heatmap = NA,
                       plot_extreme = NA,
                       type = "daily",
                       reference_window = reference_window,
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate_primary = NA,
                       aggregated_climate_control = NA,
                       significance_filter = significance_filter)

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
      if (isTRUE(significance_filter$applied) &&
          identical(significance_filter$method, "maxstat")) {

        stop(paste0("No window exceeded the max-statistic permutation ",
                    "threshold (global p = ",
                    signif(significance_filter$p_global, 3),
                    "). This indicates that the strongest partial correlation ",
                    "found anywhere in the search is no larger than expected ",
                    "when the response is unrelated to climate."))

      } else if (isTRUE(significance_filter$applied)) {
        stop(paste0("All calculations were removed by the significance filter (",
                    significance_filter$method, "). Change alpha or ",
                    "remove_insignificant_method!"))
      } else {
        stop("All calculations are insignificant! Change the alpha argument!")
      }
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

    # The fourth return element is being created: rowMeans/ apply of optimal sequence:
    # So, here we consider more options, based on the reference_winow
    # 1. reference window = "start"
    if (reference_window == 'start'){

      x1 <- data.frame(aggregate_daily_window_primary(
        env_data_primary[, as.numeric(plot_column):
                           (as.numeric(plot_column) + as.numeric(row_index) - 1),
                         drop = FALSE]
      ))

      x2 <- data.frame(aggregate_daily_window_control(
        env_data_control[, as.numeric(plot_column):
                           (as.numeric(plot_column) + as.numeric(row_index) - 1),
                         drop = FALSE]
      ))

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      x1_original <- data.frame(aggregate_daily_window_primary(
        env_data_primary_original[, as.numeric(plot_column):
                                    (as.numeric(plot_column) + as.numeric(row_index) - 1),
                                  drop = FALSE]
      ))

      x2_original <- data.frame(aggregate_daily_window_control(
        env_data_control_original[, as.numeric(plot_column):
                                    (as.numeric(plot_column) + as.numeric(row_index) - 1),
                                  drop = FALSE]
      ))

    }


    # Option 2, reference window = "end"
    if (reference_window == 'end'){

      x1 <- data.frame(aggregate_daily_window_primary(
        env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                           as.numeric(plot_column),
                         drop = FALSE]
      ))

      x2 <- data.frame(aggregate_daily_window_control(
        env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                           as.numeric(plot_column),
                         drop = FALSE]
      ))

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      x1_original <- data.frame(aggregate_daily_window_primary(
        env_data_primary_original[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                    as.numeric(plot_column),
                                  drop = FALSE]
      ))

      x2_original <- data.frame(aggregate_daily_window_control(
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

      x1 <- data.frame(aggregate_daily_window_primary(
        env_data_primary[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                           (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                         drop = FALSE]
      ))

      x2 <- data.frame(aggregate_daily_window_control(
        env_data_control[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                           (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                         drop = FALSE]
      ))

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      x1_original <- data.frame(aggregate_daily_window_primary(
        env_data_primary_original[, (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
                                    (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2),
                                  drop = FALSE]
      ))

      x2_original <- data.frame(aggregate_daily_window_control(
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
        tmp_res <- x1 - tmp_pred
        x1 <- data.frame(x1 = tmp_res/sd(tmp_res, na.rm = TRUE))

        x2 <- x2[,1]
        tmp_model <- lm(x2 ~ seq(1:length(x2)))
        tmp_pred <- predict(tmp_model)
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
        tmp_res <- x1_original - tmp_pred
        x1_original <- data.frame(x1_original = tmp_res/sd(tmp_res, na.rm = TRUE))

        x2_original <- x2_original[,1]
        tmp_model <- lm(x2_original ~ seq(1:length(x2_original)))
        tmp_pred <- predict(tmp_model)
        tmp_res <- x2_original - tmp_pred
        x2_original <- data.frame(x2_original = tmp_res/sd(tmp_res, na.rm = TRUE))

      }

    }

    # x1_full_original <- cbind(x1_original, x2_original)
    x1_full_original <- merge(x1_original, x2_original, by = 0, all = TRUE)
    colnames(x1_full_original)[1] <- "year"
    x1_full_original <- years_to_rownames(x1_full_original, "year")
    colnames(x1_full_original)[ncol(x1_full_original)-1] <- "Optimized_return_primary"
    colnames(x1_full_original)[ncol(x1_full_original)] <- "Optimized_return_control"

    # x1_full_original <- x1_original
    # colnames(x1_full_original) <- "Optimized_return"

    colnames(x1) <- "Optimized.rowNames"
    my_temporal_data <- cbind(x1_full[,1], x1_full[,2], x1_full[,3])
    colnames(my_temporal_data) <- c("x", "y", "z")
    test_calculation <- partial.r(data=my_temporal_data, x=c("x","y"), y="z", use=pcor_na_use,method = pcor_method)[2]
    test_logical <- as.numeric(max_calculation) == as.numeric(test_calculation)

    # if (test_logical == FALSE){
    #  stop(paste0("Test calculation and max calculation do not match.",
    #       "Please contact the maintainer of the dendroTools R package"))
    #}


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

    # dataset = data.frame(optimized_return =x1[,1], proxy = response)

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

    final_list <- list(calculations = temporal_matrix, method = "pcor",
                       metric = pcor_method, analysed_period = analysed_period,
                       optimized_return = x1_full,
                       optimized_return_all = x1_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = NA,
                       type = "daily",
                       reference_window = reference_window,
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate_primary = do.call(cbind, list_climate_primary),
                       aggregated_climate_control = do.call(cbind, list_climate_control),
                       significance_filter = significance_filter)


    plot_heatmapA <- plot_heatmap(final_list, reference_window = reference_window, type = "daily")
    plot_extremeA <- plot_extreme(final_list, ylimits = ylimits, reference_window = reference_window, type = "daily")


    final_list <- list(calculations = temporal_matrix, method = "pcor",
                       metric = pcor_method, analysed_period = analysed_period,
                       optimized_return = x1_full,
                       optimized_return_all = x1_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = NA,
                       plot_heatmap = plot_heatmapA,
                       plot_extreme = plot_extremeA,
                       type = "daily",
                       reference_window = reference_window,
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate_primary = do.call(cbind, list_climate_primary),
                       aggregated_climate_control = do.call(cbind, list_climate_control),
                       significance_filter = significance_filter)

    class(final_list) <- "dmrs"

    return(final_list)
  }
}
