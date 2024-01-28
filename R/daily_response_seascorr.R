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
#' @param previous_year if set to TRUE, env_data_primary, env_data_control  and
#' response variables will be rearranged in a way, that also previous year will
#' be used for calculations of selected statistical metric.
#' @param remove_insignificant if set to TRUE, removes all correlations bellow
#' the significant threshold level, based on a selected alpha.
#' @param alpha significance level used to remove insignificant calculations.
#' @param row_names_subset if set to TRUE, row.names are used to subset
#' env_data_primary, env_data_control and response data frames. Only years from
#' all three data frames are kept.
#' @param PCA_transformation if set to TRUE, all variables in the response
#' data frame will be transformed using PCA transformation.
#' @param log_preprocess if set to TRUE, variables will be transformed with
#' logarithmic transformation before used in PCA
#' @param components_selection character string specifying how to select the Principal
#' Components used as predictors.
#' There are three options: "automatic", "manual" and "plot_selection". If
#' argument is set to automatic, all scores with eigenvalues above 1 will be
#' selected. This threshold could be changed by changing the
#' eigenvalues_threshold argument. If parameter is set to "manual", user should
#' set the number of components with N_components argument. If components
#' selection is set to "plot_selection", Scree plot will be shown and a user must
#' manually enter the number of components to be used as predictors.
#' @param eigenvalues_threshold threshold for automatic selection of Principal Components
#' @param N_components number of Principal Components used as predictors
#' @param aggregate_function_env_data_primary character string specifying how the
#' daily data from env_data_primary should be aggregated. The default is 'mean',
#' the other options are 'median', 'sum', 'min' and 'max'
#' @param aggregate_function_env_data_control character string specifying how the
#' daily data from env_data_control should be aggregated. The default is 'mean',
#' the other options are 'median', 'sum', 'min' and 'max'
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
#' @param plot_specific_window integer representing window width to be displayed
#' for plot_specific
#' @param ylimits limit of the y axes for plot_extreme and plot_specific. It should be
#' given in the form of: ylimits = c(0,1)
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
#' @param day_interval a vector of two values: lower and upper time interval of
#' days that will be used to calculate statistical metrics. Negative values
#' indicate previous growing season days. This argument overwrites the calculation
#' limits defined by lower_limit and upper_limit arguments.
#' @param dc_method a character string to determine the method to detrend climate
#' (environmental) data.  Possible values are "none" (default), "SLD","Spline",
#' "ModNegExp", "Mean", "Friedman", "ModHugershoff". "SLD" refers to Simple
#' Linear Detrending
#' @param dc_nyrs a number giving the rigidity of the smoothing spline, defaults
#' to 0.67 of series length if nyrs is NULL (see dplR R package)
#' @param dc_f a number between 0 and 1 giving the frequency response or wavelength
#' cutoff. Defaults to 0.5 (see dplR R package)
#' @param dc_pos.slope a logical flag. Will allow for a positive slope to be used
#' in method "ModNegExp" and "ModHugershoff". If FALSE the line will be horizontal
#' (see dplR R package)
#' @param dc_constrain.nls a character string which controls the constraints of
#' the "ModNegExp" model and the "ModHugershoff"  (see dplR R package).
#' @param dc_span a numeric value controlling method "Friedman", or "cv" (default)
#' for automatic choice by cross-validation (see dplR R package)
#' @param dc_bass a numeric value controlling the smoothness of the fitted curve
#' in method "Friedman" (see dplR R package)
#' @param dc_difference	a logical flag. Compute residuals by subtraction if TRUE,
#' otherwise use division (see dplR R package)
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
#' @return a list with 15 elements:
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
#'  \item $plot_specific -  ggplot2 object: line plot of a row with a selected window width in a matrix of calculated metrics
#'  \item $PCA_output - princomp object: the result output of the PCA analysis
#'  \item $type - the character string describing type of analysis: daily or monthly
#'  \item $reference_window - character string, which reference window was used for calculations
#'  \item $aggregated_climate_primary - matrix with all aggregated climate series of primary data
#'  \item $aggregated_climate_control - matrix with all aggregated climate series of control data
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
#' data(LJ_daily_temperatures)
#' data(LJ_daily_precipitation)
#'
#' # 1 Basic example
#' example_basic <- daily_response_seascorr(response = data_MVA,
#'                           env_data_primary = LJ_daily_temperatures,
#'                           env_data_control = LJ_daily_precipitation,
#'                           row_names_subset = TRUE,
#'                           lower_limit = 35, upper_limit = 45,
#'                           remove_insignificant = FALSE,
#'                           aggregate_function_env_data_primary = 'median',
#'                           aggregate_function_env_data_control = 'median',
#'                           alpha = 0.05, pcor_method = "spearman",
#'                           tidy_env_data_primary = FALSE,
#'                           previous_year = FALSE, boot = TRUE,
#'                           tidy_env_data_control = TRUE, boot_n = 10,
#'                           reference_window = "end", k = 5,
#'                           dc_method = "SLD",
#'                           day_interval = c(-100, 250),
#'                           skip_window_length = 1,
#'                           skip_window_position = 1
#'                           )
#' summary(example_basic)
#' plot(example_basic, type = 1)
#' plot(example_basic, type = 2)
#' plot(example_basic, type = 3)
#' example_basic$optimized_return
#' example_basic$optimized_return_all
#' example_basic$temporal_stability
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
#'                           reference_window = "end")
#'
#' summary(example_fixed_width)
#' plot(example_fixed_width, type = 1)
#' plot(example_fixed_width, type = 2)
#' example_fixed_width$optimized_return
#' example_fixed_width$optimized_return_all
#'
#' }

daily_response_seascorr <- function(response, env_data_primary, env_data_control,
                           lower_limit = 30,
                           upper_limit = 90, fixed_width = 0,
                           previous_year = FALSE, pcor_method = "pearson",
                           remove_insignificant = TRUE,
                           alpha = .05, row_names_subset = FALSE,
                           PCA_transformation = FALSE, log_preprocess = TRUE,
                           components_selection = 'automatic',
                           eigenvalues_threshold = 1,
                           N_components = 2,
                           aggregate_function_env_data_primary = 'mean',
                           aggregate_function_env_data_control = 'mean',
                           temporal_stability_check = "sequential", k = 2,
                           k_running_window = 30,
                           subset_years = NULL, plot_specific_window = NULL,
                           ylimits = NULL, seed = NULL, tidy_env_data_primary = FALSE,
                           tidy_env_data_control = FALSE,
                           reference_window = 'start',  boot = FALSE, boot_n = 1000,
                           boot_ci_type = "norm", boot_conf_int = 0.95,
                           day_interval = ifelse(c(previous_year == TRUE,
                                                   previous_year == TRUE),
                                                 c(-1, 366), c(1, 366)),
                           dc_method = NULL,
                           dc_nyrs = NULL,
                           dc_f = 0.5,
                           dc_pos.slope = FALSE,
                           dc_constrain.nls = c("never", "when.fail", "always"),
                           dc_span = "cv",
                           dc_bass = 0,
                           dc_difference = FALSE,
                           pcor_na_use = "pairwise.complete",
                           skip_window_length = 1,
                           skip_window_position = 1){

  ##############################################################################
  # 1 day interval is organized
  offset_start <- day_interval[1]
  offset_end <- day_interval[2]

  # if both are positive but previous_year = TRUE
  if (offset_start > 0 & offset_end > 0 & previous_year == TRUE){

    previous_year <- FALSE

    warning(paste0("Previous year is not included in selected day_interval. ",
                   "The argument previous_year is set to FALSE"))
  }


  # if both are negative negative
  if (offset_start < 0 & offset_end < 0){
    offset_start <- abs(offset_start)
    offset_end <- abs(offset_end)

    # If previous_year is FALSE, we set it to TRUE
    if (previous_year == FALSE){
      previous_year = TRUE
      warning(paste0("Previous year is included in day_interval. ",
                     "The argument previous_year is set to TRUE"))
    }

    # if only offset_start is negative
  } else if (offset_start < 0 & offset_end > 0){
    offset_end <- offset_end + 366
    offset_start <- abs(offset_start)

    # If previous_year is FALSE, we set it to TRUE
    if (previous_year == FALSE){
      previous_year = TRUE
      warning(paste0("Previous year is included in day_interval. ",
                     "The argument previous_year is set to TRUE"))
    }

  }

  # Calculate the max_window allowed
  max_window <- offset_end - offset_start + 1

  # If max_window is greater then upper_limit, it must be reduced
  if (upper_limit > max_window){

    upper_limit <- max_window

    if (fixed_width == 0){
    warning(paste0("The upper_limit is outside your day_interval and",
                   " therefore reduced to the maximum allowed: ",max_window,"."))
    }
  }

  # Now, if upper_limit > max_window, we make them the same
  if (lower_limit > max_window){

    lower_limit <- max_window

    if (fixed_width == 0){
    warning(paste0("The lower_limit is outside your day_interval and",
                   " therefore reduced to the minimum allowed: ",max_window,"."))
    }
  }


  # Also correction for fixed_window approach
  if (fixed_width > max_window){

    stop(paste0("The selected fixed_width is outside your day_interval.",
                " Decrease the fixed_width argument to at least: ",max_window,"."))
  }


  if (previous_year == FALSE){

    offset_end <- 366 - offset_end

  } else {

    offset_end <- 732 - offset_end

  }

  ##############################################################################

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


 # If env_data_primary_control is given in tidy version, transformation is needed
 if (tidy_env_data_primary == TRUE){

   n_col_tidy_DF <- ncol(env_data_primary_control)
   colnames_tidy_DF <- colnames(env_data_primary_control)

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

  # PART 1 - general data arrangements, warnings and stops
  response <- data.frame(response)
  env_data_primary <- data.frame(env_data_primary)
  env_data_control <- data.frame(env_data_control)

  # Here we save the original env and response data that will be used later
  response_original <- response
  env_data_primary_original <- env_data_primary
  env_data_control_original <- env_data_control

  # For metric calculations, both objects need to have the same length,
  # with the exception, when row_names_subset is set to TRUE
  # Stop message in case both data frames do not have the same length
  if (nrow(response) !=  nrow(env_data_primary) & row_names_subset == FALSE)
    stop(paste0("Length of env_data_primary and response records differ",
                " You can use row_names_subset = TRUE"))

  if (nrow(response) !=  nrow(env_data_control) & row_names_subset == FALSE)
    stop(paste0("Length of env_data_control and response records differ",
                " You can use row_names_subset = TRUE"))


  #######################################################
  # Rules for previous_year = FALSE

  if (previous_year == FALSE){

    # Stop message if fixed_width is not between 0 and 366
    if (fixed_width < 0 | fixed_width > 366)
      stop("fixed_width should be between 1 and 366")

    if (lower_limit > upper_limit)
      stop("lower_limit can not be higher than upper_limit!")

    if (lower_limit > 366 | lower_limit < 1)
      stop("lower_limit out of bounds! It should be between 1 and 366")

    if (upper_limit > 366 | upper_limit < 1)
      stop("upper_limit out of bounds! It should be between 1 and 366")

    # Rules for previous_year = TRUE
  } else if (previous_year == TRUE){

    # Stop message if fixed_width is not between 0 and 366
    if (fixed_width < 0 | fixed_width > 732)
      stop("fixed_width should be between 1 and 732")

    if (lower_limit > upper_limit)
      stop("lower_limit can not be higher than upper_limit!")

    if (lower_limit > 732 | lower_limit < 1)
      stop("lower_limit out of bounds! It should be between 1 and 366")

    if (upper_limit > 732 | upper_limit < 1)
      stop("upper_limit out of bounds! It should be between 1 and 366")

  }

  # Make sure the selected method is appropriate
  if (!is.null(dc_method)){

    if (!(dc_method %in% c("Spline", "ModNegExp", "Mean", "Friedman", "ModHugershoff", "SLD"))){

      stop(paste0('dc_method should be one of Spline, ModNegExp, Mean, Friedman, ModHugershoff, SLD, but instead it is:',dc_method))

    }
  }

  # Warn users in case of missing values (selected threshold is 270 days)
  # A) env_data_primary
  env_temp_primary <- env_data_primary[row.names(env_data_primary) %in% row.names(response),]

  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    subset_seq <- seq(lower_subset, upper_subset)
    env_temp_primary <- subset(env_temp_primary, row.names(env_temp_primary) %in% subset_seq)
  }

  na_problem <- data.frame(na_sum = rowSums(is.na(env_temp_primary)))
  na_problem <- na_problem[na_problem$na_sum > 270, , F]
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
  na_problem <- na_problem[na_problem$na_sum > 270, , F]
  problematic_years <- paste0(row.names(na_problem), sep = "", collapse=", ")

  if (nrow(na_problem) > 0){

    warning(paste0("Problematic years with missing values are present in env_data_control: ", problematic_years))

  }

  # Data manipulation
  # If use.previous == TRUE, env_data_primary data has to be rearranged accordingly
  if (previous_year == TRUE) {

    # FIRST, both data frames need to be arranged, the most recent year is the first one
    env_data_primary$yearABC <- row.names(env_data_primary)
    env_data_primary <- dplyr::arrange(env_data_primary, desc(yearABC))
    env_data_primary <- years_to_rownames(env_data_primary, "yearABC")
    env_data_primary_previous <- env_data_primary[-1, , F]
    env_data_primary_current <- env_data_primary[-nrow(env_data_primary), ,F]
    row_names_current <- row.names(env_data_primary_current)
    env_data_primary <- cbind(env_data_primary_previous, env_data_primary_current)
    env_data_primary <- data.frame(env_data_primary)
    row.names(env_data_primary) <- row_names_current
    env_data_primary_original <- env_data_primary

    env_data_control$yearABC <- row.names(env_data_control)
    env_data_control <- dplyr::arrange(env_data_control, desc(yearABC))
    env_data_control <- years_to_rownames(env_data_control, "yearABC")
    env_data_control_previous <- env_data_control[-1, , F]
    env_data_control_current <- env_data_control[-nrow(env_data_control), ,F]
    row_names_current <- row.names(env_data_control_current)
    env_data_control <- cbind(env_data_control_previous, env_data_control_current)
    env_data_control <- data.frame(env_data_control)
    row.names(env_data_control) <- row_names_current
    env_data_control_original <- env_data_control

    }

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


  # If PCA_transformation = TRUE, PCA is performed
  if (PCA_transformation == TRUE) {

    # Logarithmic transformation before PCA
    if (log_preprocess == TRUE) {

      if (sum(response <= 0) > 1){
        stop("your response data contains negative observations. Please set the argument log_preprocess to FALSE")
      }

      response <- data.frame(log(response))
    }

    PCA_result <- princomp(response, cor = TRUE)

    if (components_selection == 'automatic'){
      subset_vector <- PCA_result$sdev > eigenvalues_threshold
      response <- as.data.frame(PCA_result$scores[, subset_vector])
    }

    if (components_selection == 'manual'){
      response <- as.data.frame(PCA_result$scores[, 1:N_components])
    }

    if (components_selection == 'plot_selection'){
      plot(PCA_result, type = 'l')

      fun <- function(){
        N_PC <- readline("What number of PC scores should be used as predictors? ")
        return(N_PC)
      }

      N_PC <- fun()
      response <- as.data.frame(PCA_result$scores[, 1:as.numeric(N_PC)])
    }

    number_PC <- ncol(response)
    df_names <-  paste( "PC_", seq(1:number_PC), sep = "")
    colnames(response) <- df_names

  } else (PCA_result <- "No PCA result avalialbe !")

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


  # PART 2 - the calculation of day-wise correlations
  # A) fixed window approach
  # B) Variable window approach

  # these are lists for climate and holder for saving mm1 and mm2
  list_climate_primary <- list()
  list_climate_control <- list()

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
      pb <- txtProgressBar(min = 0, max = (ncol(env_data_primary) - fixed_width - offset_end - offset_start + 1)/skip_window_length,
                           style = 3)
    }

    b = 0

      # An iterating loop. In each iteration x is calculated and represents
      # response (dependent) variable. X is a moving average. Window width of
      # a moving window is fixed_width. Next, partial correlation is stored in
      # temporal matrix.
    for (j in (seq((0 + offset_start -1), (ncol(env_data_primary) - max((fixed_width + offset_end), offset_end)), by = skip_window_position))) {

        b = b + 1

        if (aggregate_function_env_data_primary == 'median'){

          x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary),
                                 (1 + j): (j + fixed_width)]),1 , median, na.rm = TRUE)
        } else if (aggregate_function_env_data_primary == 'sum'){

          x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary),
                              (1 + j): (j + fixed_width)]),1 , sum, na.rm = TRUE)

        } else if (aggregate_function_env_data_primary == 'mean'){

          x1 <- rowMeans(data.frame(env_data_primary[1:nrow(env_data_primary),
                                 (1 + j): (j + fixed_width)]), na.rm = TRUE)

        } else if (aggregate_function_env_data_primary == 'min'){

          x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary),
                                                  (1 + j): (j + fixed_width)]),1 , min, na.rm = TRUE)

        } else if (aggregate_function_env_data_primary == 'max'){

          x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary),
                                                  (1 + j): (j + fixed_width)]),1 , max, na.rm = TRUE)

        } else {

          stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))

        }


        if (aggregate_function_env_data_control == 'median'){

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control),
                                       (1 + j): (j + fixed_width)]),1 , median, na.rm = TRUE)

        } else if (aggregate_function_env_data_control == 'sum'){

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control),
                                       (1 + j): (j + fixed_width)]),1 , sum, na.rm = TRUE)

        } else if (aggregate_function_env_data_control == 'mean'){

          x2 <- rowMeans(data.frame(env_data_control[1:nrow(env_data_control),
                                          (1 + j): (j + fixed_width)]), na.rm = TRUE)

        } else if (aggregate_function_env_data_control == 'min'){

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control),
                                                  (1 + j): (j + fixed_width)]),1 , min, na.rm = TRUE)

        } else if (aggregate_function_env_data_control == 'max'){

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control),
                                                  (1 + j): (j + fixed_width)]),1 , max, na.rm = TRUE)

        } else {

          stop(paste0("aggregate function for env_data_control  is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

        }


        if (!is.null(dc_method)){


          if (dc_method == "SLD"){

            tmp_model <- lm(x1 ~ seq(1:length(x1)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x1 - tmp_pred

            x1 <- data.frame(x1 = tmp_res/sd(tmp_res))

          } else {

            x1 <- dplR::detrend(data.frame(x1), method = dc_method, nyrs = dc_nyrs, f = dc_f,
                                pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                                span = dc_span, bass = dc_bass,  difference = dc_difference)}

        } else {

          x1 <- matrix(x1, nrow = nrow(env_data_primary), ncol = 1)

        }


        if (!is.null(dc_method)){

          if (dc_method == "SLD"){

            tmp_model <- lm(x2 ~ seq(1:length(x2)))
            tmp_pred <- predict(tmp_model)
            tmp_res <- x2 - tmp_pred

            x2 <- data.frame(x2 = tmp_res/sd(tmp_res))

          } else {

            x2 <- dplR::detrend(data.frame(x2), method = dc_method, nyrs = dc_nyrs, f = dc_f,
                                pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                                span = dc_span, bass = dc_bass,  difference = dc_difference)}

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

        if (fixed_width != max_window){setTxtProgressBar(pb, b)}
      }
    if (fixed_width != max_window){close(pb)}

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

      pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit)/skip_window_length, style = 3)
    }

    b = 0

  for (K in seq(lower_limit, upper_limit, by = skip_window_length)) {

    b = b + 1

    for (j in seq((0 + offset_start -1), (ncol(env_data_primary) - max((K + offset_end), offset_end)), by = skip_window_position)) {

      if (aggregate_function_env_data_primary == 'median'){

        if (K == 1){

          x1 <- env_data_primary[,K+j]

        } else {

          x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary), (1 + j) : (j + K)]),1 , median, na.rm = TRUE)}

      } else if (aggregate_function_env_data_primary == 'sum'){

        if (K == 1){

          x1 <- env_data_primary[,K+j]

        } else {

          x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary), (1 + j) : (j + K)]),1 , sum, na.rm = TRUE)}

        } else if (aggregate_function_env_data_primary == 'mean'){

        if (K == 1){

          x1 <- env_data_primary[,K+j]

        } else {

          x1 <- rowMeans(data.frame(env_data_primary[1:nrow(env_data_primary), (1 + j) : (j + K)]), na.rm = T)}

        } else if (aggregate_function_env_data_primary == 'min'){

          if (K == 1){

            x1 <- env_data_primary[,K+j]

          } else {

            x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary), (1 + j) : (j + K)]),1 , min, na.rm = TRUE)}

        } else if (aggregate_function_env_data_primary == 'max'){

          if (K == 1){

            x1 <- env_data_primary[,K+j]

          } else {

            x1 <- apply(data.frame(env_data_primary[1:nrow(env_data_primary), (1 + j) : (j + K)]),1 , max, na.rm = TRUE)}

      } else {

        stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))
      }


      if (aggregate_function_env_data_control == 'median'){

        if (K == 1){
          x2 <- env_data_control[,K+j]
        } else {

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control), (1 + j) : (j + K)]),1 , median, na.rm = TRUE)}

      } else if (aggregate_function_env_data_control == 'sum'){


        if (K == 1){
          x2 <- env_data_control[,K+j]
        } else {

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control), (1 + j) : (j + K)]),1 , sum, na.rm = TRUE)}

      } else if (aggregate_function_env_data_control == 'mean'){

        if (K == 1){
          x2 <- env_data_control[,K+j]
        } else {

          x2 <- rowMeans(data.frame(env_data_control[1:nrow(env_data_control), (1 + j) : (j + K)]), na.rm = T)}

      } else if (aggregate_function_env_data_control == 'min'){

        if (K == 1){
          x2 <- env_data_control[,K+j]
        } else {

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control), (1 + j) : (j + K)]),1 , min, na.rm = TRUE)}

      } else if (aggregate_function_env_data_control == 'max'){


        if (K == 1){
          x2 <- env_data_control[,K+j]
        } else {

          x2 <- apply(data.frame(env_data_control[1:nrow(env_data_control), (1 + j) : (j + K)]),1 , max, na.rm = TRUE)}

      } else {

        stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

      }

      if (!is.null(dc_method)){


        if (dc_method == "SLD"){

          tmp_model <- lm(x1 ~ seq(1:length(x1)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x1 - tmp_pred

          x1 <- data.frame(x1 = tmp_res/sd(tmp_res))

        } else {

          x1 <- dplR::detrend(data.frame(x1), method = dc_method, nyrs = dc_nyrs, f = dc_f,
                              pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                              span = dc_span, bass = dc_bass,  difference = dc_difference)}


      } else {

        x1 <- matrix(x1, nrow = nrow(env_data_primary), ncol = 1)

      }


      if (!is.null(dc_method)){


        if (dc_method == "SLD"){

          tmp_model <- lm(x2 ~ seq(1:length(x2)))
          tmp_pred <- predict(tmp_model)
          tmp_res <- x2 - tmp_pred

          x2 <- data.frame(x2 = tmp_res/sd(tmp_res))

        } else {

          x2 <- dplR::detrend(data.frame(x2), method = dc_method, nyrs = dc_nyrs, f = dc_f,
                              pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                              span = dc_span, bass = dc_bass,  difference = dc_difference)}

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
    if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
  }

    if (upper_limit != lower_limit){close(pb)}

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

  # To enhance the visualisation, insignificant values
  # are removed if remove_insignificant == TRUE
  if (remove_insignificant == TRUE){
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
                       plot_specific = NA,
                       PCA_output = PCA_result,
                       type = "daily",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate_primary = NA,
                       aggregated_climate_control = NA)

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

  # The fourth return element is being created: rowMeans/ apply of optimal sequence:
  # So, here we consider more options, based on the reference_winow
  # 1. reference window = "start"
  if (reference_window == 'start'){

  #x1
  if (aggregate_function_env_data_primary == 'median'){
    x1 <- data.frame(apply(data.frame(env_data_primary[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1), drop = FALSE]),1 , median, na.rm = TRUE))

  } else if (aggregate_function_env_data_primary == 'sum'){
    x1 <- data.frame(apply(data.frame(env_data_primary[, as.numeric(plot_column):
                                         (as.numeric(plot_column) +
                                            as.numeric(row_index) - 1), drop = FALSE]),1 , sum, na.rm = TRUE))

  } else if (aggregate_function_env_data_primary == 'mean'){
    x1 <- data.frame(rowMeans(data.frame(env_data_primary[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1), drop = FALSE]),
                                 na.rm = TRUE))

  } else if (aggregate_function_env_data_primary == 'min'){
    x1 <- data.frame(apply(data.frame(env_data_primary[, as.numeric(plot_column):
                                                         (as.numeric(plot_column) +
                                                            as.numeric(row_index) - 1), drop = FALSE]),1 , min, na.rm = TRUE))

  } else if (aggregate_function_env_data_primary == 'max'){
    x1 <- data.frame(apply(data.frame(env_data_primary[, as.numeric(plot_column):
                                                         (as.numeric(plot_column) +
                                                            as.numeric(row_index) - 1), drop = FALSE]),1 , max, na.rm = TRUE))
  } else {

    stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))

  }



 # x2
 if (aggregate_function_env_data_control == 'median'){
   x2 <- data.frame(apply(data.frame(env_data_control[, as.numeric(plot_column):
                                             (as.numeric(plot_column) +
                                                as.numeric(row_index) - 1), drop = FALSE]),1 , median, na.rm = TRUE))

 } else if (aggregate_function_env_data_control == 'sum'){
   x2 <- data.frame(apply(data.frame(env_data_control[, as.numeric(plot_column):
                                             (as.numeric(plot_column) +
                                                as.numeric(row_index) - 1), drop = FALSE]),1 , sum, na.rm = TRUE))

 } else if (aggregate_function_env_data_control == 'mean'){
   x2 <- data.frame(rowMeans(data.frame(env_data_control[, as.numeric(plot_column):
                                                (as.numeric(plot_column) +
                                                   as.numeric(row_index) - 1), drop = FALSE]), na.rm = TRUE))

 } else if (aggregate_function_env_data_control == 'min'){
   x2 <- data.frame(apply(data.frame(env_data_control[, as.numeric(plot_column):
                                                        (as.numeric(plot_column) +
                                                           as.numeric(row_index) - 1), drop = FALSE]),1 , min, na.rm = TRUE))

 } else if (aggregate_function_env_data_control == 'max'){
   x2 <- data.frame(apply(data.frame(env_data_control[, as.numeric(plot_column):
                                                        (as.numeric(plot_column) +
                                                           as.numeric(row_index) - 1), drop = FALSE]),1 , max, na.rm = TRUE))

 } else {
   stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))
 }


  ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
  # for the analysed period.

  if (aggregate_function_env_data_primary == 'median'){
    x1_original <- data.frame(apply(data.frame(env_data_primary_original[, as.numeric(plot_column):
                                         (as.numeric(plot_column) +
                                            as.numeric(row_index) - 1), drop = FALSE]),1 , median, na.rm = TRUE))
  } else if (aggregate_function_env_data_primary == 'sum'){
    x1_original <- data.frame(apply(data.frame(env_data_primary_original[, as.numeric(plot_column):
                                                           (as.numeric(plot_column) +
                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , sum, na.rm = TRUE))
  } else if (aggregate_function_env_data_primary == 'mean'){

     x1_original <- data.frame(apply(data.frame(env_data_primary_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , mean, na.rm = TRUE))

  } else if (aggregate_function_env_data_primary == 'min'){
    x1_original <- data.frame(apply(data.frame(env_data_primary_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , min, na.rm = TRUE))
  } else if (aggregate_function_env_data_primary == 'max'){

    x1_original <- data.frame(apply(data.frame(env_data_primary_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , max, na.rm = TRUE))
  } else {

    stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))
  }


  if (aggregate_function_env_data_control == 'median'){
    x2_original <- data.frame(apply(data.frame(env_data_control_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , median, na.rm = TRUE))
  } else if (aggregate_function_env_data_control == 'sum'){
    x2_original <- data.frame(apply(data.frame(env_data_control_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , sum, na.rm = TRUE))
  } else if (aggregate_function_env_data_control == 'mean'){

    x2_original <- data.frame(apply(data.frame(env_data_control_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , mean, na.rm = TRUE))

  } else if (aggregate_function_env_data_control == 'min'){
    x2_original <- data.frame(apply(data.frame(env_data_control_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , min, na.rm = TRUE))
  } else if (aggregate_function_env_data_control == 'max'){

    x2_original <- data.frame(apply(data.frame(env_data_control_original[, as.numeric(plot_column):
                                                                           (as.numeric(plot_column) +
                                                                              as.numeric(row_index) - 1), drop = FALSE]),1 , max, na.rm = TRUE))
  } else {

    stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

  }

}



  # Option 2, reference window = "end"
  if (reference_window == 'end'){

    if (aggregate_function_env_data_primary == 'median'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                           (as.numeric(plot_column)), drop = FALSE]),1 , median, na.rm = TRUE))
    } else if (aggregate_function_env_data_primary == 'sum'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                           (as.numeric(plot_column)), drop = FALSE]),1 , sum, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'mean'){

      x1 <- data.frame(apply(data.frame(env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                           (as.numeric(plot_column)), drop = FALSE]),1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'min'){

        x1 <- data.frame(apply(data.frame(env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                             (as.numeric(plot_column)), drop = FALSE]),1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'max'){

      x1 <- data.frame(apply(data.frame(env_data_primary[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                           (as.numeric(plot_column)), drop = FALSE]),1 , max, na.rm = TRUE))

    } else {
      stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))
    }


    if (aggregate_function_env_data_control == 'median'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                (as.numeric(plot_column)), drop = FALSE]),1 , median, na.rm = TRUE))
    } else if (aggregate_function_env_data_control == 'sum'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                (as.numeric(plot_column)), drop = FALSE]),1 , sum, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'mean'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                (as.numeric(plot_column)), drop = FALSE]),1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'min'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                           (as.numeric(plot_column)), drop = FALSE]),1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'max'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                                           (as.numeric(plot_column)), drop = FALSE]),1 , max, na.rm = TRUE))
    } else {

      stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

    }

    ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
    # for the analysed period.

    if (aggregate_function_env_data_primary == 'median'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column)), drop = FALSE]),1 , median, na.rm = TRUE))
    } else if (aggregate_function_env_data_primary == 'sum'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column)), drop = FALSE]),1 , sum, na.rm = TRUE))
    } else if (aggregate_function_env_data_primary == 'mean'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column)), drop = FALSE]),1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'min'){

      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column)), drop = FALSE]),1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'max'){

      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column)), drop = FALSE]),1 , max, na.rm = TRUE))

    } else {
      stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))
    }


    if (aggregate_function_env_data_control == 'median'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                                                as.numeric(plot_column)), drop = FALSE]),1 , median, na.rm = TRUE))
    } else if (aggregate_function_env_data_control == 'sum'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                                                as.numeric(plot_column)), drop = FALSE]),1 , sum, na.rm = TRUE))
    } else if (aggregate_function_env_data_control == 'mean'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                                                as.numeric(plot_column)), drop = FALSE]),1 , mean, na.rm = TRUE))
    } else if (aggregate_function_env_data_control == 'min'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                                                as.numeric(plot_column)), drop = FALSE]),1 , min, na.rm = TRUE))
    } else if (aggregate_function_env_data_control == 'max'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                                                as.numeric(plot_column)), drop = FALSE]),1 , max, na.rm = TRUE))
    } else {

      stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

    }
}


#  Third option: reference window = "middle"
  if (reference_window == 'middle'){

    if (as.numeric(row_index)%%2 == 0){
      adjustment_1 = 0
      adjustment_2 = 1
    } else {
      adjustment_1 = 1
      adjustment_2 = 2
    }

    if (aggregate_function_env_data_primary == 'median'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                          1 , median, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'sum'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                1 , sum, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'mean'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'min'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'max'){
      x1 <- data.frame(apply(data.frame(env_data_primary[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , max, na.rm = TRUE))

    } else {
      stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))
    }



    if (aggregate_function_env_data_control == 'median'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , median, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'sum'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , sum, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'mean'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'min'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'max'){
      x2 <- data.frame(apply(data.frame(env_data_control[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                             1 , max, na.rm = TRUE))
    } else {

      stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

    }


    ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
    # for the analysed period.

    if (aggregate_function_env_data_primary == 'median'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                1 , median, na.rm = TRUE))
    } else if (aggregate_function_env_data_primary == 'sum'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                         1 , sum, na.rm = TRUE))
    } else if (aggregate_function_env_data_primary == 'mean'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                         1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'min'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_primary == 'max'){
      x1_original <- data.frame(apply(data.frame(env_data_primary_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , max, na.rm = TRUE))

    } else {
      stop(paste0("aggregate function for env_data_primary is ", aggregate_function_env_data_primary, ". Instead it should be mean, median, sum, min or max."))
    }


    if (aggregate_function_env_data_control == 'median'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , median, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'sum'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , sum, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'mean'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , mean, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'min'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , min, na.rm = TRUE))

    } else if (aggregate_function_env_data_control == 'max'){
      x2_original <- data.frame(apply(data.frame(env_data_control_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2), drop = FALSE]),
                                      1 , max, na.rm = TRUE))
    } else {

      stop(paste0("aggregate function for env_data_control is ", aggregate_function_env_data_control, ". Instead it should be mean, median, sum, min or max."))

    }

  }

# Final output list

  if (!is.null(dc_method)){




    if (dc_method == "SLD"){

      x1 <- x1[,1]
      tmp_model <- lm(x1 ~ seq(1:length(x1)))
      tmp_pred <- predict(tmp_model)
      tmp_res <- x1 - tmp_pred
      x1 <- data.frame(x1 = tmp_res/sd(tmp_res))

      x2 <- x2[,1]
      tmp_model <- lm(x2 ~ seq(1:length(x2)))
      tmp_pred <- predict(tmp_model)
      tmp_res <- x2 - tmp_pred
      x2 <- data.frame(x2 = tmp_res/sd(tmp_res))

    } else {

      x1 <- dplR::detrend(x1, method = dc_method, nyrs = dc_nyrs, f = dc_f,
                          pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                          span = dc_span, bass = dc_bass,  difference = dc_difference)

      x2 <- dplR::detrend(x2, method = dc_method, nyrs = dc_nyrs, f = dc_f,
                          pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                          span = dc_span, bass = dc_bass,  difference = dc_difference)

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
      x1_original <- data.frame(x1_original = tmp_res/sd(tmp_res))

      x2_original <- x2_original[,1]
      tmp_model <- lm(x2_original ~ seq(1:length(x2_original)))
      tmp_pred <- predict(tmp_model)
      tmp_res <- x2_original - tmp_pred
      x2_original <- data.frame(x2_original = tmp_res/sd(tmp_res))

    } else {

      x1_original <- dplR::detrend(x1_original, method = dc_method, nyrs = dc_nyrs, f = dc_f,
                                   pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                                   span = dc_span, bass = dc_bass,  difference = dc_difference)

      x2_original <- dplR::detrend(x2_original, method = dc_method, nyrs = dc_nyrs, f = dc_f,
                                   pos.slope = dc_pos.slope, constrain.nls = dc_constrain.nls,
                                   span = dc_span, bass = dc_bass,  difference = dc_difference)

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
                       aggregated_climate_primary = do.call(cbind, list_climate_primary),
                       aggregated_climate_control = do.call(cbind, list_climate_control))


    plot_heatmapA <- plot_heatmap(final_list, reference_window = reference_window, type = "daily")
    plot_extremeA <- plot_extreme(final_list, ylimits = ylimits, reference_window = reference_window, type = "daily")

    width_sequence = seq(lower_limit, upper_limit)

    if (is.null(plot_specific_window)){
      (plot_specificA <- "plot_specific_window is not available. No plot_specific is made!")
    } else if (fixed_width != 0){

      if (fixed_width != plot_specific_window){
        warning(paste0("plot_specific_window and fixed_width differ!",
                       " fixed_wdith will be used to generate plot_specific!"))
      }

      plot_specific_window = fixed_width
      plot_specificA <- plot_specific(final_list, window_width = plot_specific_window, ylimits = ylimits,
                                      reference_window = reference_window)
    } else if (plot_specific_window %in% width_sequence){
      plot_specificA <- plot_specific(final_list, window_width = plot_specific_window, ylimits = ylimits,
                                      reference_window = reference_window)
    } else (plot_specificA <- "Selected plot_specific_window is not available. No plot_specific is made!")


      final_list <- list(calculations = temporal_matrix, method = "pcor",
                         metric = pcor_method, analysed_period = analysed_period,
                         optimized_return = x1_full,
                         optimized_return_all = x1_full_original,
                         transfer_function = p1, temporal_stability = temporal_stability,
                         cross_validation = NA,
                         plot_heatmap = plot_heatmapA,
                         plot_extreme = plot_extremeA,
                         plot_specific = plot_specificA,
                         PCA_output = PCA_result,
                         type = "daily",
                         reference_window = reference_window,
                         boot_lower = temporal_matrix_lower,
                         boot_upper = temporal_matrix_upper,
                         aggregated_climate_primary = do.call(cbind, list_climate_primary),
                         aggregated_climate_control = do.call(cbind, list_climate_control))

    class(final_list) <- "dmrs"

  return(final_list)
}
}
