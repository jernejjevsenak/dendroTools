#' monthly_response
#'
#' Function calculates all possible values of a selected statistical metric
#' between one or more response variables and monthly sequences of environmental
#' data. Calculations are based on moving window which slides through monthly
#' environmental data. All calculated metrics are stored in a matrix. The
#' location of stored calculated metric in the matrix is indicating a window
#' width (row names) and a location in a matrix of monthly sequences of
#' environmental data (column names).
#'
#' @param response a data frame with tree-ring proxy variables as columns and
#' (optional) years as row names. Row.names should be matched with those from a
#' env_data data frame. If not, set row_names_subset = TRUE.
#' @param env_data a data frame of monthly sequences of environmental
#' data as columns and years as row names. Each row represents a year and
#' each column represents a day of a year (or month). Row.names should be matched with
#' those from a response data frame. If not, set row_names_subset = TRUE.
#' Alternatively, env_data could be a tidy data with three columns,
#' i.e. Year, DOY (Month) and third column representing values of mean temperatures,
#' sum of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data to TRUE.
#' @param lower_limit lower limit of window width (i.e. number of consecutive months
#' to be used for calculations)
#' @param upper_limit upper limit of window width (i.e. number of consecutive months
#' to be used for calculations)
#' @param fixed_width fixed width used for calculations (i.e. number of consecutive
#' months to be used for calculations)
#' @param method a character string specifying which method to use. Current
#' possibilities are "cor", "lm" and "brnn".
#' @param metric a character string specifying which metric to use. Current
#' possibilities are "r.squared" and "adj.r.squared". If method = "cor",
#' metric is not relevant.
#' @param cor_method a character string indicating which correlation
#' coefficient is to be computed. One of "pearson" (default), "kendall", or
#' "spearman".
#' @param previous_year if set to TRUE, env_data and response variables will be
#' rearranged in a way, that also previous year will be used for calculations of
#' selected statistical metric.
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
#' @param aggregate_function character string specifying how the monthly data should be
#' aggregated. The default is 'mean', the two other options are 'median' and 'sum'
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
#' @param plot_specific_window integer representing window width to be displayed
#' for plot_specific
#' @param ylimits limit of the y axes for plot_extreme and plot_specific. It should be
#' given in the form of: ylimits = c(0,1)
#' @param seed optional seed argument for reproducible results
#' @param tidy_env_data if set to TRUE, env_data should be inserted as a data frame with three
#' columns: "Year", "Month", "Precipitation/Temperature/etc."
#' @param boot logical, if TRUE, bootstrap procedure will be used to calculate
#' estimates correlation coefficients, R squared or adjusted R squared metrices
#' @param boot_n The number of bootstrap replicates
#' @param boot_ci_type A character string representing the type of bootstrap intervals
#' required. The value should be any subset of the values c("norm","basic", "stud",
#' "perc", "bca").
#' @param boot_conf_int A scalar or vector containing the confidence level(s) of
#' the required interval(s)
#' @param month_interval a vector of two values: lower and upper time interval of
#' months that will be used to calculate statistical metrics. Negative values indicate
#' previous growing season months. This argument overwrites the calculation
#' limits defined by lower_limit and upper_limit arguments.
#'
#' @return a list with 17 elements:
#' \enumerate{
#'  \item $calculations - a matrix with calculated metrics
#'  \item $method - the character string of a method
#'  \item $metric - the character string indicating the metric used for calculations
#'  \item $analysed_period - the character string specifying the analysed period based on the information from row names. If there are no row names, this argument is given as NA
#'  \item $optimized_return - data frame with two columns, response variable and aggregated (averaged) monthly data that return the optimal results. This data.frame could be directly used to calibrate a model for climate reconstruction
#'  \item $optimized_return_all - a data frame with aggregated monthly data, that returned the optimal result for the entire env_data (and not only subset of analysed years)
#'  \item $transfer_function - a ggplot object: scatter plot of optimized return and a transfer line of the selected method
#'  \item $temporal_stability - a data frame with calculations of selected metric for different temporal subsets
#'  \item $cross_validation - a data frame with cross validation results
#'  \item $plot_heatmap - ggplot2 object: a heatmap of calculated metrics
#'  \item $plot_extreme - ggplot2 object: line or bar plot of a row with the highest value in a matrix of calculated metrics
#'  \item $plot_specific - not available for monthly_response()
#'  \item $PCA_output - princomp object: the result output of the PCA analysis
#'  \item $type - the character string describing type of analysis: daily or monthly
#'  \item $reference_window - character string, which reference window was used for calculations
#'  \item $boot_lower - matrix with lower limit of confidence intervals of bootstrap calculations
#'  \item $boot_upper - matrix with upper limit of confidence intervals of bootstrap calculations
#'}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the dendroTools R package
#' library(dendroTools)
#'
#' # Load data used for examples
#' data(data_MVA)
#' data(data_TRW)
#' data(data_TRW_1)
#' data(example_proxies_individual)
#' data(example_proxies_1)
#' data(LJ_monthly_temperatures)
#' data(LJ_monthly_precipitation)
#'
#' # 1 Example with tidy precipitation data
#' example_tidy_data <- monthly_response(response = data_MVA,
#'     lower_limit = 1, upper = 12,
#'     env_data = LJ_monthly_precipitation, fixed_width = 0,
#'     method = "cor", row_names_subset = TRUE, metric = "adj.r.squared",
#'     remove_insignificant = TRUE, previous_year = FALSE,
#'     alpha = 0.05, aggregate_function = 'sum', boot = TRUE,
#'     tidy_env_data = TRUE, boot_n = 100)
#'
#' summary(example_tidy_data)
#' plot(example_tidy_data, type = 1)
#' plot(example_tidy_data, type = 2)
#'
#' # 2 Example with splited data for past and present
#' example_MVA_past <- monthly_response(response = data_MVA,
#'     env_data = LJ_monthly_temperatures,
#'     method = "cor", row_names_subset = TRUE, previous_year = TRUE,
#'     remove_insignificant = TRUE, alpha = 0.05,
#'     subset_years = c(1940, 1980), aggregate_function = 'mean')
#'
#' example_MVA_present <- monthly_response(response = data_MVA,
#'     env_data = LJ_monthly_temperatures,
#'     method = "cor", row_names_subset = TRUE, alpha = 0.05,
#'     previous_year = TRUE, remove_insignificant = TRUE,
#'     subset_years = c(1981, 2010), aggregate_function = 'mean')
#'
#' summary(example_MVA_present)
#' plot(example_MVA_past, type = 1)
#' plot(example_MVA_present, type = 1)
#' plot(example_MVA_past, type = 2)
#' plot(example_MVA_present, type = 2)
#'
#'
#' # 3 Example with principal component analysis
#' example_PCA <- monthly_response(response = example_proxies_individual,
#'    env_data = LJ_monthly_temperatures, method = "lm",
#'    row_names_subset = TRUE, remove_insignificant = TRUE,
#'    alpha = 0.01, PCA_transformation = TRUE, previous_year = TRUE,
#'    components_selection = "manual", N_components = 2, boot = TRUE)
#'
#' summary(example_PCA$PCA_output)
#' plot(example_PCA, type = 1)
#' plot(example_PCA, type = 2)
#'
#' # 4 Example negative correlations
#' example_neg_cor <- monthly_response(response = data_TRW_1, alpha = 0.05,
#'    env_data = LJ_monthly_temperatures,
#'    method = "cor", row_names_subset = TRUE,
#'    remove_insignificant = TRUE, boot = TRUE)
#'
#' summary(example_neg_cor)
#' plot(example_neg_cor, type = 1)
#' plot(example_neg_cor, type = 2)
#' example_neg_cor$temporal_stability
#'
#' # 5 Example of multiproxy analysis
#' summary(example_proxies_1)
#' cor(example_proxies_1)
#'
#' example_multiproxy <- monthly_response(response = example_proxies_1,
#'    env_data = LJ_monthly_temperatures,
#'    method = "lm", metric = "adj.r.squared",
#'    row_names_subset = TRUE, previous_year = FALSE,
#'    remove_insignificant = TRUE, alpha = 0.05)
#'
#' summary(example_multiproxy)
#' plot(example_multiproxy, type = 1)
#'
#' # 6 Example to test the temporal stability
#' example_MVA_ts <- monthly_response(response = data_MVA,
#'    env_data = LJ_monthly_temperatures,
#'    method = "lm", metric = "adj.r.squared", row_names_subset = TRUE,
#'    remove_insignificant = TRUE, alpha = 0.05,
#'    temporal_stability_check = "running_window", k_running_window = 10)
#'
#' summary(example_MVA_ts)
#' example_MVA_ts$temporal_stability
#'
#' }

monthly_response <- function(response, env_data, method = "cor",
                           metric = "r.squared", cor_method = "pearson",
                           previous_year = FALSE, neurons = 1,
                           lower_limit = 1, upper_limit = 12,
                           fixed_width = 0, brnn_smooth = TRUE,
                           remove_insignificant = TRUE,
                           alpha = .05, row_names_subset = FALSE,
                           PCA_transformation = FALSE, log_preprocess = TRUE,
                           components_selection = 'automatic',
                           eigenvalues_threshold = 1,
                           N_components = 2, aggregate_function = 'mean',
                           temporal_stability_check = "sequential", k = 2,
                           k_running_window = 30, cross_validation_type = "blocked",
                           subset_years = NULL, plot_specific_window = NULL,
                           ylimits = NULL, seed = NULL, tidy_env_data = FALSE,
                           boot = FALSE, boot_n = 1000, boot_ci_type = "norm",
                           boot_conf_int = 0.95,
                           month_interval = ifelse(c(previous_year == TRUE,
                                                     previous_year == TRUE),
                                                   c(-1, 12), c(1, 12))) {

  ##############################################################################
  # 1 day interval is organized
  offset_start <- month_interval[1]
  offset_end <- month_interval[2]

  # if both are positive but previous_year = TRUE
  if (offset_start > 0 & offset_end > 0 & previous_year == TRUE){

    previous_year <- FALSE

    warning(paste0("Previous year is not included in selected month_interval. ",
                   "The argument previous_year is set to FALSE"))
  }


  # if both are negative negative
  if (offset_start < 0 & offset_end < 0){
    offset_start <- abs(offset_start)
    offset_end <- abs(offset_end)

    # If previous_year is FALSE, we set it to TRUE
    if (previous_year == FALSE){
      previous_year = TRUE
      warning(paste0("Previous year is included in month_interval. ",
                     "The argument previous_year is set to TRUE"))
    }

    # if only offset_start is negative
  } else if (offset_start < 0 & offset_end > 0){
    offset_end <- offset_end + 12
    offset_start <- abs(offset_start)

    # If previous_year is FALSE, we set it to TRUE
    if (previous_year == FALSE){
      previous_year = TRUE
      warning(paste0("Previous year is included in month_interval. ",
                     "The argument previous_year is set to TRUE"))
    }

  }

  # Calculate the max_window allowed
  max_window <- offset_end - offset_start + 1

  # If max_window is greater then upper_limit, it must be reduced
  if (upper_limit > max_window){

    upper_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The upper_limit is outside your month_interval and",
                     " therefore reduced to the maximum allowed: ",max_window,"."))
    }

  }

  # Now, if lower_limit > max_window, we make them the same
  if (lower_limit > max_window){
    lower_limit <- max_window

    if (fixed_width == 0){
    warning(paste0("The lower_limit is outside your month_interval and",
                   " therefore reduced to the minimum allowed: ",max_window,"."))
    }
  }


  # Also correction for fixed_window approach - can only be ERROR
  if (fixed_width > max_window){

    stop(paste0("The selected fixed_width is outside your month_interval.",
                " Decrease the fixed_width argument to at least: ",max_window,"."))
  }

  if (previous_year == FALSE){

    offset_end <- 12 - offset_end

  } else {

    offset_end <- 24 - offset_end

  }

################################################################################

if (!is.null(seed)) {
    set.seed(seed)
  }

if (fixed_width != 0){
    lower_limit = 1
    upper_limit = 12
}


  if (fixed_width > 12 & previous_year == FALSE){
    stop(paste0("fixed_width argument can not be greater than 12! Instead, it is ", fixed_width, "!"))
  }

  reference_window = 'start'


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


 # If you insert daily data into monthly response, you get an error
 if ((ncol(env_data) != 12) && (tidy_env_data == FALSE)){
   stop(paste0("You must insert env_data with 12 columns (months)! Instead, you have ", ncol(env_data),
               " columns!"))
 }

 if ((upper_limit > 24 & previous_year == TRUE) | (upper_limit > 12 & previous_year == FALSE))
   stop("upper_limit out of bounds!")

 if (lower_limit < 1)
   stop("lower_limit must be positive!")

 if (upper_limit > 24 | upper_limit < 1)
   stop("upper_limit out of bounds! It should be between 1 and 12 (24)")



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

   if (colnames_tidy_DF[2] != "Month"){
     stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                "The second column name of the env_data should be 'DOY', but it is",
                colnames_tidy_DF[2], "instead!"))
   }

   value_variable = colnames(env_data)[3]
   env_data <- dcast(env_data, Year~Month, value.var = value_variable)
   env_data <- years_to_rownames(env_data, "Year")



 }


  # PART 1 - general data arrangements, warnings and stops
  # Both bojects (response and env_data) are converted to data frames
  response <- data.frame(response)
  env_data <- data.frame(env_data)

  # Here we save the original env and response data that will be used later
  response_original <- response
  env_data_original <- env_data


  # For metric calculations, both objects need to have the same length,
  # with the exception, when row_names_subset is set to TRUE
  # Stop message in case both data frames do not have the same length
  if (nrow(response) !=  nrow(env_data) & row_names_subset == FALSE)
    stop("Length of env_data and response records differ")

  # Stop message if fixed_width is not between 0 and 24
  if (fixed_width < 0 | fixed_width > 24)
    stop("fixed_width should be between 1 and 24")

  # Stop in case of method == "cor" and ncol(proxies) > 1
  # Correlations could be calculated only for one variable
  if (method == "cor" & ncol(response) > 1)
    stop(paste("More than 1 variable in response data frame not suitable ",
  "for 'cor' method. Use 'lm' or 'brnn'"))


  if (previous_year == FALSE){

    # Stop message if fixed_width is not between 0 and 366
    if (fixed_width < 0 | fixed_width > 12)
      stop("fixed_width should be between 1 and 12")

    if (lower_limit > upper_limit)
      stop("lower_limit can not be higher than upper_limit!")

    if (lower_limit > 12 | lower_limit < 1)
      stop("lower_limit out of bounds! It should be between 1 and 12")

    if (upper_limit > 12 | upper_limit < 1)
      stop("upper_limit out of bounds! It should be between 1 and 12")

    # Rules for previous_year = TRUE
  } else if (previous_year == TRUE){

    # Stop message if fixed_width is not between 0 and 366
    if (fixed_width < 0 | fixed_width > 24)
      stop("fixed_width should be between 1 and 24")

    if (lower_limit > upper_limit)
      stop("lower_limit can not be higher than upper_limit!")

    if (lower_limit > 24 | lower_limit < 1)
      stop("lower_limit out of bounds! It should be between 1 and 24")

    if (upper_limit > 24 | upper_limit < 1)
      stop("upper_limit out of bounds! It should be between 1 and 24")

  }


  # Data manipulation
  # If use.previous == TRUE, env_data data has to be rearranged accordingly
  if (previous_year == TRUE) {

    # FIRST, both data frames need to be arranged, the most recent year is the first one
    env_data$yearABC <- row.names(env_data)
    env_data <- dplyr::arrange(env_data, desc(yearABC))
    env_data <- years_to_rownames(env_data, "yearABC")
    env_data_previous <- env_data[-1, , F]
    env_data_current <- env_data[-nrow(env_data), ,F]
    row_names_current <- row.names(env_data_current)
    env_data <- cbind(env_data_previous, env_data_current)
    env_data <- data.frame(env_data)
    row.names(env_data) <- row_names_current
    env_data_original <- env_data

    response$yearABC <- row.names(response)
    response <- dplyr::arrange(response, desc(yearABC))
    response <- years_to_rownames(response, "yearABC")
    response <- data.frame(response[-nrow(response),,F ])
    response <- data.frame(response)
    response_original <- response

    }

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

    if (any(!(subset_seq %in% row.names(env_data)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  " in the env_data data frame. Change the subset_years argument"))
    }

    response <- subset(response, row.names(response) %in% subset_seq)
    env_data <- subset(env_data, row.names(env_data) %in% subset_seq)
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
        pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width - offset_end - offset_start + 1),
                           style = 3)}

      b = 0

      # An iterating loop. In each itteration x is calculated and represents
      # response (dependent) variable. X is a moving average. Window width of
      # a moving window is fixed_width. Next, statistical metric is calculated
      # based on a selected method (cor, lm or brnn). Calculation is stored in
      # temporal matrix.
        for (j in (0 + offset_start -1): (ncol(env_data) - max((fixed_width + offset_end), offset_end))) {

        b = b + 1

        if (aggregate_function == 'median'){

          if (fixed_width == 1){
            x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
          } else {
            x <- apply(env_data[1:nrow(env_data),
                                (1 + j): (j + fixed_width)],1 , median, na.rm = TRUE)
            }

        } else if (aggregate_function == 'sum'){

          if (fixed_width == 1){
            x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
          } else {
            x <- apply(env_data[1:nrow(env_data),
                                (1 + j): (j + fixed_width)],1 , sum, na.rm = TRUE)
          }

        } else if (aggregate_function == 'mean'){

          if (fixed_width == 1){
            x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
          } else {
            x <- rowMeans(env_data[1:nrow(env_data),
                                   (1 + j): (j + fixed_width)], na.rm = TRUE)
          }

        } else {
          stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
        }

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)

        if (boot == FALSE){

          temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method)
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

        if (fixed_width != max_window){setTxtProgressBar(pb, b)}
      }

      if (fixed_width != max_window){close(pb)}

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
      pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width - offset_end - offset_start + 1),
                           style = 3)}

    b = 0

    for (j in (0 + offset_start -1): (ncol(env_data) - max((fixed_width + offset_end), offset_end))) {

      b = b + 1

      if (aggregate_function == 'median'){

        if (fixed_width == 1){
          x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
        } else {

        x <- apply(env_data[1:nrow(env_data),
                               (1 + j) : (j + fixed_width)],1 , median, na.rm = TRUE) }
      } else if (aggregate_function == 'sum'){

        if (fixed_width == 1){
          x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
        } else {
        x <- apply(env_data[1:nrow(env_data),
                            (1 + j) : (j + fixed_width)],1 , median, na.rm = TRUE)

      }
        } else if (aggregate_function == 'mean'){

          if (fixed_width == 1){
            x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
          } else {

        x <- rowMeans(env_data[1:nrow(env_data),
                               (1 + j) : (j + fixed_width)], na.rm = TRUE)
          }
        } else {
        stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
      }

      x <- matrix(x, nrow = nrow(env_data), ncol = 1)

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

      if (fixed_width != max_window){setTxtProgressBar(pb, b)}

    }

    if (fixed_width != max_window){close(pb)}

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
      pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width - offset_end - offset_start + 1),
                           style = 3)}

    b = 0

    for (j in (0 + offset_start -1): (ncol(env_data) - max((fixed_width + offset_end), offset_end))) {

       b = b + 1

        if (aggregate_function == 'median'){

          if (fixed_width == 1){
            x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
          } else {

         x <- apply(env_data[1:nrow(env_data),
                                (1 + j): (j + fixed_width)],1 , median, na.rm = TRUE)
          }
        } else if (aggregate_function == 'sum'){

          if (fixed_width == 1){
            x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
          } else {

          x <- apply(env_data[1:nrow(env_data),
                              (1 + j): (j + fixed_width)],1 , sum, na.rm = TRUE)
          }

       } else if (aggregate_function == 'mean') {

         if (fixed_width == 1){
           x <- env_data[1:nrow(env_data), (1 + j): (j + fixed_width)]
         } else {

         x <- rowMeans(env_data[1:nrow(env_data),
                                (1 + j): (j + fixed_width)], na.rm = TRUE)
         }

       } else {
         stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
       }

      x <- matrix(x, nrow = nrow(env_data), ncol = 1)

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
                                              ncol(as.data.frame(x))
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

      if (fixed_width != max_window){setTxtProgressBar(pb, b)}
     }

    if (fixed_width != max_window){close(pb)}

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
      pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)
    }

  b = 0


  for (K in lower_limit:upper_limit) {

    b = b + 1

    for (j in (0 + offset_start -1): (ncol(env_data) - max((K + offset_end), offset_end))) {

      if (aggregate_function == 'median'){
        if (K == 1){
          x <- env_data[,K+j]
        } else {
           x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , median, na.rm = TRUE)}
      } else if (aggregate_function == 'sum'){

        if (K == 1){
          x <- env_data[,K+j]
        } else {

        x <- apply(data.frame(env_data[1:nrow(env_data), (1 + j) : (j + K)]),1 , sum, na.rm = TRUE)}
        }
      else if (aggregate_function == 'mean'){

        if (K == 1){
          x <- env_data[,K+j]
        } else {

        x <- rowMeans(data.frame(env_data[1:nrow(env_data), (1 + j) : (j + K)]), na.rm = T)}
      } else {
        stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
      }

      x <- matrix(x, nrow = nrow(env_data), ncol = 1)


      if (boot == FALSE){
        temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method)
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
    if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
  }

  if (upper_limit != lower_limit){close(pb)}

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
      pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)
    }

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1

      for (j in (0 + offset_start -1): (ncol(env_data) - max((K + offset_end), offset_end))) {

        if (aggregate_function == 'median'){

          if (K == 1){
            x <- env_data[,K+j]
          } else {

          x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , median, na.rm = TRUE)}
        } else if(aggregate_function == 'sum'){
          if (K == 1){
            x <- env_data[,K+j]
          } else {
          x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , sum, na.rm = TRUE)}
        } else if (aggregate_function == 'mean'){

          if (K == 1){
            x <- env_data[,K+j]
          } else {

          x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + K)], na.rm = T)}

        } else {
          stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
        }

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)

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
      if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
    }

    if (upper_limit != lower_limit){close(pb)}

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
      pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)
    }

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1


      for (j in (0 + offset_start -1): (ncol(env_data) - max((K + offset_end), offset_end))) {

        if (aggregate_function == 'median'){

          if (K == 1){
            x <- env_data[,K+j]
          } else {

          x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , median, na.rm = TRUE)}
        } else if (aggregate_function == 'sum'){

          if (K == 1){
            x <- env_data[,K+j]
          } else {

          x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , sum, na.rm = TRUE)}

        } else if (aggregate_function == 'mean'){

          if (K == 1){
            x <- env_data[,K+j]
          } else {

          x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + K)], na.rm = T)}
        } else {
          stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
        }

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)

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
                                                ncol(as.data.frame(x))
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
      if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
    }

    if (upper_limit != lower_limit){close(pb)}

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

    if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){
      stop("All calculations are insignificant! Please change the alpha argument.")
    }


  }

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
  # 1. reference window = "start"
  if (reference_window == 'start'){


  if (aggregate_function == 'median'){
    dataf <- data.frame(apply(data.frame(env_data[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1)]),1 , median, na.rm = TRUE))

  } else if (aggregate_function == 'sum'){
    dataf <- data.frame(apply(data.frame(env_data[, as.numeric(plot_column):
                                         (as.numeric(plot_column) +
                                            as.numeric(row_index) - 1)]),1 , sum, na.rm = TRUE))

  } else if (aggregate_function == 'mean'){
    dataf <- data.frame(rowMeans(data.frame(env_data[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1)]),
                                 na.rm = TRUE))
  } else {
    stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
  }

  dataf_full <- cbind(response, dataf)
  colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
  colnames(dataf) <- "Optimized.rowNames"

  ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
  # for the analysed period.

  if (aggregate_function == 'median'){
    dataf_original <- data.frame(apply(data.frame(env_data_original[, as.numeric(plot_column):
                                         (as.numeric(plot_column) +
                                            as.numeric(row_index) - 1)]),1 , median, na.rm = TRUE))
  } else if (aggregate_function == 'sum'){
    dataf_original <- data.frame(apply(data.frame(env_data_original[, as.numeric(plot_column):
                                                           (as.numeric(plot_column) +
                                                              as.numeric(row_index) - 1)]),1 , sum, na.rm = TRUE))
  } else if (aggregate_function == 'mean'){
    dataf_original <- data.frame(rowMeans(data.frame(env_data_original[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1)]),
                                 na.rm = TRUE))
  } else {
    stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
  }

  dataf_full_original <- dataf_original
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
    optimized_result <- cor(dataf, response, method = cor_method)
  }

  # Just give a nicer colname
  colnames(dataf) <- "Optimized return"

  }

  # Option 2, reference window = "end"
    if (reference_window == 'end'){

    if (aggregate_function == 'median'){
      dataf <- data.frame(apply(env_data[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                           (as.numeric(plot_column))],1 , median, na.rm = TRUE))
    } else if (aggregate_function == 'sum'){
      dataf <- data.frame(apply(env_data[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                           (as.numeric(plot_column))],1 , sum, na.rm = TRUE))

    } else if (aggregate_function == 'mean'){
      dataf <- data.frame(apply(env_data[, (as.numeric(plot_column) - as.numeric(row_index) + 1):
                                           (as.numeric(plot_column))],1 , mean, na.rm = TRUE))
    } else {
      stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
    }

    dataf_full <- cbind(response, dataf)
    colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
    colnames(dataf) <- "Optimized.rowNames"

    ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
    # for the analysed period.

    if (aggregate_function == 'median'){
      dataf_original <- data.frame(apply(env_data_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column))],1 , median, na.rm = TRUE))
    } else if (aggregate_function == 'sum'){
      dataf_original <- data.frame(apply(env_data_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column))],1 , sum, na.rm = TRUE))
    } else if (aggregate_function == 'mean'){
      dataf_original <- data.frame(apply(env_data_original[, (as.numeric(plot_column) - (as.numeric(row_index) + 1):
                                                             as.numeric(plot_column))],1 , mean, na.rm = TRUE))
    } else {
      stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
    }

    dataf_full_original <- dataf_original
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
      optimized_result <- cor(dataf, response, method = cor_method)
    }

    # Just give a nicer colname
    colnames(dataf) <- "Optimized return"

  }

  # 1. reference window = "middle"
  if (reference_window == 'middle'){

    if (as.numeric(row_index)%%2 == 0){
      adjustment_1 = 0
      adjustment_2 = 1
    } else {
      adjustment_1 = 1
      adjustment_2 = 2
    }

    if (aggregate_function == 'median'){
      dataf <- data.frame(apply(env_data[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2)],
                                          1 , median, na.rm = TRUE))

    } else if (aggregate_function == 'sum'){
      dataf <- data.frame(apply(env_data[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2)],
                                1 , sum, na.rm = TRUE))

    } else if (aggregate_function == 'mean'){
      dataf <- data.frame(apply(env_data[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2)],
                                1 , mean, na.rm = TRUE))
    } else {
      stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
    }

    dataf_full <- cbind(response, dataf)
    colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
    colnames(dataf) <- "Optimized.rowNames"

    ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
    # for the analysed period.

    if (aggregate_function == 'median'){
      dataf_original <- data.frame(apply(env_data_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                           (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2)],
                                1 , median, na.rm = TRUE))
    } else if (aggregate_function == 'sum'){
      dataf_original <- data.frame(apply(env_data_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2)],
                                         1 , sum, na.rm = TRUE))
    } else if (aggregate_function == 'mean'){
      dataf_original <- data.frame(apply(env_data_original[, (round2((as.numeric(plot_column) - (as.numeric(row_index))/2)) - adjustment_1):
                                                             (round2((as.numeric(plot_column) + as.numeric(row_index)/2)) - adjustment_2)],
                                         1 , mean, na.rm = TRUE))
    } else {
      stop(paste0("aggregate function is ", aggregate_function, ". Instead it should be mean, median or sum."))
    }

    dataf_full_original <- dataf_original
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
      optimized_result <- cor(dataf, response, method = cor_method)
    }

    # Just give a nicer colname
    colnames(dataf) <- "Optimized return"

  }


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
  full_range = data.frame(proxy = seq(from = min(response[,1]), to = max(response[,1]), length.out = 100))

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

  dataset = data.frame(optimized_return =dataf[,1], proxy = response)

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
        calculation <- cor(dataset_temp[,1], dataset_temp[,2])
        sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exaxt= FALSE)$p.value
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
          calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method)
          sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exaxt= FALSE)$p.value
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
    calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method)
    sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exaxt= FALSE)$p.value
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
                       cross_validation = cross_validation)
  }

  if (method == "cor"){
    final_list <- list(calculations = temporal_matrix, method = method,
                       metric = cor_method, analysed_period = analysed_period,
                       optimized_return = dataf_full,
                       optimized_return_all = dataf_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = cross_validation)
  }



    plot_heatmapA <- plot_heatmap(final_list, reference_window = reference_window, type = "monthly")
    plot_extremeA <- plot_extreme(final_list, ylimits = ylimits, reference_window = reference_window, type = "monthly")

    width_sequence = seq(lower_limit, upper_limit)

    if (is.null(plot_specific_window)){
      (plot_specificA <- "plot_specific_window is not available for monthly_response function!")
    } else if (fixed_width != 0){

      if (fixed_width != plot_specific_window){
        warning(paste0("plot_specific_window and fixed_width differ!",
                       " fixed_wdith will be used to generate plot_specific!"))
      }

      plot_specific_window = fixed_width
      plot_specificA <- plot_specific(final_list, window_width = plot_specific_window, ylimits = ylimits,
                                      reference_window = reference_window, type = "monthly")
    } else if (plot_specific_window %in% width_sequence){
      plot_specificA <- plot_specific(final_list, window_width = plot_specific_window, ylimits = ylimits,
                                      reference_window = reference_window, type = "monthly")
    } else (plot_specificA <- "Selected plot_specific_window is not available. No plot_specific is made!")

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
                         plot_specific = plot_specificA,
                         PCA_output = PCA_result,
                         type = "monthly",
                         reference_window = reference_window,
                         boot_lower = temporal_matrix_lower,
                         boot_upper = temporal_matrix_upper)
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
                         plot_specific = plot_specificA,
                         PCA_output = PCA_result,
                         type = "monthly",
                         reference_window = reference_window,
                         boot_lower = temporal_matrix_lower,
                         boot_upper = temporal_matrix_upper)
    }

    class(final_list) <- 'dmrs'

  return(final_list)
}
