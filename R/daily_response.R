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
#' possibilities are "cor", "lm" and "brnn".
#' @param metric a character string specifying which metric to use. Current
#' possibilities are "r.squared" and "adj.r.squared". If method = "cor",
#' metric is not relevant.
#' @param lower_limit lower limit of window width
#' @param upper_limit upper limit of window width
#' @param fixed_width fixed width used for calculation. If fixed_width is
#' assigned a value, upper_limit and lower_limit will be ignored
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
#' @param use_median if set to TRUE, median will be used instead of mean to calculate
#' means of various ranges of env_data.
#' @param temporal_stability_check character string, specifying, how temporal stability
#' between the optimal selection and response variables will be analysed. Current
#' possibilities are "sequential" and "progressive". Sequential check will split data into
#' k splits and calculate selected metric for each split. Progressive check will split data
#' into k splits, calculate metric for the first split and then progressively add 1 split at
#' a time and calculate selected metric.
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
#' columns: "Year", "DOY", "Precipitation/Temperature/etc."
#'
#' @return a list with 13 elements:
#' \tabular{rll}{
#'  1 \tab $calculations   \tab a matrix with calculated metrics\cr
#'  2 \tab $method \tab the character string of a method \cr
#'  3 \tab $metric   \tab the character string indicating the metric used for calculations \cr
#'  4 \tab $analysed_period    \tab the character string specifying the analysed period based on the information from row names. If there are no row names, this argument is given as NA \cr
#'  5 \tab $optimized_return   \tab data frame with two columns, response variable and aggregated (averaged) daily data that return the optimal results. This data.frame could be directly used to calibrate a model for climate reconstruction \cr
#'  6 \tab $optimized_return_all    \tab a data frame with aggregated daily data, that returned the optimal result for the entire env_data (and not only subset of analysed years) \cr
#'  7 \tab $transfer_function    \tab a ggplot object: scatter plot of optimized return and a transfer line of the selected method \cr
#'  8 \tab $temporal_stability    \tab a data frame with calculations of selected metric for different temporal subsets\cr
#'  9\tab $cross_validation   \tab a data frame with cross validation results \cr
#'  10 \tab $plot_heatmap    \tab ggplot2 object: a heatmap of calculated metrics\cr
#'  11 \tab $plot_extreme    \tab ggplot2 object: line plot of a row with the highest value in a matrix of calculated metrics\cr
#'  12 \tab $plot_specific    \tab ggplot2 object: line plot of a row with a selected window width in a matrix of calculated metrics\cr
#'  13 \tab $PCA_output    \tab princomp object: the result output of the PCA analysis
#'}
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' # 1 Example with fixed width
#' example_fixed_width <- daily_response(response = data_MVA, env_data = LJ_daily_temperatures,
#'                                      method = "brnn", fixed_width = 60,
#'                                      row_names_subset = TRUE, remove_insignificant = TRUE,
#'                                      alpha = 0.05)
#' example_fixed_width$plot_extreme
#'
#' # 2 Example for past and present
#' example_MVA_past <- daily_response(response = data_MVA, env_data = LJ_daily_temperatures,
#' method = "cor", lower_limit = 21, upper_limit = 180,
#' row_names_subset = TRUE, previous_year = TRUE,
#' remove_insignificant = TRUE, alpha = 0.05,
#' plot_specific_window = 60, subset_years = c(1940, 1980))
#'
#' example_MVA_present <- daily_response(response = data_MVA, env_data = LJ_daily_temperatures,
#'                                       method = "cor", lower_limit = 21, upper_limit = 180,
#'                                       row_names_subset = TRUE, previous_year = TRUE,
#'                                       remove_insignificant = TRUE, alpha = 0.05,
#'                                       plot_specific_window = 60, subset_years = c(1981, 2010))
#'
#' example_MVA_past$plot_heatmap
#' example_MVA_present$plot_heatmap
#' example_MVA_past$plot_specific
#' example_MVA_present$plot_specific
#'
#' # 3 Example PCA
#' example_PCA <- daily_response(response = example_proxies_individual,
#'                               env_data = LJ_daily_temperatures, method = "lm",
#'                               lower_limit = 21, upper_limit = 180,
#'                               row_names_subset = TRUE, remove_insignificant = TRUE,
#'                               alpha = 0.01, PCA_transformation = TRUE,
#'                               components_selection = "manual", N_components = 2)
#'
#' summary(example_PCA$PCA_output)
#' example_PCA$plot_heatmap
#'
#' # 4 Example negative correlations
#' example_neg_cor <- daily_response(response = data_TRW_1, env_data = LJ_daily_temperatures,
#'                                   method = "cor", lower_limit = 21, upper_limit = 180,
#'                                   row_names_subset = TRUE, remove_insignificant = TRUE,
#'                                   alpha = 0.05)
#'
#' example_neg_cor$plot_heatmap
#' example_neg_cor$plot_extreme
#' example_neg_cor$temporal_stability
#'
#' # 5 Example of multiproxy analysis
#' summary(example_proxies_1)
#' cor(example_proxies_1)
#'
#' example_multiproxy <- daily_response(response = example_proxies_1,
#'                                      env_data = LJ_daily_temperatures,
#'                                      method = "lm", metric = "adj.r.squared",
#'                                      lower_limit = 21, upper_limit = 180,
#'                                      row_names_subset = TRUE, previous_year = FALSE,
#'                                      remove_insignificant = TRUE, alpha = 0.05)
#'
#' example_multiproxy$plot_heatmap
#' }

daily_response <- function(response, env_data, method = "lm",
                           metric = "r.squared", lower_limit = 30,
                           upper_limit = 270, fixed_width = 0,
                           previous_year = FALSE, neurons = 1,
                           brnn_smooth = TRUE, remove_insignificant = TRUE,
                           alpha = .05, row_names_subset = FALSE,
                           PCA_transformation = FALSE, log_preprocess = TRUE,
                           components_selection = 'automatic',
                           eigenvalues_threshold = 1,
                           N_components = 2, use_median = FALSE,
                           temporal_stability_check = "sequential", k = 2,
                           cross_validation_type = "blocked",
                           subset_years = NULL, plot_specific_window = NULL,
                           ylimits = NULL, seed = NULL, tidy_env_data = FALSE) {


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

  # Here we save the original env and response data that will be used later
  response_original <- response
  env_data_original <- env_data


    # For metric calculations, both objects need to have the same length,
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
  # never be analysed, there is no such environmental data avaliable)
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

  # In case of selected window size is less than 14 (2 weeks) or greater than 270 (9 months)
  if (lower_limit < 14) {
    warning("Selected lower_limit is less than 14. Consider increasing it!")
  }

  if (upper_limit > 270) {
    warning("Selected upper_limit is greater than 270. Consider using lower upper_limit!")
  }

  if (fixed_width < 14 & fixed_width > 0) {
    warning("Selected fixed_width is less than 14. Consider increasing it!")
  }

  if (fixed_width > 270) {
    warning("Selected fixed_width is greater than 270. Consider using lower fixed_width!")
  }


  # If PCA_transformation = TRUE, PCA is performed
  if (PCA_transformation == TRUE) {

    # Logarithmic transformation before PCA
    if (log_preprocess == TRUE) {
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
      temporal_matrix <- matrix(NA, nrow = 1,
        ncol = (ncol(env_data) - fixed_width) + 1)

      pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width),
                           style = 3)

      b = 0

      # An iterating loop. In each itteration x is calculated and represents
      # response (dependent) variable. X is a moving average. Window width of
      # a moving window is fixed_width. Next, statistical metric is calculated
      # based on a selected method (cor, lm or brnn). Calculation is stored in
      # temporal matrix.
      for (j in 0: (ncol(env_data) - fixed_width)) {

        b = b + 1

        if (use_median == TRUE){
          x <- apply(env_data[1:nrow(env_data),
                                 (1 + j): (j + fixed_width)],1 , median, na.rm = TRUE)
        } else {
          x <- rowMeans(env_data[1:nrow(env_data),
                                 (1 + j): (j + fixed_width)], na.rm = TRUE)
        }

        # print(paste(j, fixed_width), sep = "")

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)
        temporal_correlation <- cor(response[, 1], x[, 1])

        # Each calculation is printed. Reason: usually it takes several minutes
        # to go through all loops and therefore, users might think that R is
        # not responding. But if each calculation is printed, user could be
        # confident, that R is responding.


        #print (temporal_correlation)
        temporal_matrix[1, j + 1] <- temporal_correlation

        setTxtProgressBar(pb, b)
      }
      close(pb)

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

    pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width),
                         style = 3)

    b = 0

    for (j in 0:(ncol(env_data) - fixed_width)) {

      b = b + 1

      if (use_median == TRUE){
        x <- apply(env_data[1:nrow(env_data),
                               (1 + j) : (j + fixed_width)],1 , median, na.rm = TRUE)
      } else {
        x <- rowMeans(env_data[1:nrow(env_data),
                               (1 + j) : (j + fixed_width)], na.rm = TRUE)
      }

      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      temporal_df <- data.frame(cbind(x, response))
      temporal_model <- lm(x ~ ., data = temporal_df)
      temporal_summary <- summary(temporal_model)
      temporal_r_squared <- temporal_summary$r.squared
      temporal_adj_r_squared <- temporal_summary$adj.r.squared

      if (metric == "r.squared"){
        temporal_matrix[1, j + 1] <- temporal_r_squared
        # print(temporal_r_squared)
      }

      if (metric == "adj.r.squared"){
        temporal_matrix[1, j + 1] <- temporal_adj_r_squared
        # print(temporal_adj_r_squared)
      }


      setTxtProgressBar(pb, b)
    }
    close(pb)

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

    pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width),
                         style = 3)

    b = 0

     for (j in 0: (ncol(env_data) - fixed_width)) {

       b = b + 1

        if (use_median == TRUE){
         x <- apply(env_data[1:nrow(env_data),
                                (1 + j): (j + fixed_width)],1 , median, na.rm = TRUE)
       } else {
         x <- rowMeans(env_data[1:nrow(env_data),
                                (1 + j): (j + fixed_width)], na.rm = TRUE)
       }

      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
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
                                            ncol(as.data.frame(response[, 1]))
                                          -  1))

        if (metric == "r.squared"){
          temporal_matrix[1, j + 1] <- temporal_r_squared
           # print(temporal_r_squared)
        } else if (metric == "adj.r.squared"){
          temporal_matrix[1, j + 1] <- temporal_adj_r_squared
          # print(temporal_adj_r_squared)
        } else {

          temporal_matrix[1, j + 1] <- NA
        }


      }

      setTxtProgressBar(pb, b)
     }
    close(pb)

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
  # variable. x is a moving average, based on rowMeans/apply function.
  # Next, statistical metric is calculated based on a selected method (cor,
  # lm or brnn). Calculation is stored in temporal matrix in a proper place.
  # The position of stored calculation is informative later used for
  # indiciating optimal values.

  pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit),
                       style = 3)

  b = 0


  for (K in lower_limit:upper_limit) {

    b = b + 1

    for (j in 0: (ncol(env_data) - K)) {

      if (use_median == TRUE){
        x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , median, na.rm = TRUE)
      } else {
        x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + K)], na.rm = T)
      }

      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      temporal_correlation <- cor(response[, 1], x[, 1])
      # print(temporal_correlation)
      temporal_matrix[(K - lower_limit) + 1, j + 1] <- temporal_correlation
    }
    setTxtProgressBar(pb, b)
  }

  close(pb)

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


    pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit),
                         style = 3)

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1

      for (j in 0: (ncol(env_data) - K)) {
        if (use_median == TRUE){
          x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , median, na.rm = TRUE)
        } else {
          x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + K)], na.rm = T)
        }

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)
        temporal_df <- data.frame(cbind(x, response))
        temporal_model <- lm(x ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        temporal_r_squared <- temporal_summary$r.squared
        temporal_adj_r_squared <- temporal_summary$adj.r.squared

        if (metric == "r.squared"){
          temporal_matrix[(K - lower_limit) + 1, j + 1]  <-
            temporal_r_squared
          # print(temporal_r_squared)
        }

        if (metric == "adj.r.squared"){
          temporal_matrix[(K - lower_limit) + 1, j + 1]  <-
            temporal_adj_r_squared
          # print(temporal_adj_r_squared)
        }
      }
      setTxtProgressBar(pb, b)
    }

    close(pb)
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

    pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit),
                         style = 3)

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1


      for (j in 0: (ncol(env_data) - K)) {

        if (use_median == TRUE){
          x <- apply(env_data[1:nrow(env_data), (1 + j) : (j + K)],1 , median, na.rm = TRUE)
        } else {
          x <- rowMeans(env_data[1:nrow(env_data), (1 + j) : (j + K)], na.rm = T)
        }

        x <- matrix(x, nrow = nrow(env_data), ncol = 1)
        temporal_df <- data.frame(cbind(x, response))
        capture.output(temporal_model <- try(brnn(x ~ ., data = temporal_df, neurons = neurons,
                                   tol = 1e-6), silent = TRUE))
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

          if (metric == "r.squared"){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <- temporal_r_squared
           # print(temporal_r_squared)
          }

          if (metric == "adj.r.squared"){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-
              temporal_adj_r_squared
          # print(temporal_adj_r_squared)
          }

        } else {
          temporal_matrix[(K - lower_limit) + 1, j + 1] <- NA

        }
      }
      setTxtProgressBar(pb, b)
    }

    close(pb)
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

    if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){
      stop("All calculations are insignificant!")
    }

  # 1 Method is correlation
   if (method == "cor") {
     temporal_matrix[abs(temporal_matrix) < abs(critical_threshold_cor)] <- NA
  # 2 lm and brnn method
      } else if (method == "lm" | method == "brnn") {
    temporal_matrix[abs(temporal_matrix) < abs(critical_threshold_cor2)] <- NA
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

  # The fourth return element is being created: rowMeans/ apply of optimal sequence:
  if (use_median == TRUE){
    dataf <- data.frame(apply(env_data[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1)],1 , median, na.rm = TRUE))
  } else {
    dataf <- data.frame(rowMeans(env_data[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1)],
                                 na.rm = TRUE))
  }

  dataf_full <- cbind(response, dataf)
  colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
  colnames(dataf) <- "Optimized.rowNames"

  ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
  # for the analysed period.

  if (use_median == TRUE){
    dataf_original <- data.frame(apply(env_data_original[, as.numeric(plot_column):
                                         (as.numeric(plot_column) +
                                            as.numeric(row_index) - 1)],1 , median, na.rm = TRUE))
  } else {
    dataf_original <- data.frame(rowMeans(env_data_original[, as.numeric(plot_column):
                                            (as.numeric(plot_column) +
                                               as.numeric(row_index) - 1)],
                                 na.rm = TRUE))
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
    optimized_result <- cor(dataf, response)
  }

  # Just give a nicer colname
  colnames(dataf) <- "Optimized return"

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

analysed_period
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
        empty_list[[m]] <- calculation
        colname = "correlation"
      } else if (method == "lm" & metric == "r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        colname = "r.squared"
        empty_list[[m]] <- summary(MLR)$r.squared
      } else if (method == "lm" & metric == "adj.r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        empty_list[[m]] <- summary(MLR)$adj.r.squared
        colname = "adj.r.squared"
      } else if (method == "brnn" & metric == "r.squared"){
        capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
        if (class(BRNN)[[1]] != "try-error"){
          predictions <- predict(BRNN, dataset_temp, neurons = neurons)
          r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                              sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))
          empty_list[[m]] <- r_squared
          colname = "r.squared"
        } else {
          empty_list[[m]] <- NA
          colname = "r.squared"
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
          colname = "adj.r.squared"
        } else {
          empty_list[[m]] <- NA
          colname = "adj.r.squared"
        }
      }
      }
    m1 <- do.call(rbind, empty_list)
    m2 <- do.call(rbind, empty_list_period)

    temporal_stability <- data.frame(cbind(m2, round(m1, 3)))
    colnames(temporal_stability) <-c("Period", colname)
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
          calculation <- cor(dataset_temp[,1], dataset_temp[,2])
          empty_list[[m]] <- calculation
          colname = "correaltion"
        } else if (method == "lm" & metric == "r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          colname = "r.squared"
          empty_list[[m]] <- summary(MLR)$r.squared
        } else if (method == "lm" & metric == "adj.r.squared"){
          MLR <- lm(optimized_return ~ ., data = dataset_temp)
          empty_list[[m]] <- summary(MLR)$adj.r.squared
          colname = "adj.r.squared"
        } else if (method == "brnn" & metric == "r.squared"){
          capture.output(BRNN <- try(brnn(optimized_return ~ ., data = dataset_temp, neurons = neurons), silent = TRUE))
          if (class(BRNN)[[1]] != "try-error"){
          predictions <- predict(BRNN, dataset_temp, neurons = neurons)
          r_squared <- 1 - (sum((dataset_temp[, 1] - predictions) ^ 2) /
                              sum((dataset_temp[, 1] - mean(dataset_temp[, 1])) ^ 2))
          empty_list[[m]] <- r_squared
          colname = "r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "r.squared"
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
          colname = "adj.r.squared"
          } else {
            empty_list[[m]] <- NA
            colname = "adj.r.squared"
          }
        }
      }
      m1 <- do.call(rbind, empty_list)
      m2 <- do.call(rbind, empty_list_period)

      temporal_stability <- data.frame(cbind(m2, round(m1, 3)))
      colnames(temporal_stability) <-c("Period", colname)
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
                       metric = method, analysed_period = analysed_period,
                       optimized_return = dataf_full,
                       optimized_return_all = dataf_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = cross_validation)
  }


    final_list[[4]]


    plot_heatmapA <- plot_heatmap(final_list)
    plot_extremeA <- plot_extreme(final_list, ylimits = ylimits)

    width_sequence = seq(lower_limit, upper_limit)

    if (is.null(plot_specific_window)){
      (plot_specificA <- "plot_specific_window is not avaliable. No plot_specific is made!")
    } else if (fixed_width != 0){
      plot_specific_window = fixed_width
      plot_specificA <- plot_specific(final_list, window_width = plot_specific_window, ylimits = ylimits)
    } else if (plot_specific_window %in% width_sequence){
      plot_specificA <- plot_specific(final_list, window_width = plot_specific_window, ylimits = ylimits)
    } else (plot_specificA <- "Selected plot_specific_window is not avaliable. No plot_specific is made!")






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
                         PCA_output = PCA_result)
    }

    if (method == "cor"){
      final_list <- list(calculations = temporal_matrix, method = method,
                         metric = method, analysed_period = analysed_period,
                         optimized_return = dataf_full,
                         optimized_return_all = dataf_full_original,
                         transfer_function = p1, temporal_stability = temporal_stability,
                         cross_validation = cross_validation,
                         plot_heatmap = plot_heatmapA,
                         plot_extreme = plot_extremeA,
                         plot_specific = plot_specificA,
                         PCA_output = PCA_result)
    }



  return(final_list)
}
