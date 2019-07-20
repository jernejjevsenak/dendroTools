#' compare_methods
#'
#' Calculates performance metrics for calibration (train) and validation (test)
#' data of different regression methods: multiple linear regression (MLR),
#' artificial neural networks with Bayesian regularization training
#' algorithm (BRNN), (ensemble of) model trees (MT) and random forest of regression
#' trees (RF). With the subset argument, specific methods of interest could be
#' specified. Calculated performance metrics are the correlation coefficient (r),
#' the root mean squared error (RMSE), the root relative squared error (RRSE),
#' the index of agreement (d), the reduction of error (RE), the coefficient of
#' efficiency (CE), the detrended efficiency (DE) and mean bias. For each of the
#' considered methods, there are also residual diagnostic plots available,
#' separately for calibration, holdout and edge data, if applicable.
#'
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param dataset a data frame with dependent and independent variables as
#' columns and (optional) years as row names.
#' @param k number of folds for cross-validation
#' @param repeats number of cross-validation repeats. Should be equal or more
#' than 1
#' @param optimize if set to TRUE (default), the optimal values for the tuning
#' parameters will be selected in a preliminary cross-validation procedure
#' @param dataset_complete optional, a data frame with the full length of tree-ring
#' parameter, which will be used to reconstruct the climate variable specified
#' with the formula argument
#' @param BRNN_neurons number of neurons to be used for the brnn method
#' @param MT_committees an integer:  how many committee models (e.g.  boosting
#' iterations) should be used?
#' @param MT_neighbors how many, if any, neighbors should be used to correct the
#' model predictions
#' @param MT_rules an integer (or NA): define an explicit limit to the number of
#' rules used (NA let’s Cubist decide).
#' @param MT_unbiased a logical: should unbiased rules be used?
#' @param MT_extrapolation a number between 0 and 100: since Cubist uses linear models,
#' predictions can be outside of the outside of the range seen the training set. This
#' parameter controls how much rule predictions are adjusted to be consistent with the
#' training set.
#' @param MT_sample a number between 0 and 99.9:  this is the percentage of the dataset
#' to be randomly selected for model building (not for out-of-bag type evaluation)
#' @param RF_mtry number of variables randomly sampled as candidates at each
#' split
#' @param RF_ntree number of trees to grow. This should not be set to too small
#' a number, to ensure that every input row gets predicted at least a few times
#' @param RF_maxnodes maximum number of terminal nodes trees in the forest can
#' have
#' @param RF_nodesize minimum size of terminal nodes. Setting this number larger
#' causes smaller trees to be grown (and thus take less time).
#' @param seed_factor an integer that will be used to change the seed options
#' for different repeats.
#' @param digits integer of number of digits to be displayed in the final
#' result tables
#' @param edge_share the share of the data to be considered as the edge (extreme) data.
#' This argument could be between 0.10 and 0.50. If the argument is set to 0.10, then
#' the 5 % of the maximal extreme values and 5 % of the minimal extreme values are
#' considered to be the edge data.
#' @param blocked_CV default is FALSE, if changed to TRUE, blocked cross-validation
#' will be used to compare regression methods.
#' @param MLR_stepwise if set to TRUE, stepwise selection of predictors will be used
#' for the MLR method
#' @param stepwise_direction the mode of stepwise search, can be one of "both",
#' "backward", or "forward", with a default of "backward".
#' @param PCA_transformation if set to TRUE, all independent variables will be
#' transformed using PCA transformation.
#' @param log_preprocess if set to TRUE, variables will be transformed with
#' logarithmic transformation before used in PCA
#' @param components_selection character string specifying how to select the Principal
#' Components used as predictors.
#' There are three options: "automatic", "manual" and "plot_selection". If
#' parameter is set to automatic, all scores with eigenvalues above 1 will be
#' selected. This threshold could be changed by changing the
#' eigenvalues_threshold argument. If parameter is set to "manual", user should
#' set the number of components with N_components argument. If component
#' selection is se to "plot_selection", Scree plot will be shown and user must
#' manually enter the number of components used as predictors.
#' @param eigenvalues_threshold threshold for automatic selection of Principal Components
#' @param N_components number of Principal Components used as predictors
#' @param round_bias_cal number of digits for bias in calibration period. Effects
#' the outlook of the final ggplot  of mean bias for calibration data (element 3 of
#' the output list)
#' @param round_bias_val number of digits for bias in validation period. Effects
#' the outlook of the final ggplot of mean bias for validation data (element 4 of
#' the output list)
#' @param n_bins number of bins used for the histograms of mean bias
#' @param methods a vector of strings related to methods that will be compared. A full
#' method vector is methods = c("MLR", "BRNN", "MT", "RF").
#' To use only a subset of methods, pass a vector of methods that you would like to compare.
#' @param tuning_metric a string that specifies what summary metric will be used to select
#' the optimal value of tuning parameters. By default, the argument is set to "RMSE". It is
#' also possible to use "RSquared".
#' @param BRNN_neurons_vector a vector of possible values for BRNN_neurons argument optimization
#' @param MT_committees_vector a vector of possible values for MT_committees argument optimization
#' @param MT_neighbors_vector a vector of possible values for MT_neighbors argument optimization
#' @param MT_rules_vector a vector of possible values for MT_rules argument optimization
#' @param MT_unbiased_vector a vector of possible values for MT_unbiased argument optimization
#' @param MT_extrapolation_vector a vector of possible values for MT_extrapolation argument optimization
#' @param MT_sample_vector a vector of possible values for MT_sample argument optimization
#' @param RF_ntree_vector a vector of possible values for RF_ntree argument optimization
#' @param RF_maxnodes_vector a vector of possible values for RF_maxnodes argument optimization
#' @param RF_mtry_vector a vector of possible values for RF_mtry argument optimization
#' @param RF_nodesize_vector a vector of possible values for RF_nodesize argument optimization
#' @param holdout this argument is used to define observations, which are excluded
#' from the cross-validation and hyperparameters optimization. The holdout argument must be
#' a character with one of the following inputs: “early”, “late” or “manual”. If
#' "early" or "late" characters are specified, then the early or late years will be
#' used as a holdout data. How many of the "early" or "late" years are used as a holdout
#' is specified with the argument holdout_share. If the argument holdout is set to “manual”,
#' then supply a vector of years (or row names) to the argument holdout_manual. Defined
#' years will be used as a holdout. For the holdout data, the same statistical measures are
#' calculated as for the cross-validation. The results for holdout metrics are given in the
#' output element $holdout_results.
#' @param holdout_share the share of the whole dataset to be used as a holdout.
#' Default is 0.10.
#' @param holdout_manual a vector of years (or row names) which will be used as a holdout.
#' calculated as for the cross-validation.
#' @param total_reproducibility logical, default is FALSE. This argument ensures total
#' reproducibility despite the inclusion/exclusion of different methods. By default, the
#' optimization is done only for the methods, that are included in the methods vector. If
#' one method is absent or added, the optimization phase is different, and this affects
#' all the final cross-validation results. By setting the total_reproducibility = TRUE,
#' all methods will be optimized, even though they are not included in the methods vector
#' and the final results will be subset based on the methods vector. Setting the
#' total_reproducibility to TRUE will result in longer optimization phase as well.
#'
#' @return a list with twelve elements:
#'\tabular{rll}{
#'  1 \tab $mean_std   \tab data frame with calculated metrics for the selected regression methods. For each regression method and each calculated metric, mean and standard deviation are given\cr
#'  2 \tab $ranks \tab data frame with ranks of calculated metrics: mean rank and  share of rank_1 are given \cr
#'  3 \tab $edge_results   \tab data frame with calculated performance metrics for the central-edge test. The central part of the data represents the calibration data, while the edge data, i.e. extreme values, represent the test/validation data. Different regression models are calibrated using the central data and validated for the edge (extreme) data. This test is particularly important to assess the performance of models for the predictions of the extreme data. The share of the edge (extreme) data is defined with the edge_share argument \cr
#'  4 \tab $holdout_results    \tab calculated metrics for the holdout data \cr
#'  5 \tab $bias_cal   \tab ggplot object of mean bias for calibration data \cr
#'  6 \tab $bias_val    \tab ggplot object of mean bias for validation data \cr
#'  7 \tab $transfer_functions   \tab ggplot or plotly object with transfer functions of methods \cr
#'  8 \tab $transfer_functions_together    \tab ggplot or plotly object with transfer functions of methods plotted together \cr
#'  9 \tab $parameter_values    \tab a data frame with specifications of parameters used for different regression methods \cr
#'  10 \tab $PCA_output    \tab princomp object: the result output of the PCA analysis \cr
#'  11 \tab $reconstructions    \tab ggplot object: reconstructed dependent variable based on the dataset_complete argument, facet is used to split plots by methods  \cr
#'  12 \tab $reconstructions_together    \tab ggplot object: reconstructed dependent variable based on the dataset_complete argument, all reconstructions are on the same plot \cr
#'  13 \tab $normal_QQ_cal    \tab normal q-q plot for calibration data \cr
#'  14 \tab $normal_QQ_holdout    \tab normal q-q plot for holdout data \cr
#'  15 \tab $normal_QQ_edge    \tab normal q-q plot for edge data \cr
#'  16 \tab $residuals_vs_fitted_cal    \tab residuals vs fitted values plot for calibration data \cr
#'  17 \tab $residuals_vs_fitted_holdout    \tab residuals vs fitted values plot for holdout data \cr
#'  18 \tab $residuals_vs_fitted_edge    \tab residuals vs fitted values plot for edge data
#'}
#'
#' @export
#'
#' @references
#' Bishop, C.M., 1995. Neural Networks for Pattern Recognition. Oxford
#' University Press, Inc. 482 pp.
#'
#' Breiman, L., 1996. Bagging predictors. Machine Learning 24, 123-140.
#'
#' Breiman, L., 2001. Random forests. Machine Learning 45, 5-32.
#'
#' Burden, F., Winkler, D., 2008. Bayesian Regularization of Neural Networks,
#' in: Livingstone, D.J. (ed.), Artificial Neural Networks: Methods and
#' Applications, vol. 458. Humana Press, Totowa, NJ, pp. 23-42.
#'
#' Hastie, T., Tibshirani, R., Friedman, J.H., 2009. The Elements of
#' Statistical Learning : Data Mining, Inference, and Prediction, 2nd ed.
#' Springer, New York xxii, 745 p. pp.
#'
#' Ho, T.K., 1995. Random decision forests, Proceedings of the Third
#' International Conference on Document Analysis and Recognition Volume 1.
#' IEEE Computer Society, pp. 278-282.
#'
#' Hornik, K., Buchta, C., Zeileis, A., 2009. Open-source machine learning: R
#' meets Weka. Comput. Stat. 24, 225-232.
#'
#' Perez-Rodriguez, P., Gianola, D., 2016. Brnn: Brnn (Bayesian Regularization
#' for Feed-forward Neural Networks). R package version 0.6.
#'
#' Quinlan, J.R., 1992. Learning with Continuous Classes, Proceedings of the
#' 5th Australian Joint Conference on Artificial Intelligence (AI '92). World
#' Scientific, Hobart, pp. 343-348.
#'
#' @examples
#' \dontrun{
#'
#' # An example with default settings of machine learning algorithms
#' library(dendroTools)
#' data(example_dataset_1)
#' example_1 <- compare_methods(formula = MVA~., dataset = example_dataset_1,
#' edge_share = 0, holdout = "late")
#' example_1$mean_std
#' example_1$holdout_results
#' example_1$edge_results
#' example_1$ranks
#' example_1$bias_cal
#' example_1$bias_val
#' example_1$transfer_functions
#' example_1$transfer_functions_together
#' example_1$PCA_output
#' example_1$parameter_values
#' example_1$residuals_vs_fitted_cal
#' example_1$residuals_vs_fitted_edge
#' example_1$residuals_vs_fitted_holdout
#' example_1$normal_QQ_cal
#' example_1$normal_QQ_edge
#' example_1$normal_QQ_holdout
#'
#'
#' example_2 <- compare_methods(formula = MVA ~  T_APR,
#' dataset = example_dataset_1, k = 5, repeats = 10, BRNN_neurons = 1,
#' RF_ntree = 100, RF_mtry = 2, RF_maxnodes = 35, seed_factor = 5)
#' example_2$mean_std
#' example_2$ranks
#' example_2$bias_cal
#' example_2$transfer_functions
#' example_2$transfer_functions_together
#' example_2$PCA_output
#' example_2$parameter_values
#'
#' example_3 <- compare_methods(formula = MVA ~ .,
#' dataset = example_dataset_1, k = 2, repeats = 5,
#' methods = c("MLR", "BRNN", "MT"),
#' optimize = TRUE, MLR_stepwise = TRUE)
#' example_3$mean_std
#' example_3$ranks
#' example_3$bias_val
#' example_3$transfer_functions
#' example_3$transfer_functions_together
#' example_3$parameter_values
#'
#' library(dendroTools)
#' library(ggplot2)
#' data(dataset_TRW)
#' comparison_TRW <- compare_methods(formula = T_Jun_Jul ~ TRW, dataset = dataset_TRW,
#' k = 3, repeats = 10, optimize = FALSE, methods = c("MLR", "BRNN", "RF", "MT"),
#' seed_factor = 5, dataset_complete = dataset_TRW_complete, MLR_stepwise = TRUE,
#' stepwise_direction = "backward")
#' comparison_TRW$mean_std
#' comparison_TRW$bias_val
#' comparison_TRW$transfer_functions + xlab(expression(paste('TRW'))) +
#' ylab("June-July Mean Temperature [°C]")
#' comparison_TRW$reconstructions
#' comparison_TRW$reconstructions_together
#' comparison_TRW$edge_results
#' }

compare_methods <- function(formula, dataset, k = 10, repeats = 2,
                            optimize = TRUE, dataset_complete = NULL,
                            BRNN_neurons = 1, MT_committees = 1, MT_neighbors = 5,
                            MT_rules = 200, MT_unbiased = TRUE, MT_extrapolation = 100,
                            MT_sample = 0,
                            RF_ntree = 500,
                            RF_maxnodes = 5,
                            RF_mtry  = 1,
                            RF_nodesize = 1,
                            seed_factor = 5,
                            digits = 3, blocked_CV = FALSE,
                            PCA_transformation = FALSE, log_preprocess = TRUE,
                            components_selection = 'automatic',
                            eigenvalues_threshold = 1, N_components = 2,
                            round_bias_cal = 15, round_bias_val = 4,
                            n_bins = 30, edge_share = 0.10,
                            MLR_stepwise = FALSE, stepwise_direction = "backward",
                            methods = c("MLR", "BRNN", "MT", "RF"),
                            tuning_metric = "RMSE",
                            BRNN_neurons_vector = c(1, 2, 3),
                            MT_committees_vector = c(1, 5, 10),
                            MT_neighbors_vector = c(0, 5),
                            MT_rules_vector = c(100, 200),
                            MT_unbiased_vector = c(TRUE, FALSE),
                            MT_extrapolation_vector = c(100),
                            MT_sample_vector = c(0),
                            RF_ntree_vector = c(100, 250, 500),
                            RF_maxnodes_vector = c(5, 10, 20, 25),
                            RF_mtry_vector  = c(1),
                            RF_nodesize_vector = c(1, 5, 10),
                            holdout = NULL, holdout_share = 0.10,
                            holdout_manual = NULL, total_reproducibility = FALSE){

if (k < 2 | k > 26){
  stop(paste0("Selected k is ", k,", but it should be between 2 and 26"))
}

if (repeats > 1 & blocked_CV == TRUE){
  warning("blocked_CV is set to TRUE, repeats argument is ignored")
}

dataset <- data.frame(dataset) # dataset needs to be of class data.frame!

full_methods <- c("MLR", "BRNN", "MT", "RF")

if ("BMT" %in% methods){
  warning(paste0("BMT and MT from RWeka package are now replaced with cubist R function and refered as MT.",
                 "To fit ensemble model tree, use committees argument (MT_committees)."))
}

methods <- sort(methods)

set.seed(seed_factor) # We ensure that the optimization results are always the same

#############################################################################

# Here, empty lists are defined, where calculations will be stored. Empty lists
# for bias are defined separately, since bias should not be averaged. It is
# later given as density plots
list_MLR <- list()
list_BRNN <- list()
list_MT <- list()
list_RF <- list()

# Here, idex of dependent variable is extracted and later used to locate the
# observed values
DepIndex <- grep(as.character(formula[[2]]), colnames(dataset))
DepName <- as.character(formula[[2]])

# If PCA_transformation = TRUE, PCA is performed
if (PCA_transformation == TRUE) {

  # Logarithmic transformation before PCA
  if (log_preprocess == TRUE) {

    dataset_temp <- dataset[,-DepIndex]
    dataset_temp <- data.frame(log(dataset_temp))
  }

  PCA_result <- princomp(dataset_temp, cor = TRUE)

  if (components_selection == 'automatic'){
    subset_vector <- PCA_result$sdev > eigenvalues_threshold
    dataset_temp <- as.data.frame(PCA_result$scores[, subset_vector])
  }

  if (components_selection == 'manual'){
    dataset_temp <- as.data.frame(PCA_result$scores[, 1:N_components])
  }

  if (components_selection == 'plot_selection'){
    plot(PCA_result, type = 'l')

    fun <- function(){
      N_PC <- readline("What number of PC scores should be used as predictors? ")
      return(N_PC)
    }

    N_PC <- fun()
    dataset_temp <- as.data.frame(PCA_result$scores[, 1:as.numeric(N_PC)])
  }

  dataset <- data.frame(dataset[, DepIndex], dataset_temp)
  colnames(dataset)[1] = DepName
  for (i in 2:ncol(dataset)){
    colnames(dataset)[i] <- paste("PC_", i-1, sep = "")
    formula = as.formula(paste(DepName, "~ ."))
  }
} else (PCA_result <- "No PCA result avalialbe!")

# Here we fit a lm model, just to get information about the number of independent variables;
# when formula is used in the form of: y~., we don't know the number of independent variables
# this information is used later
DepIndex <- grep(as.character(formula[[2]]), colnames(dataset))
quazi_mod <- lm(formula, data = dataset)
numIND <- length(quazi_mod[[1]]) - 1
indep_names <- colnames(quazi_mod[[7]][[1]])[-1]
indep_name_1 <- colnames(quazi_mod[[7]][[1]])[-1][1]
indep_name_2 <- colnames(quazi_mod[[7]][[1]])[-1][2]
allnames <- c(as.character(formula[[2]]), indep_names)

# So, sometimes there are variables in dataframe, which are not used as predictors.
# Let's remove them
dataset <- dataset[ names(dataset)[names(dataset) %in% allnames] ]


# Here, mtry is checked. mtry should not be larger than the unmber of independet variables

if(sum(RF_mtry_vector > numIND) > 0){
critical_values <- as.numeric(which(RF_mtry_vector > numIND))
critical_values_inf <- RF_mtry_vector[critical_values]
RF_mtry_vector <- RF_mtry_vector[-critical_values]
warning(paste0("mtry argument for the RF should not be larger than the number of independent variables."),
        " Values ", paste((critical_values_inf), collapse=", "), " for the argument RF_mtry_vector are depreciated!")
}

# Here, we define stop, if edge_share is 0
if (!edge_share > 0 & numIND < 2){
  stop("edge_share should be greater than 0!")
}

# Here, if holdout data is defined, we exclude those observations from cross-validation
# and hyperparameters optimization
if (!is.null(holdout)){
  if (holdout == "early"){
  dataset <- dataset[ order(as.numeric(row.names(dataset)), decreasing = TRUE),]
  dataset_holdout = dataset[1:round((nrow(dataset)*holdout_share),0),]
  dataset = dataset[round(((nrow(dataset)*holdout_share)+1),0):nrow(dataset),]

  } else if (holdout == "late"){
    dataset <- dataset[ order(as.numeric(row.names(dataset)), decreasing = TRUE),]
    dataset_holdout = dataset[(nrow(dataset) - (round2(nrow(dataset)*holdout_share, 0)) + 1):nrow(dataset), ]
    dataset = dataset[1:((nrow(dataset) - (round2(nrow(dataset)*holdout_share, 0)))), ]

  } else if (holdout == "manual"){
    dataset_holdout = dataset[row.names(dataset) %in% holdout_manual, ]
    dataset = dataset[!row.names(dataset) %in% holdout_manual, ]

  } else {
      stop(paste0("The argument holdout is set to ", holdout,". Instead, it should be early, late or manual. ",
                 "Please, see the compare_methods() documentation."))
    }
  }


# Here we use caret package to tune parameters of different methods
# Optimization for randomized procedure

if (optimize == TRUE & blocked_CV == FALSE){

  print("Optimization phase...")

  model = NULL

  if ("BRNN" %in% methods | total_reproducibility == TRUE){

    print("Tuning parameters for the BRNN method...")
    # 2 Optimization for BRNN
    BRNN_neurons <- BRNN_neurons_vector
    hyper_grid <- expand.grid(neurons = BRNN_neurons)

    # Number of potential models in the grid
    num_models <- nrow(hyper_grid)

    # Create an empty list to store models
    grade_models <- list()

    # Write a loop over the rows of hyper_grid to train the grid of models
    for (i in 1:num_models) {

      # Get minsplit, maxdepth values at row i
      neurons <- hyper_grid$neurons[i]

      # cross_validation
      foldi <- seq(1:k)
      foldi <- paste("fold_", foldi)

      #Randomly shuffle the data
      dataset <- dataset[sample(nrow(dataset)), ]

      #Create 10 equally size folds
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      #Perform k fold cross validation
      tuning_vector <- c()

      for (j in 1:k){

        #Segement your data by fold using the which() function
        testIndexes <- which(folds == j, arr.ind = TRUE)
        test <- dataset[testIndexes, ]
        train <- dataset[-testIndexes, ]

        #MLR MODEL
        capture.output(model_temp <- brnn(formula, data = train, neurons = neurons, verbose = FALSE,
                                          tol = 1e-6))
        test_observed <- test[, DepIndex]
        test_predicted <- predict(model_temp, test)

        if (tuning_metric == "RMSE"){
          tuning_vector[j] <- MLmetrics::RMSE(test_predicted, test_observed)
        } else if (tuning_metric == "RSquared"){
          tuning_vector[j] <- cor(test_predicted, test_observed)^2
        } else {
          stop(paste0("tuning_metric argument should be RMSE or RSquared! Instead it is ", tuning_metric))
        }

      }

      grade_models[i] <- mean(tuning_vector)

    }

    grade_list <- unlist(grade_models)

    # Identify the model with smallest validation set RMSE
    if (tuning_metric == "RMSE"){
        best_model <- which.min(grade_list)}
    if (tuning_metric == "RSquared"){
      best_model <- which.max(grade_list)
      }

    best_parameters <- hyper_grid[best_model, ]

    BRNN_neurons <- as.numeric(best_parameters[1])
  }

  if ("MT" %in% methods | total_reproducibility == TRUE){

  print("Tuning parameters for the MT method...")

  # 2 Optimization for MT
    MT_committees <- MT_committees_vector
    MT_neighbors <- MT_neighbors_vector
    MT_rules <- MT_rules_vector
    MT_unbiased <- MT_unbiased_vector
    MT_extrapolation <- MT_extrapolation_vector
    MT_sample <- MT_sample_vector

    hyper_grid <- expand.grid(
    committees = MT_committees,
    neighbors = MT_neighbors,
    rules = MT_rules,
    unbiased = MT_unbiased,
    extrapolation = MT_extrapolation,
    sample = MT_sample)

  # Number of potential models in the grid
  num_models <- nrow(hyper_grid)

  # Create an empty list to store models
  grade_models <- list()

  # Write a loop over the rows of hyper_grid to train the grid of models
  for (i in 1:num_models) {

    # Get minsplit, maxdepth values at row i
    committees <- hyper_grid$committees[i]
    neighbors <- hyper_grid$neighbors[i]
    rules <- hyper_grid$rules[i]
    unbiased <- hyper_grid$unbiased[i]
    extrapolation <- hyper_grid$extrapolation[i]
    sample <- hyper_grid$sample[i]

    # cross_validation
    foldi <- seq(1:k)
    foldi <- paste("fold_", foldi)

    #Randomly shuffle the data
    dataset <- dataset[sample(nrow(dataset)), ]

    #Create 10 equally size folds
    folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

    #Perform k fold cross validation
    tuning_vector <- c()

    for (j in 1:k){

      #Segement your data by fold using the which() function
      testIndexes <- which(folds == j, arr.ind = TRUE)
      test <- dataset[testIndexes, ]
      train <- dataset[-testIndexes, ]

      #MT MODEL
      train_ind <- data.frame(train[, -DepIndex])
      colnames(train_ind) <- indep_names
      test_ind <- data.frame(test[, -DepIndex])
      colnames(test_ind) <- indep_names

      train_dep <- train[, DepIndex]
      test_dep <- test[, DepIndex]

      model_temp <- cubist(x = train_ind, y = train_dep, committees = committees,
                           neighbors = neighbors, cubistControl(rules = rules, unbiased = unbiased,
                                                                extrapolation = extrapolation,
                                                                sample = sample))

      test_observed <- test[, DepIndex]
      test_predicted <- predict(model_temp, test_ind)

      if (tuning_metric == "RMSE"){
        tuning_vector[j] <- MLmetrics::RMSE(test_predicted, test_observed)
      } else if (tuning_metric == "RSquared"){
        tuning_vector[j] <- cor(test_predicted, test_observed)^2
      } else {
        stop(paste0("tuning_metric argument should be RMSE or RSquared! Instead it is ", tuning_metric))
      }

    }

    grade_models[i] <- mean(tuning_vector)

  }

  grade_list <- unlist(grade_models)

  # Identify the model with smallest validation set RMSE
  if (tuning_metric == "RMSE"){
    best_model <- which.min(grade_list)}
  if (tuning_metric == "RSquared"){
    best_model <- which.max(grade_list)
  }

  best_parameters <- hyper_grid[best_model, ]

  MT_committees <- as.numeric(best_parameters[1])
  MT_neighbors <- as.numeric(best_parameters[2])
  MT_rules <- as.numeric(best_parameters[3])
  MT_unbiased <- as.logical(best_parameters[4])
  MT_extrapolation <- as.numeric(best_parameters[5])
  MT_sample <- as.numeric(best_parameters[6])

  }


  if ("RF" %in% methods | total_reproducibility == TRUE){
    # 2 Optimization for RF

    print("Tuning parameters for the RF method...")

    RF_mtry <- RF_mtry_vector
    RF_maxnodes <- RF_maxnodes_vector
    RF_ntree  <- RF_ntree_vector
    RF_nodesize  <- RF_nodesize_vector

    hyper_grid <- expand.grid(mtry = RF_mtry, maxnodes = RF_maxnodes,
                              ntree = RF_ntree,  nodesize = RF_nodesize)

    # Number of potential models in the grid
    num_models <- nrow(hyper_grid)

    # Create an empty list to store models
    grade_models <- list()

    # Write a loop over the rows of hyper_grid to train the grid of models
    for (i in 1:num_models) {

      # Get minsplit, maxdepth values at row i
      mtry <- hyper_grid$mtry[i]
      maxnodes <- hyper_grid$maxnodes[i]
      ntree <- hyper_grid$ntree[i]
      nodesize <- hyper_grid$nodesize[i]

      # cross_validation
      foldi <- seq(1:k)
      foldi <- paste("fold_", foldi)

      #Randomly shuffle the data
      dataset <- dataset[sample(nrow(dataset)), ]

      #Create 10 equally size folds
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      #Perform k fold cross validation
      tuning_vector <- c()

      for (j in 1:k){

        #Segement your data by fold using the which() function
        testIndexes <- which(folds == j, arr.ind = TRUE)
        test <- dataset[testIndexes, ]
        train <- dataset[-testIndexes, ]

        #MLR MODEL


        model_temp <- randomForest(formula, data = train, mtry = mtry, maxnodes = maxnodes,
                                   ntree = ntree, nodesize = nodesize)

        test_observed <- test[, DepIndex]
        test_predicted <- predict(model_temp, test)

        if (tuning_metric == "RMSE"){
          tuning_vector[j] <- MLmetrics::RMSE(test_predicted, test_observed)
        } else if (tuning_metric == "RSquared"){
          tuning_vector[j] <- cor(test_predicted, test_observed)^2
        } else {
          stop(paste0("tuning_metric argument should be RMSE or RSquared! Instead it is ", tuning_metric))
        }

      }

      grade_models[i] <- mean(tuning_vector)

    }

    grade_list <- unlist(grade_models)

    # Identify the model with smallest validation set RMSE
    if (tuning_metric == "RMSE"){
      best_model <- which.min(grade_list)}
    if (tuning_metric == "RSquared"){
      best_model <- which.max(grade_list)
    }

    best_parameters <- hyper_grid[best_model, ]

    RF_mtry <- as.numeric(best_parameters[1])
    RF_maxnodes <- as.numeric(best_parameters[2])
    RF_ntree  <- as.numeric(best_parameters[3])
    RF_nodesize  <- as.numeric(best_parameters[4])

  }

}


###############################################################################
# Optimization for blocked_CV


if (optimize == TRUE & blocked_CV == TRUE){

  print("Optimization phase...")

  model = NULL

  if ("BRNN" %in% methods | total_reproducibility == TRUE){

    print("Tuning parameters for the BRNN method...")
    # 2 Optimization for BRNN
    BRNN_neurons <- BRNN_neurons_vector
    hyper_grid <- expand.grid(neurons = BRNN_neurons)

    # Number of potential models in the grid
    num_models <- nrow(hyper_grid)

    # Create an empty list to store models
    grade_models <- list()

    # Write a loop over the rows of hyper_grid to train the grid of models
    for (i in 1:num_models) {

      # Get minsplit, maxdepth values at row i
      neurons <- hyper_grid$neurons[i]

      # cross_validation
      foldi <- seq(1:k)
      foldi <- paste("fold_", foldi)

      #Randomly shuffle the data
      # dataset <- dataset[sample(nrow(dataset)), ]

      #Create 10 equally size folds
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      #Perform k fold cross validation
      tuning_vector <- c()

      for (j in 1:k){

        #Segement your data by fold using the which() function
        testIndexes <- which(folds == j, arr.ind = TRUE)
        test <- dataset[testIndexes, ]
        train <- dataset[-testIndexes, ]

        #MLR MODEL
        capture.output(model_temp <- brnn(formula, data = train, neurons = neurons, verbose = FALSE))
        test_observed <- test[, DepIndex]
        test_predicted <- predict(model_temp, test)

        if (tuning_metric == "RMSE"){
          tuning_vector[j] <- MLmetrics::RMSE(test_predicted, test_observed)
        } else if (tuning_metric == "RSquared"){
          tuning_vector[j] <- cor(test_predicted, test_observed)^2
        } else {
          stop(paste0("tuning_metric argument should be RMSE or RSquared! Instead it is ", tuning_metric))
        }

      }

      grade_models[i] <- mean(tuning_vector)

    }

    grade_list <- unlist(grade_models)

    # Identify the model with smallest validation set RMSE
    if (tuning_metric == "RMSE"){
      best_model <- which.min(grade_list)}
    if (tuning_metric == "RSquared"){
      best_model <- which.max(grade_list)
    }

    best_parameters <- hyper_grid[best_model, ]

    BRNN_neurons <- as.numeric(best_parameters[1])
  }




  if ("MT" %in% methods | total_reproducibility == TRUE){

    print("Tuning parameters for the MT method...")

    # 2 Optimization for MT
    MT_committees <- MT_committees_vector
    MT_neighbors <- MT_neighbors_vector
    MT_rules <- MT_rules_vector
    MT_unbiased <- MT_unbiased_vector
    MT_extrapolation <- MT_extrapolation_vector
    MT_sample <- MT_sample_vector

    hyper_grid <- expand.grid(
      committees = MT_committees,
      neighbors = MT_neighbors,
      rules = MT_rules,
      unbiased = MT_unbiased,
      extrapolation = MT_extrapolation,
      sample = MT_sample)

    # Number of potential models in the grid
    num_models <- nrow(hyper_grid)

    # Create an empty list to store models
    grade_models <- list()

    # Write a loop over the rows of hyper_grid to train the grid of models
    for (i in 1:num_models) {

      # Get minsplit, maxdepth values at row i
      committees <- hyper_grid$committees[i]
      neighbors <- hyper_grid$neighbors[i]
      rules <- hyper_grid$rules[i]
      unbiased <- hyper_grid$unbiased[i]
      extrapolation <- hyper_grid$extrapolation[i]
      sample <- hyper_grid$sample[i]

      # cross_validation
      foldi <- seq(1:k)
      foldi <- paste("fold_", foldi)

      #Randomly shuffle the data
      #dataset <- dataset[sample(nrow(dataset)), ]

      #Create 10 equally size folds
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      #Perform k fold cross validation
      tuning_vector <- c()

      for (j in 1:k){

        #Segement your data by fold using the which() function
        testIndexes <- which(folds == j, arr.ind = TRUE)
        test <- dataset[testIndexes, ]
        train <- dataset[-testIndexes, ]

        train_ind <- data.frame(train[, -DepIndex])
        colnames(train_ind) <- indep_names
        test_ind <- data.frame(test[, -DepIndex])
        colnames(test_ind) <- indep_names

        train_dep <- train[, DepIndex]
        test_dep <- test[, DepIndex]

        #MT MODEL
        model_temp <- cubist(x = train_ind, y = train_dep, committees = committees,
                             neighbors = neighbors, cubistControl(rules = rules, unbiased = unbiased,
                                                                  extrapolation = extrapolation,
                                                                  sample = sample))

        test_observed <- test[, DepIndex]
        test_predicted <- predict(model_temp, test_ind)

        if (tuning_metric == "RMSE"){
          tuning_vector[j] <- MLmetrics::RMSE(test_predicted, test_observed)
        } else if (tuning_metric == "RSquared"){
          tuning_vector[j] <- cor(test_predicted, test_observed)^2
        } else {
          stop(paste0("tuning_metric argument should be RMSE or RSquared! Instead it is ", tuning_metric))
        }

      }

      grade_models[i] <- mean(tuning_vector)

    }

    grade_list <- unlist(grade_models)

    # Identify the model with smallest validation set RMSE
    if (tuning_metric == "RMSE"){
      best_model <- which.min(grade_list)}
    if (tuning_metric == "RSquared"){
      best_model <- which.max(grade_list)
    }

    best_parameters <- hyper_grid[best_model, ]

    MT_committees <- as.numeric(best_parameters[1])
    MT_neighbors <- as.numeric(best_parameters[2])
    MT_rules <- as.numeric(best_parameters[3])
    MT_unbiased <- as.logical(best_parameters[4])
    MT_extrapolation <- as.numeric(best_parameters[5])
    MT_sample <- as.numeric(best_parameters[6])

  }

  if ("RF" %in% methods | total_reproducibility == TRUE){
    # 2 Optimization for RF

    print("Tuning parameters for the RF method...")

    RF_mtry <- RF_mtry_vector
    RF_maxnodes <- RF_maxnodes_vector
    RF_ntree  <- RF_ntree_vector
    RF_nodesize  <- RF_nodesize_vector

    hyper_grid <- expand.grid(mtry = RF_mtry, maxnodes = RF_maxnodes,
                              ntree = RF_ntree,  nodesize = RF_nodesize)

    # Number of potential models in the grid
    num_models <- nrow(hyper_grid)

    # Create an empty list to store models
    grade_models <- list()

    # Write a loop over the rows of hyper_grid to train the grid of models
    for (i in 1:num_models) {

      # Get minsplit, maxdepth values at row i
      mtry <- hyper_grid$mtry[i]
      maxnodes <- hyper_grid$maxnodes[i]
      ntree <- hyper_grid$ntree[i]
      nodesize <- hyper_grid$nodesize[i]

      # cross_validation
      foldi <- seq(1:k)
      foldi <- paste("fold_", foldi)

      #Randomly shuffle the data
      #dataset <- dataset[sample(nrow(dataset)), ]

      #Create 10 equally size folds
      folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

      #Perform k fold cross validation
      tuning_vector <- c()

      for (j in 1:k){

        #Segement your data by fold using the which() function
        testIndexes <- which(folds == j, arr.ind = TRUE)
        test <- dataset[testIndexes, ]
        train <- dataset[-testIndexes, ]

        #MLR MODEL
        model_temp <- randomForest(formula, data = train, mtry = mtry, maxnodes = maxnodes,
                                   ntree = ntree, nodesize = nodesize)

        test_observed <- test[, DepIndex]
        test_predicted <- predict(model_temp, test)

        if (tuning_metric == "RMSE"){
          tuning_vector[j] <- MLmetrics::RMSE(test_predicted, test_observed)
        } else if (tuning_metric == "RSquared"){
          tuning_vector[j] <- cor(test_predicted, test_observed)^2
        } else {
          stop(paste0("tuning_metric argument should be RMSE or RSquared! Instead it is ", tuning_metric))
        }

      }

      grade_models[i] <- mean(tuning_vector)

    }

    grade_list <- unlist(grade_models)

    # Identify the model with smallest validation set RMSE
    if (tuning_metric == "RMSE"){
      best_model <- which.min(grade_list)}
    if (tuning_metric == "RSquared"){
      best_model <- which.max(grade_list)
    }

    best_parameters <- hyper_grid[best_model, ]

    RF_mtry <- as.numeric(best_parameters[1])
    RF_maxnodes <- as.numeric(best_parameters[2])
    RF_ntree  <- as.numeric(best_parameters[3])
    RF_nodesize  <- as.numeric(best_parameters[4])
  }

}


##################################################################################
##################################################################################
##################################################################################
print("Evaluation of methods...")

# Normal cross-validation with repeats.

if (blocked_CV == FALSE){

# create progress bar
pb <- txtProgressBar(min = 0, max = repeats, style = 3)

b = 0 # place holder for saving results

for (m in 1:repeats){

foldi <- seq(1:k)
foldi <- paste("fold_", foldi)

#ly shuffle the data
set.seed(seed_factor * m)
dataset <- dataset[sample(nrow(dataset)), ]

#Create 10 equally size folds
folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)


#Perform k fold cross validation

for (j in 1:k){

  b <- b + 1

  #Segement your data by fold using the which() function
  testIndexes <- which(folds == j, arr.ind = TRUE)
  test <- dataset[testIndexes, ]
  train <- dataset[-testIndexes, ]

  #MLR MODEL
  if (MLR_stepwise == FALSE) {
    MLR <- lm(formula, data = train)
    train_predicted <- predict(MLR, train)
    test_predicted <- predict(MLR, test)
    train_observed <- train[, DepIndex]
    test_observed <- test[, DepIndex]
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                      train_observed, test_observed, digits = 15,
                                      formula = formula, test = test)
    list_MLR[[b]] <- calculations

  } else {
    capture.output(MLR <- step(lm(formula = formula, data = train), direction = stepwise_direction))
    train_predicted <- predict(MLR, train)
    test_predicted <- predict(MLR, test)
    train_observed <- train[, DepIndex]
    test_observed <- test[, DepIndex]
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                      train_observed, test_observed, digits = 15,
                                      formula = formula, test = test)
    list_MLR[[b]] <- calculations
  }


  #BRNN Model
  capture.output(BRNN <- brnn(formula, data = train, neurons = BRNN_neurons, verbose = FALSE,
                              tol = 1e-6))
  train_predicted <- predict(BRNN, train)
  test_predicted <- predict(BRNN, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15,
                                    formula = formula, test = test)
  list_BRNN[[b]] <- calculations

  # Model Trees
  train_ind <- data.frame(train[, -DepIndex])
  colnames(train_ind) <- indep_names
  test_ind <- data.frame(test[, -DepIndex])
  colnames(test_ind) <- indep_names

  train_dep <- train[, DepIndex]
  test_dep <- test[, DepIndex]

  MT_model <- cubist(x = train_ind, y = train_dep, committees = MT_committees,
                     neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                      unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
  train_predicted <- predict(MT_model, train_ind)
  test_predicted <- predict(MT_model, test_ind)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                    train_observed, test_observed, digits = 15,
                                    formula = formula, test = test)
  list_MT[[b]] <- calculations

  # Random Forest, Regression Tree with random forest, WEKA
  RF_model <- randomForest(formula = formula, data = train, ntree = RF_ntree,
                           maxnodes = RF_maxnodes, mtry  = RF_mtry,
                           nodesize = RF_nodesize)
  train_predicted <- predict(RF_model, train)
  test_predicted <- predict(RF_model, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15,
                                    formula = formula, test = test)
  list_RF[[b]] <- calculations
}
  setTxtProgressBar(pb, m)

  } # repeats zaključek

close(pb)

position <- k * repeats

}




###################################################################################
##### And now the second option: Blocked cross-validation #########################
###################################################################################

if (blocked_CV == TRUE){

  # create progress bar
  pb <- txtProgressBar(min = 0, max = k, style = 3)

  b = 0 # place holder for saving results

  # Here, idex of dependent variable is extracted and later used to locate the
  # observed values
  DepIndex <- grep(as.character(formula[[2]]), colnames(dataset))

  foldi <- seq(1:k)
  #foldi <- paste("fold_", foldi)
  folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

  #Perform k fold cross validation

  for (m in 1:k){

    b <- b + 1
    #Segement your data by fold using the which() function
    testIndexes <- which(folds == m, arr.ind = TRUE)
    test <- dataset[testIndexes, ]
    train <- dataset[-testIndexes, ]

    #MLR MODEL
    if (MLR_stepwise == FALSE) {
      MLR <- lm(formula, data = train)
      train_predicted <- predict(MLR, train)
      test_predicted <- predict(MLR, test)
      train_observed <- train[, DepIndex]
      test_observed <- test[, DepIndex]
      calculations <- calculate_metrics(train_predicted, test_predicted,
                                        train_observed, test_observed, digits = 15,
                                        formula = formula, test = test)
      list_MLR[[b]] <- calculations

    } else {
      capture.output(MLR <- step(lm(formula = formula, data = train), direction = stepwise_direction))
      train_predicted <- predict(MLR, train)
      test_predicted <- predict(MLR, test)
      train_observed <- train[, DepIndex]
      test_observed <- test[, DepIndex]
      calculations <- calculate_metrics(train_predicted, test_predicted,
                                        train_observed, test_observed, digits = 15,
                                        formula = formula, test = test)
      list_MLR[[b]] <- calculations
    }

    #BRNN Model
    capture.output(BRNN <- brnn(formula, data = train, neurons = BRNN_neurons, verbose = FALSE,
                                tol = 1e-6))
    train_predicted <- predict(BRNN, train)
    test_predicted <- predict(BRNN, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15,
                                      formula = formula, test = test)
    list_BRNN[[b]] <- calculations

    # Model Trees
    train_ind <- data.frame(train[, -DepIndex])
    colnames(train_ind) <- indep_names
    test_ind <- data.frame(test[, -DepIndex])
    colnames(test_ind) <- indep_names

    train_dep <- train[, DepIndex]
    test_dep <- test[, DepIndex]

    MT_model <- cubist(x = train_ind, y = train_dep, committees = MT_committees,
                       neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
    train_predicted <- predict(MT_model, train_ind)
    test_predicted <- predict(MT_model, test_ind)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                      train_observed, test_observed, digits = 15,
                                      formula = formula, test = test)
    list_MT[[b]] <- calculations


    ##Random Forest
    RF_model <- randomForest(formula = formula, data = train, ntree = RF_ntree,
                             maxnodes = RF_maxnodes, mtry  = RF_mtry,
                             nodesize = RF_nodesize)
    train_predicted <- predict(RF_model, train)
    test_predicted <- predict(RF_model, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                      train_observed, test_observed, digits = 15,
                                      formula = formula, test = test)
    list_RF[[b]] <- calculations

    setTxtProgressBar(pb, m)
  }

  close(pb)

  position <- k
}

###########################################################################################
###########################################################################################
###########################################################################################
# Now the proces of extraction starts
# Here, lists are rearranged and metrics are extracted

listVec <- lapply(list_MLR, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean, na.rm = TRUE)
std <- apply(m, 1, sd, na.rm = TRUE)
m <- cbind(m, averages, std)
df_MLR <- data.frame(m)
df_MLR_bias <- df_MLR[c(15, 16), c(1: position)]
df_MLR_rank <- df_MLR[-c(15, 16), c(1: position)]
df_MLR_avg <- df_MLR[-c(15, 16), c(position + 1, position + 2)]
rownames(df_MLR_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RRSE_cal",
                      "RRSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val", "DE_cal", "DE_val")

listVec <- lapply(list_BRNN, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean, na.rm = TRUE)
std <- apply(m, 1, sd, na.rm = TRUE)
m <- cbind(m, averages, std)
df_BRNN <- data.frame(m)
df_BRNN_bias <- df_BRNN[c(15, 16), c(1: position)]
df_BRNN_rank <- df_BRNN[-c(15, 16), c(1: position)]
df_BRNN_avg <- df_BRNN[-c(15, 16), c(position + 1, position + 2)]
rownames(df_BRNN_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RRSE_cal",
                      "RRSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val", "DE_cal", "DE_val")

listVec <- lapply(list_MT, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean, na.rm = TRUE)
std <- apply(m, 1, sd, na.rm = TRUE)
m <- cbind(m, averages, std)
df_MT <- data.frame(m)
df_MT_bias <- df_MT[c(15, 16), c(1: position)]
df_MT_rank <- df_MT[-c(15, 16), c(1: position)]
df_MT_avg <- df_MT[-c(15, 16), c(position + 1, position + 2)]
rownames(df_MT_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RRSE_cal",
                      "RRSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val", "DE_cal", "DE_val")


listVec <- lapply(list_RF, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean, na.rm = TRUE)
std <- apply(m, 1, sd, na.rm = TRUE)
m <- cbind(m, averages, std)
df_RF <- data.frame(m)
df_RF_bias <- df_RF[c(15, 16), c(1: position)]
df_RF_rank <- df_RF[-c(15, 16), c(1: position)]
df_RF_avg <- df_RF[-c(15, 16), c(position + 1, position + 2)]
rownames(df_RF_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RRSE_cal",
                      "RRSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val", "DE_cal", "DE_val")


##################################################################################################
##################################################################################################
############################## Here, all data frames are binded together #########################
##################################################################################################

start_position = 0
method_list <- list()
for (i in methods){
  start_position <- start_position + 1
  method_list[[start_position]] <- paste0("df_",i,"_avg")
}

method_vector <- unlist(method_list, use.names=FALSE)

empty_LIST <- list()

for (i in 1:length(method_vector)){
  temp_DF <- get(method_vector[i])
  empty_LIST[[i]] <- temp_DF
}

names(empty_LIST) <- method_vector

df_all_avg <- round(do.call(cbind, empty_LIST), 8)


############################################################################################
############################################################################################
# Calculation of ranks

# First, let's subset the methods, similar to the previous

start_position = 0
method_list <- list()
for (i in methods){
  start_position <- start_position + 1
  method_list[[start_position]] <- paste0("df_",i,"_rank")
}

method_vector <- unlist(method_list, use.names = FALSE)

empty_LIST <- list()

for (i in 1:length(method_vector)){
  temp_DF <- get(method_vector[i])
  empty_LIST[[i]] <- temp_DF
}

names(empty_LIST) <- method_vector

df_all <- round(do.call(rbind, empty_LIST), 8)

# Now, all metrics (except bias) are extracted for calibration and validation
# data.
r_cal <- df_all[c(seq(1, nrow(df_all), by = 14)), ]
r_val <- df_all[c(seq(2, nrow(df_all), by = 14)), ]

RMSE_cal <- df_all[c(seq(3, nrow(df_all), by = 14)), ]
RMSE_val <- df_all[c(seq(4, nrow(df_all), by = 14)), ]

RRSE_cal <- df_all[c(seq(5, nrow(df_all), by = 14)), ]
RRSE_val <- df_all[c(seq(6, nrow(df_all), by = 14)), ]

d_cal <- df_all[c(seq(7, nrow(df_all), by = 14)), ]
d_val <- df_all[c(seq(8, nrow(df_all), by = 14)), ]

RE_cal <- df_all[c(seq(9, nrow(df_all), by = 14)), ]
RE_val <- df_all[c(seq(10, nrow(df_all), by = 14)), ]

CE_cal <- df_all[c(seq(11, nrow(df_all), by = 14)), ]
CE_val <- df_all[c(seq(12, nrow(df_all), by = 14)), ]

DE_cal <- df_all[c(seq(13, nrow(df_all), by = 14)), ]
DE_val <- df_all[c(seq(14, nrow(df_all), by = 14)), ]


# Mean rank and share of rank 1 is calculated
AVG_rank <- data.frame(rowMeans(apply(-r_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-r_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
r_cal_ranks <- cbind(AVG_rank, shareOne)
names(r_cal_ranks) <- c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-r_val, 2, rank, ties.method =  "min")))
shareOne <- data.frame(apply(apply(-r_val, 2, rank, ties.method =  "min"), 1,
                             count_ones) /  position)
r_val_ranks <- cbind(AVG_rank, shareOne)
names(r_val_ranks) <-  c("Mean Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RMSE_cal, 2, rank,
                                      ties.method =  "min")))
shareOne <- data.frame(apply(apply(RMSE_cal, 2, rank, ties.method =  "min"),
                             1, count_ones) /  position)
RMSE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RMSE_cal_ranks) <-  c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RMSE_val, 2, rank,
                                      ties.method = "min")))
shareOne <- data.frame(apply(apply(RMSE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RMSE_val_ranks <- cbind(AVG_rank, shareOne)
names(RMSE_val_ranks) <-  c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RRSE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(RRSE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RRSE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RRSE_cal_ranks) <-  c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RRSE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(RRSE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RRSE_val_ranks <- cbind(AVG_rank, shareOne)
names(RRSE_val_ranks) <-  c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-d_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-d_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
d_cal_ranks <- cbind(AVG_rank, shareOne)
names(d_cal_ranks) <-  c("Mean Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-d_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-d_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
d_val_ranks <- cbind(AVG_rank, shareOne)
names(d_val_ranks) <-  c("Mean Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-RE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-RE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RE_cal_ranks) <-  c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-RE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-RE_val, 2, rank, ties.method = "min"),
                             1, count_ones) /  position)
RE_val_ranks <- cbind(AVG_rank, shareOne)
names(RE_val_ranks) <-  c("Mean Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-CE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-CE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
CE_cal_ranks <- cbind(AVG_rank, shareOne)
names(CE_cal_ranks) <- c("Mean Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-CE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-CE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
CE_val_ranks <- cbind(AVG_rank, shareOne)
names(CE_val_ranks) <-  c("Mean Rank",  "%rank_1")


AVG_rank <- data.frame(rowMeans(apply(-DE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-DE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
DE_cal_ranks <- cbind(AVG_rank, shareOne)
names(DE_cal_ranks) <- c("Mean Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-DE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-DE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
DE_val_ranks <- cbind(AVG_rank, shareOne)
names(DE_val_ranks) <-  c("Mean Rank",  "%rank_1")




# Results are rbinded together
ranks_together <- rbind(r_cal_ranks, r_val_ranks,
                       RMSE_cal_ranks, RMSE_val_ranks,
                       RRSE_cal_ranks, RRSE_val_ranks,
                       d_cal_ranks, d_val_ranks,
                       RE_cal_ranks, RE_val_ranks,
                       CE_cal_ranks, CE_val_ranks,
                       DE_cal_ranks, DE_val_ranks)

# Those variables have to be defined, solution suggest on Stackoverflow.com
BRNN <- NULL
BRNN_AR <- NULL
BRNN_M <- NULL
BRNN_S1 <- NULL
BRNN_SD <- NULL
MLR <- NULL
MLR_AR <- NULL
MLR_M <- NULL
MLR_S1 <- NULL
MLR_SD <- NULL
MT <- NULL
MT_AR <- NULL
MT_S1 <- NULL
MT_SD <- NULL
Metric <- NULL
Period <- NULL
RF <- NULL
RF_AR <- NULL
RF_M <- NULL
RF_S1 <- NULL
RF_SD <- NULL
bias <- NULL
Method <- NULL
value <- NULL
pred <- NULL
method <- NULL
Year <- NULL
reconstruction <- NULL
edge_observation <- NULL
edge_prediction <- NULL
type <- NULL
data_edge_central <- NULL
DE <- NULL
data_holdout_calibration <- NULL
residuals <- NULL
predicted <- NULL

ranks_together$Method <- methods
ranks_together$Period <- c(rep("cal", length(methods)), rep("val", length(methods)))
ranks_together$Metric <- c(rep("r", length(methods) * 2),
                            rep("RMSE", length(methods) * 2),
                            rep("RRSE", length(methods) * 2),
                            rep("d", length(methods) *2),
                            rep("RE", length(methods) * 2),
                            rep("CE", length(methods) * 2),
                            rep("DE", length(methods) * 2))

colnames(ranks_together)[1] <- "Avg_rank"
togeter_AVG_rank <- reshape2::dcast(ranks_together,
                                  formula = Metric + Period ~ Method,
                                  value.var = c("Avg_rank"))
togeter_AVG_rank$Metric  <- factor(togeter_AVG_rank$Metric,
                                    levels = c("r", "RMSE", "RRSE", "d",
                                               "RE", "CE", "DE"))
togeter_AVG_rank <- togeter_AVG_rank[order(togeter_AVG_rank$Metric), ]
togeter_AVG_rank <- dplyr::select(togeter_AVG_rank, Metric, Period, methods)

colnames(ranks_together)[2] <- "Share_rank1"
together_share1 <- reshape2::dcast(ranks_together,
                                 formula = Metric + Period ~ Method,
                                 value.var = c("Share_rank1"))

together_share1$Metric  <- factor(together_share1$Metric,
                                   levels = c("r", "RMSE", "RRSE", "d",
                                              "RE", "CE", "DE"))

together_share1 <- together_share1[order(together_share1$Metric), ]
together_share1 <- dplyr::select(together_share1, Metric, Period, methods)

###############################################################################

temp_string_mean <- paste("Mean", methods)
temp_string_std <- paste("Std", methods)

emptly_list = list()

sp_position = 1
for (i in 1:length(temp_string_mean)){
  emptly_list[sp_position] <- temp_string_mean[i]
  sp_position <- sp_position + 1

  emptly_list[sp_position] <- temp_string_std[i]
  sp_position <- sp_position + 1
}

df_all_avg_colnames <- as.vector(do.call(cbind, emptly_list))
colnames(df_all_avg) <- df_all_avg_colnames

df_all_avg$Period <- c("cal", "val")
df_all_avg$Metric <- c("r", "r", "RMSE", "RMSE", "RRSE", "RRSE",
                            "d", "d", "RE", "RE", "CE", "CE", "DE", "DE")
row.names(df_all_avg) <- NULL

Results_mean_std <- dplyr::select(df_all_avg, Metric, Period, df_all_avg_colnames)

# Here, we organize the rank and share1 data frame



#1 create string for column names
temp_string_rank <- paste("Avg rank", methods)
colnames(togeter_AVG_rank) <- c("Metric", "Period", temp_string_rank)

temp_string_share1 <- paste("%rank_1", methods)
together_share1 <- together_share1[, -c(1,2)]
colnames(together_share1) <- temp_string_share1

ranks <- cbind(togeter_AVG_rank, together_share1)

# We have to create a vector of strings, cobmination avg rank and %rank 1
emptly_list = list()

sp_position = 1
for (i in 1:length(temp_string_rank)){
  emptly_list[sp_position] <- temp_string_rank[i]
  sp_position <- sp_position + 1

  emptly_list[sp_position] <- temp_string_share1[i]
  sp_position <- sp_position + 1
}

Results_ranks <- dplyr::select(ranks, Metric, Period,
                               as.vector(do.call(cbind, emptly_list)))

##################################################################
# Here is a function to calculate bias
df_MLR_bias$Period <- c("Calibration", "Validation")
df_MLR_bias$Method <- "MLR"

df_BRNN_bias$Period <- c("Calibration", "Validation")
df_BRNN_bias$Method <- "BRNN"

df_MT_bias$Period <- c("Calibration", "Validation")
df_MT_bias$Method <- "MT"

df_RF_bias$Period <- c("Calibration", "Validation")
df_RF_bias$Method <- "RF"

# Subset of methods
start_position = 0
method_list <- list()
for (i in methods){
  start_position <- start_position + 1
  method_list[[start_position]] <- paste0("df_",i,"_bias")
}

method_vector <- unlist(method_list, use.names = FALSE)

empty_LIST <- list()

for (i in 1:length(method_vector)){
  temp_DF <- get(method_vector[i])
  empty_LIST[[i]] <- temp_DF
}

names(empty_LIST) <- method_vector

bias_together <- do.call(rbind, empty_LIST)

bias_together <- melt(bias_together, id.vars = c("Period", "Method"))

bias_together_calibration <- dplyr::filter(bias_together, Period == "Calibration")
bias_together_validation <- dplyr::filter(bias_together, Period == "Validation")

bias_together_calibration$value <- round(bias_together_calibration$value, round_bias_cal)

gg_object_cal <- ggplot(bias_together_calibration, aes(value)) +
  geom_histogram(aes(group = Method), bins = n_bins) +
  # geom_density(aes(group = Method)) +
  geom_vline(xintercept = 0) +
  facet_grid(Method ~ ., scales = "free") +
  theme_bw() +
  xlab("bias") +
  theme(legend.position = "NONE", legend.title = element_blank(),
        text = element_text(size = 15))

bias_together_validation$value <- round(bias_together_validation$value, round_bias_val)

gg_object_val <- ggplot(bias_together_validation, aes(value)) +
  geom_histogram(aes(group = Method), bins = n_bins) +
  # geom_density(aes(group = Method)) +
  geom_vline(xintercept = 0) +
  facet_grid(Method ~ .) +
  theme_bw() +
  xlab("bias") +
  theme(legend.position = "NONE", legend.title = element_blank(),
        text = element_text(size = 15))


########################################################################################
# Here, we create transfer functions for each method
if (numIND == 1) {
  Ind_name <- colnames(dataset)[-DepIndex]
  Dep_name <- colnames(dataset)[DepIndex]

  range_max <- max(dataset[,-DepIndex])
  range_min <- min(dataset[,-DepIndex])

  # I want to increase the range for transfer function, to see greater
  # difference among transfer functions
  diff2 <- (range_max - range_min)/4

  full_range <- data.frame(c1 = NA, c2 = seq(range_min - abs(diff2), range_max + abs(diff2), diff2/1000))
  colnames(full_range) <- c(Dep_name, Ind_name)
  #full_range <- select(full_range, colnames(dataset))

#MLR MODEL
MLR <- lm(formula, data = dataset)
predicted_MLR <- data.frame(pred = predict(MLR, full_range), method = "MLR")

#BRNN Model
capture.output(BRNN <- brnn(formula, data = dataset, BRNN_neurons = BRNN_neurons, verbose = FALSE))
predicted_BRNN <- data.frame(pred = predict(BRNN, full_range), method = "BRNN")

# Model Trees
dataset_ind <- data.frame(dataset[, -DepIndex])
colnames(dataset_ind) <- indep_names
full_range_ind <- data.frame(full_range[, -DepIndex])
colnames(full_range_ind) <- indep_names

dataset_dep <- dataset[, DepIndex]
full_range_dep <- full_range[, DepIndex]

MT_model <- cubist(x = dataset_ind, y = dataset_dep, committees = MT_committees,
                   neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                   unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
predicted_MT <- data.frame(pred = predict(MT_model, full_range_ind), method = "MT")

# Random Forest
RF_model <- randomForest(formula = formula, data = dataset, ntree = RF_ntree,
                         maxnodes = RF_maxnodes, mtry  = RF_mtry,
                         nodesize = RF_nodesize)
predicted_RF <- data.frame(pred = predict(RF_model, full_range), method = "RF")

# Subset of methods
start_position = 0
method_list <- list()
for (i in methods){
  start_position <- start_position + 1
  method_list[[start_position]] <- paste0("predicted_",i)
}

method_vector <- unlist(method_list, use.names = FALSE)

empty_LIST <- list()

for (i in 1:length(method_vector)){
  temp_DF <- get(method_vector[i])
  empty_LIST[[i]] <- temp_DF
}

names(empty_LIST) <- method_vector
predictions <- do.call(rbind, empty_LIST)
predictions$range <- full_range[, 2]
colnames(predictions)

journal_theme <- theme_bw() +
  theme(axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18), text = element_text(size = 18),
        plot.title = element_text(size = 16,  face = "bold"),
        legend.position="bottom", legend.title=element_blank())

dataset1 <-dataset
colnames(dataset1)[DepIndex] <- "pred"
colnames(dataset1)[-DepIndex] <- "range"

plot_1 <- ggplot(predictions, aes(x = range, y = pred)) + geom_line() +
  geom_point(data = dataset1, aes(x = range, y = pred)) +
  facet_grid(method~.) +
  xlab("Range of Independent Variable") +
  ylab("Dependent Variable") +
  ggtitle("Transfer Functions") +
  journal_theme

dataset1$method = NA

plot_2 <- ggplot(predictions, aes(x = range, y = pred, group = method, colour = method)) +
  geom_line(aes(x = range, y = pred, group = method, colour = method, linetype = method), size = 1.05) +
  geom_point(data = dataset1, aes(x = range, y = pred)) +
  xlab("Range of Independent Variable") +
  ylab("Dependent Variable") +
  ggtitle("Transfer Functions Plotted Together") +
  journal_theme

} else if (numIND == 2){

  modelList = list()
  plotList = list()
  plotListNames = methods

  if ("BRNN" %in% methods){
    modelList$BRNN = brnn(formula, dataset, neurons = BRNN_neurons)
  }

  if ("MLR" %in% methods){
    modelList$MLR = lm(formula, dataset)
  }

  if ("MT" %in% methods){
    dataset_ind <- data.frame(dataset[, -DepIndex])
    colnames(dataset_ind) <- indep_names
    dataset_dep <- dataset[, DepIndex]

    modelList$MT = cubist(x = dataset_ind, y = dataset_dep, committees = MT_committees,
                          neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                          unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
  }

  if ("RF" %in% methods){

    modelList$RF = randomForest(formula = formula, data = dataset, ntree = RF_ntree,
                                maxnodes = RF_maxnodes, mtry  = RF_mtry,
                                nodesize = RF_nodesize)
  }

  for (i in 1:length(modelList)){

    model = modelList[[i]]

    x1 = dataset[[indep_name_1]]
    x2 = dataset[[indep_name_2]]
    y = dataset[[DepName]]

    grd <- data.frame(X1 = seq(range(x1)[1],range(x1)[2],len=40),
                      X2 = seq(range(x2)[1],range(x2)[2],len=40))

    names(grd) <- c(indep_names)
    grd$pred <- predict(model, newdata=grd)
    grd <- grd[order(grd[[indep_name_1]],grd[[indep_name_2]]),]


    x1 <- unique(grd[[indep_name_1]])
    x2 <- unique(grd[[indep_name_2]])   # shouldn't have used y
    z=matrix(grd$pred,length(x1),length(x2))

    color_workaround <- rep(0,length(grd$pred))
    dim(color_workaround) <- dim(grd$pred)

    thisName = plotListNames[i]
    plotList[[thisName]] = plot_ly(x = x1, y = x2, z = ~z, name =thisName,
                                   scene = paste0('scene',i)) %>%
      add_surface(opacity = 0.5,  showscale=FALSE, hoverinfo = "name", surfacecolor=color_workaround,
                  hoverlabel = list(color = "black")) %>%
      layout(scene = list(xaxis = list(title = indep_name_1),
                          yaxis = list(title = indep_name_2),
                          zaxis = list(title = DepName))) %>%
      add_trace(x = dataset[[indep_name_1]], y = dataset[[indep_name_2]], z = dataset[[DepName]],   mode = "markers", type = "scatter3d",
                marker = list(size = 5, color = "red", symbol = 104), hoverinfo = "name",
                hoverlabel = list(color = "black")) %>%
      #add_annotations(x = 0.5, y = 0.8, text = paste(plotListNames[i]), showarrow = FALSE) %>%
      layout(showlegend = FALSE) %>% plotly_build()

  }

  if (length (plotList) == 1){
    plot_1 <- subplot(plotList)%>%
      layout(scene = list(aspectmode='cube',
                          yaxis = list(title = indep_name_2),
                          xaxis = list(title = indep_name_1),
                          zaxis = list(title = DepName)))
  }

  if (length (plotList) == 2){
    plot_1 <- subplot(plotList)%>%
      layout(scene = list(domain=list(x=c(0,0.5),y=c(0,1)),
                          aspectmode='cube',
                          yaxis = list(title = indep_name_2),
                          xaxis = list(title = indep_name_1),
                          zaxis = list(title = DepName)),
             scene2 = list(domain=list(x=c(0.5,1),y=c(0,1)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)))
  }

  if (length (plotList) == 3){
    plot_1 <- subplot(plotList)%>%
      layout(scene = list(domain=list(x=c(0,0.5),y=c(0.5,1)),
                          aspectmode='cube',
                          yaxis = list(title = indep_name_2),
                          xaxis = list(title = indep_name_1),
                          zaxis = list(title = DepName)),
             scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)),
             scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)))
  }

  if (length (plotList) == 4){
    plot_1 <- subplot(plotList)%>%
      layout(scene = list(domain=list(x=c(0,0.5),y=c(0.5,1)),
                          aspectmode='cube',
                          yaxis = list(title = indep_name_2),
                          xaxis = list(title = indep_name_1),
                          zaxis = list(title = DepName)),
             scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)),
             scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)),
             scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)))
  }

  if (length (plotList) == 5){
    plot_1 <- subplot(plotList)%>%
      layout(scene = list(domain=list(x=c(0,0.33),y=c(0.5,1)),
                          aspectmode='cube',
                          yaxis = list(title = indep_name_2),
                          xaxis = list(title = indep_name_1),
                          zaxis = list(title = DepName)),
             scene2 = list(domain=list(x=c(0.33,0.66),y=c(0.5,1)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)),
             scene3 = list(domain=list(x=c(0.66,1),y=c(0.5,1)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)),
             scene4 = list(domain=list(x=c(0,0.33),y=c(0,0.5)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName)),
             scene5 = list(domain=list(x=c(0.33,0.66),y=c(0,0.5)),
                           aspectmode='cube',
                           yaxis = list(title = indep_name_2),
                           xaxis = list(title = indep_name_1),
                           zaxis = list(title = DepName))
             )
  }



  for (i in 1:length(modelList)){

    model = modelList[[i]]

    x1 = dataset[[indep_name_1]]
    x2 = dataset[[indep_name_2]]
    y = dataset[[DepName]]

    grd <- data.frame(X1 = seq(range(x1)[1],range(x1)[2],len=40),
                      X2 = seq(range(x2)[1],range(x2)[2],len=40))

    names(grd) <- c(indep_names)
    grd$pred <- predict(model, newdata=grd)
    grd <- grd[order(grd[[indep_name_1]],grd[[indep_name_2]]),]


    x1 <- unique(grd[[indep_name_1]])
    x2 <- unique(grd[[indep_name_2]])   # shouldn't have used y
    z <- matrix(grd$pred,length(x1),length(x2))

    color_workaround <- rep(i,length(grd$pred))
    dim(color_workaround) <- dim(grd$pred)

    color_vector = c("red", "blue", "green", "orange", "brown")

    thisName = plotListNames[i]

    z1 <- list()
    z1[[1]] <- z
    names(z1) <- thisName


    plotList[[thisName]] = plot_ly(x = x1, y = x2, z = ~z1[[1]], name = thisName,
                                   scene = 'scene') %>%
      add_surface(opacity = 0.5, hoverinfo = "name", surfacecolor = color_workaround, colors = color_vector[i],
                  colorbar =list(title=methods[i], dtick = "tick0", tickvals = c()),
                  hoverlabel = list(color = "black")) %>%
      layout(scene = list(xaxis = list(title = indep_name_1),
                          yaxis = list(title = indep_name_2),
                          zaxis = list(title = DepName))) %>%
      add_trace(x = dataset[[indep_name_1]], y = dataset[[indep_name_2]], z = dataset[[DepName]],   mode = "markers", type = "scatter3d",
                marker = list(size = 5, color= "red", symbol = 104), hoverinfo = "none",
                hoverlabel = list(color = "black")) %>%
      layout(showlegend = FALSE) %>% plotly_build()

  }


  plot_2 <- subplot(plotList)%>%
    layout(scene = list(aspectmode='cube',
                        yaxis = list(title = indep_name_2),
                        xaxis = list(title = indep_name_1),
                        zaxis = list(title = DepName)))





  } else {

  plot_1 <- "transfer functions are not available for regression problems with 2 or more independent variables!"
  plot_2 <- "Transfer functions are not available for regression problems with 2 or more independent variables!"

  }



##### Here both data frames are subset with round_df function #############

Results_mean_std <- round_df(Results_mean_std, digits = digits)
Results_ranks <- round_df(Results_ranks, digits = digits)

# Here, all optimized parameters are saved in a data frame, which will be saved as
# a fifth elemnt of the final_list
parameters <- data.frame(
  Method = c("BRNN", "MT", "MT", "MT", "MT","MT" , "MT", "RF", "RF", "RF", "RF"),
  Parameter = c("BRNN_neurons", "MT_committees", "MT_neighbors", "MT_rules", "MT_unbiased",
                "MT_extrapolation", "MT_sample", "RF_mtry", "RF_maxnodes", "RF_ntree", "RF_nodesize"),
  Considered_values = c(
    as.character(paste0(BRNN_neurons_vector, collapse=", ")),
    as.character(paste0(MT_committees_vector, collapse=", ")),
    as.character(paste0(MT_neighbors_vector, collapse=", ")),
    as.character(paste0(MT_rules_vector, collapse=", ")),
    as.character(paste0(MT_unbiased_vector, collapse=", ")),
    as.character(paste0(MT_extrapolation_vector, collapse=", ")),
    as.character(paste0(MT_sample_vector, collapse=", ")),
    as.character(paste0(RF_mtry_vector, collapse=", ")),
    as.character(paste0(RF_maxnodes_vector, collapse=", ")),
    as.character(paste0(RF_ntree_vector, collapse=" ")),
    as.character(paste0(RF_nodesize_vector, collapse=" "))
    ),
  Value = c(BRNN_neurons, MT_committees, MT_neighbors, MT_rules,
            ifelse(MT_unbiased == 1, as.character("TRUE"), as.character("FALSE")),
            MT_extrapolation, MT_sample, RF_mtry, RF_maxnodes, RF_ntree, RF_nodesize))

colnames(parameters) = c("Method", "Parameter", "Considered values", "Selected value")

parameters <- parameters[parameters$Method %in% methods,]

if (optimize == FALSE){
  parameters[, 3] <- NULL

}

                 ############################################################
                 # The reconstruction of climate based on all methods used  #
                 ############################################################
if (is.null(dataset_complete) == FALSE){
MLR_model <- lm(formula, data = dataset)
reconstruction_MLR <- data.frame(reconstruction = predict(MLR_model, dataset_complete), method = "MLR")

#BRNN Model
capture.output(BRNN_model <- brnn(formula, data = dataset, BRNN_neurons = BRNN_neurons, verbose = FALSE))
reconstruction_BRNN <- data.frame(reconstruction = predict(BRNN_model, dataset_complete), method = "BRNN")

# Model Trees
dataset_ind <- data.frame(dataset[, -DepIndex])
colnames(dataset_ind) <- indep_names
dataset_complete_ind <- data.frame(dataset_complete)
colnames(dataset_complete_ind) <- indep_names

dataset_dep <- dataset[, DepIndex]

MT_model <- cubist(x = dataset_ind, y = dataset_dep, committees = MT_committees,
                   neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                   unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
reconstruction_MT <- data.frame(reconstruction = predict(MT_model, dataset_complete_ind), method = "MT")


# Random Forest
RF_model <- randomForest(formula = formula, data = dataset, ntree = RF_ntree,
                         maxnodes = RF_maxnodes, mtry  = RF_mtry,
                         nodesize = RF_nodesize)

reconstruction_RF <- data.frame(reconstruction = predict(RF_model, dataset_complete), method = "RF")



# Subset of methods
start_position = 0
method_list <- list()
for (i in methods){
  start_position <- start_position + 1
  method_list[[start_position]] <- paste0("reconstruction_",i)
}

method_vector <- unlist(method_list, use.names = FALSE)

empty_LIST <- list()

for (i in 1:length(method_vector)){
  temp_DF <- get(method_vector[i])
  empty_LIST[[i]] <- temp_DF
}

names(empty_LIST) <- method_vector
reconstructions <- do.call(rbind, empty_LIST)
reconstructions$Year = as.numeric(row.names(dataset_complete))

# Add plots
plot_3 <- ggplot(reconstructions, aes(x = Year, y = reconstruction, group = method, colour = method)) +
  geom_line() +
  ylab("Reconstruction") +
  ggtitle("Reconstruction of dependent variable") +
  journal_theme

plot_4 <- ggplot(reconstructions, aes(x = Year, y = reconstruction, group = method)) +
  geom_line() +
  ylab("Reconstruction") +
  ggtitle("Reconstruction of dependent variable") +
  facet_wrap(~method, ncol = 1) +
  journal_theme +
  theme(legend.position = 'none')
} else {
  plot_3 <- "The dataset full argument is not supplied, therefore reconstrucions are not calculated."
  plot_4 <- "The dataset full argument is not supplied, therefore reconstrucions are not calculated."
}


##########################################################
#   The comparison of methods, only the edge instances     #
############################################################
# Here, I extract the "edge data", I calibrate models using the central part of the data and
# afterwards use it on the edge data to see, how methods perform in modeling
if (numIND == 1 & edge_share > 0) {
edge_factor <- round2(nrow(dataset)*(edge_share/2),0)

dataset_max <- dplyr::arrange(dataset, desc(dataset[, -DepIndex]))[1:edge_factor,]
dataset_min <- dplyr::arrange(dataset, desc(dataset[, -DepIndex]))[(nrow(dataset)-edge_factor+1):nrow(dataset),]
dataset_edges <- rbind(dataset_max, dataset_min)
dataset_central <- dplyr::arrange(dataset, desc(dataset[, -DepIndex]))[(edge_factor+1):(nrow(dataset)-edge_factor),]

   MLR_model <- lm(formula, data = dataset_central)
   edge_prediction <- data.frame(value = predict(MLR_model, dataset_edges), data_edge_central = "edge", type = "predicted", method = "MLR")
   central_prediction <- data.frame(value = predict(MLR_model, dataset_central), data_edge_central = "central", type = "predicted", method = "MLR")
   edge_observed <- data.frame(value = dataset_edges[, DepIndex], data_edge_central = "edge", type = "observed", method = "MLR")
   central_observed <- data.frame(value = dataset_central[, DepIndex], data_edge_central = "central", type = "observed", method = "MLR")

  MLR_central_edge <- rbind(edge_prediction, central_prediction, edge_observed, central_observed)

  residuals_edge_MLR <- data.frame(residuals = edge_prediction$value - edge_observed$value,
                                      predicted = edge_prediction$value, method = "MLR")

  test_1 <- data.frame(dataset_edges[, DepIndex])
  colnames(test_1) <- "DE_TRICK"
  test_1$sequence <- seq(1:nrow(dataset_edges))
  MLR_DE <- lm(formula = DE_TRICK~., data = test_1)
  DE_predicted <- predict(MLR_DE, test_1)

  #BRNN Model
  capture.output(BRNN_model <- brnn(formula, data = dataset_central, BRNN_neurons = BRNN_neurons, verbose = FALSE))
  edge_prediction <- data.frame(value = predict(BRNN_model, dataset_edges), data_edge_central = "edge", type = "predicted", method = "BRNN")
  central_prediction <- data.frame(value = predict(BRNN_model, dataset_central), data_edge_central = "central", type = "predicted", method = "BRNN")
  edge_observed <- data.frame(value = dataset_edges[, DepIndex], data_edge_central = "edge", type = "observed", method = "BRNN")
  central_observed <- data.frame(value = dataset_central[, DepIndex], data_edge_central = "central", type = "observed", method = "BRNN")
  BRNN_central_edge <- rbind(edge_prediction, central_prediction, edge_observed, central_observed)

  residuals_edge_BRNN <- data.frame(residuals = edge_prediction$value - edge_observed$value,
                                   predicted = edge_prediction$value, method = "BRNN")

  # Model Trees
  dataset_central_ind <- data.frame(dataset_central[, -DepIndex])
  colnames(dataset_central_ind) <- indep_names
  dataset_edges_ind <- data.frame(dataset_edges[, -DepIndex])
  colnames(dataset_edges_ind) <- indep_names

  dataset_central_dep <- dataset_central[, DepIndex]
  dataset_edges_dep <- dataset_edges[, DepIndex]

  MT_model <- cubist(x = dataset_central_ind, y = dataset_central_dep, committees = MT_committees,
                     neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                     unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
  edge_prediction <- data.frame(value = predict(MT_model, dataset_edges_ind), data_edge_central = "edge", type = "predicted", method = "MT")
  central_prediction <- data.frame(value = predict(MT_model, dataset_central_ind), data_edge_central = "central", type = "predicted", method = "MT")
  edge_observed <- data.frame(value = dataset_edges[, DepIndex], data_edge_central = "edge", type = "observed", method = "MT")
  central_observed <- data.frame(value = dataset_central[, DepIndex], data_edge_central = "central", type = "observed", method = "MT")
  MT_central_edge <- rbind(edge_prediction, central_prediction, edge_observed, central_observed)

  residuals_edge_MT <- data.frame(residuals = edge_prediction$value - edge_observed$value,
                                   predicted = edge_prediction$value, method = "MT")

  # Random Forest
  RF_model <- randomForest(formula = formula, data = dataset_central, ntree = RF_ntree, maxnodes = RF_maxnodes, mtry  = RF_mtry,
                           nodesize = RF_nodesize)

  edge_prediction <- data.frame(value = predict(RF_model, dataset_edges), data_edge_central = "edge", type = "predicted", method = "RF")
  central_prediction <- data.frame(value = predict(RF_model, dataset_central), data_edge_central = "central", type = "predicted", method = "RF")
  edge_observed <- data.frame(value = dataset_edges[, DepIndex], data_edge_central = "edge", type = "observed", method = "RF")
  central_observed <- data.frame(value = dataset_central[, DepIndex], data_edge_central = "central", type = "observed", method = "RF")
  RF_central_edge <- rbind(edge_prediction, central_prediction, edge_observed, central_observed)

  residuals_edge_RF <- data.frame(residuals = edge_prediction$value - edge_observed$value,
                                   predicted = edge_prediction$value, method = "RF")

  edge_central_data <- rbind(MLR_central_edge,
                             BRNN_central_edge,
                             MT_central_edge,
                             RF_central_edge)

  # Here i define the residual dataframe for the holdout
  residuals_edge = rbind(residuals_edge_MLR,
                            residuals_edge_BRNN,
                            residuals_edge_MT,
                            residuals_edge_RF)

  edge_central_data$DE = 0

  DE_data_MLR <- data.frame(value = DE_predicted, data_edge_central = "edge", type = "predicted", method = "MLR", DE = 1)
  DE_data_BRNN <- data.frame(value = DE_predicted, data_edge_central = "edge", type = "predicted", method = "BRNN", DE = 1)
  DE_data_MT <- data.frame(value = DE_predicted, data_edge_central = "edge", type = "predicted", method = "MT", DE = 1)
  DE_data_RF <- data.frame(value = DE_predicted, data_edge_central = "edge", type = "predicted", method = "RF", DE = 1)

  DE_data <- rbind(DE_data_MLR,DE_data_BRNN,DE_data_MT,DE_data_RF)

  edge_central_data <- rbind(edge_central_data, DE_data)

  edge_central_data <- dplyr::filter(edge_central_data, method %in% methods)

  edge_results <- edge_central_data %>%
                    group_by(method) %>%
                    summarise(r_central = cor(value[data_edge_central == "central" & type == "predicted"  & DE == 0 ],
                                          value[data_edge_central == "central" & type == "observed" & DE == 0 ]),
                             r_edge = cor(value[data_edge_central == "edge" & type == "predicted" & DE == 0 ],
                                      value[data_edge_central == "edge" & type == "observed" & DE == 0 ]),

                             RMSE_central = MLmetrics::RMSE(value[data_edge_central == "central" & type == "predicted" & DE == 0 ],
                                                            value[data_edge_central == "central" & type == "observed" & DE == 0 ]),
                             RMSE_edge =  MLmetrics::RMSE(value[data_edge_central == "edge" & type == "predicted" & DE == 0 ],
                                                          value[data_edge_central == "edge" & type == "observed" & DE == 0 ]),

                             RRSE_central = MLmetrics::RRSE(value[data_edge_central == "central" & type == "predicted" & DE == 0 ],
                                                            value[data_edge_central == "central" & type == "observed" & DE == 0 ]),
                             RRSE_edge = MLmetrics::RRSE(value[data_edge_central == "edge" & type == "predicted" & DE == 0 ],
                                                         value[data_edge_central == "edge" & type == "observed" & DE == 0 ]),
                             d_central  = 1 - (sum((value[data_edge_central == "central" & type == "observed" & DE == 0 ] -
                                                      value[data_edge_central == "central" & type == "predicted" & DE == 0 ]) ^ 2)) /
                               sum((abs(value[data_edge_central == "central" & type == "predicted" & DE == 0 ] -
                                          mean(value[data_edge_central == "central" & type == "observed" & DE == 0 ])) +
                                      abs(value[data_edge_central == "central" & type == "observed" & DE == 0 ] -
                                            mean(value[data_edge_central == "central" & type == "observed" & DE == 0 ]))) ^ 2),
                             d_edge = 1 - (sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                  value[data_edge_central == "edge" & type == "predicted" & DE == 0 ]) ^ 2)) /
                               sum((abs(value[data_edge_central == "edge" & type == "predicted" & DE == 0 ] -
                                          mean(value[data_edge_central == "edge" & type == "observed" & DE == 0 ])) +
                                      abs(value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                            mean(value[data_edge_central == "edge" & type == "observed" & DE == 0 ]))) ^ 2),

                             RE_edge = 1 - (sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                   value[data_edge_central == "edge" & type == "predicted" & DE == 0 ]) ^ 2) /
                                              sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                     mean(value[data_edge_central == "central" & type == "observed" & DE == 0 ])) ^ 2)),

                             CE_edge  = 1 - (sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                    value[data_edge_central == "edge" & type == "predicted" & DE == 0 ]) ^ 2) /
                                               sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                      mean(value[data_edge_central == "edge" & type == "observed" & DE == 0 ])) ^ 2)),
                             DE_edge = 1 - (sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                   value[data_edge_central == "edge" & type == "predicted" & DE == 0 ]) ^ 2) /
                                              sum((value[data_edge_central == "edge" & type == "observed" & DE == 0 ] -
                                                     value[data_edge_central == "edge" & type == "predicted" & DE == 1 ]) ^ 2))


                               )


# summary(edge_central_data)
# test_observed: edge_central_data$value[edge_central_data$data_edge_central == "edge" & edge_central_data$type == "predicted" & edge_central_data$DE == 1 & edge_central_data$method == "MLR"]
# test_predicted: value[data_edge_central == "edge" & type == "predicted"]
# train_observed: value[data_edge_central == "central" & type == "observed"]
# train_predicted: value[data_edge_central == "central" & type == "predicted"]



  edge_results_t <- t(edge_results[,-1])
  colnames(edge_results_t) <- edge_results$method
  edge_results_t <- round(edge_results_t, digits)
  edge_results_t <- data.frame(edge_results_t, Period =  c("cal_central", "val_edge",
                                                           "cal_central", "val_edge",
                                                           "cal_central", "val_edge",
                                                           "cal_central", "val_edge",
                                                           "val_edge","val_edge",
                                                           "val_edge"
                                                            ),
                               Metric = c("r", "r", "RMSE", "RMSE", "RRSE", "RRSE",
                             "d", "d", "RE", "CE", "DE"))
  edge_results_t <- select(edge_results_t, Metric, Period, methods)
  row.names(edge_results_t) <- NULL


} else {
  edge_results_t <- "No edge data is available for regression problems with more than 1 independent variable."
}

# Here, we once again subset the data frames with preformance metrics, to exclude
# RE, CE and DE for the calibration data
Results_mean_std <- Results_mean_std[-c(9, 11, 13), ]
Results_ranks <- Results_ranks[-c(9, 11, 13), ]




##########################################################
#   holdout calculations    #
############################################################

if (!is.null(holdout)) {

  MLR_model <- lm(formula, data = dataset)
  holdout_prediction <- data.frame(value = predict(MLR_model, dataset_holdout), data_holdout_calibration = "holdout", type = "predicted", method = "MLR")
  calibration_prediction <- data.frame(value = predict(MLR_model, dataset), data_holdout_calibration = "calibration", type = "predicted", method = "MLR")
  holdout_observed <- data.frame(value = dataset_holdout[, DepIndex], data_holdout_calibration = "holdout", type = "observed", method = "MLR")
  calibration_observed <- data.frame(value = dataset[, DepIndex], data_holdout_calibration = "calibration", type = "observed", method = "MLR")

  residuals_holdout_MLR <- data.frame(residuals = holdout_prediction$value - holdout_observed$value,
                                      predicted = holdout_prediction$value, method = "MLR")

  MLR_calibration_holdout <- rbind(holdout_prediction, calibration_prediction, holdout_observed, calibration_observed)

  test_1 <- data.frame(dataset_holdout[, DepIndex])
  colnames(test_1) <- "DE_TRICK"
  test_1$sequence <- seq(1:nrow(dataset_holdout))
  MLR_DE <- lm(formula = DE_TRICK~., data = test_1)
  DE_predicted <- predict(MLR_DE, test_1)

  #BRNN Model
  capture.output(BRNN_model <- brnn(formula, data = dataset, BRNN_neurons = BRNN_neurons, verbose = FALSE))
  holdout_prediction <- data.frame(value = predict(BRNN_model, dataset_holdout), data_holdout_calibration = "holdout", type = "predicted", method = "BRNN")
  calibration_prediction <- data.frame(value = predict(BRNN_model, dataset), data_holdout_calibration = "calibration", type = "predicted", method = "BRNN")
  holdout_observed <- data.frame(value = dataset_holdout[, DepIndex], data_holdout_calibration = "holdout", type = "observed", method = "BRNN")
  calibration_observed <- data.frame(value = dataset[, DepIndex], data_holdout_calibration = "calibration", type = "observed", method = "BRNN")
  BRNN_calibration_holdout <- rbind(holdout_prediction, calibration_prediction, holdout_observed, calibration_observed)

  residuals_holdout_BRNN <- data.frame(residuals = holdout_prediction$value - holdout_observed$value,
                                       predicted = holdout_prediction$value, method = "BRNN")

  # Model Trees
  dataset_ind <- data.frame(dataset[, -DepIndex])
  colnames(dataset_ind) <- indep_names
  dataset_holdout_ind <- data.frame(dataset_holdout[, -DepIndex])
  colnames(dataset_holdout_ind) <- indep_names

  dataset_dep <- dataset[, DepIndex]
  dataset_holdout_dep <- dataset_holdout[, DepIndex]

  MT_model <- cubist(x = dataset_ind, y = dataset_dep, committees = MT_committees,
                     neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                     unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
  holdout_prediction <- data.frame(value = predict(MT_model, dataset_holdout_ind), data_holdout_calibration = "holdout", type = "predicted", method = "MT")
  calibration_prediction <- data.frame(value = predict(MT_model, dataset_ind), data_holdout_calibration = "calibration", type = "predicted", method = "MT")
  holdout_observed <- data.frame(value = dataset_holdout[, DepIndex], data_holdout_calibration = "holdout", type = "observed", method = "MT")
  calibration_observed <- data.frame(value = dataset[, DepIndex], data_holdout_calibration = "calibration", type = "observed", method = "MT")
  MT_calibration_holdout <- rbind(holdout_prediction, calibration_prediction, holdout_observed, calibration_observed)

  residuals_holdout_MT <- data.frame(residuals = holdout_prediction$value - holdout_observed$value,
                                     predicted = holdout_prediction$value, method = "MT")

  # Random Forest
  RF_model <- randomForest(formula = formula, data = dataset, mtry = RF_mtry,
                           maxnodes = RF_maxnodes, ntree = RF_ntree, nodesize = RF_nodesize)

  holdout_prediction <- data.frame(value = predict(RF_model, dataset_holdout), data_holdout_calibration = "holdout", type = "predicted", method = "RF")
  calibration_prediction <- data.frame(value = predict(RF_model, dataset), data_holdout_calibration = "calibration", type = "predicted", method = "RF")
  holdout_observed <- data.frame(value = dataset_holdout[, DepIndex], data_holdout_calibration = "holdout", type = "observed", method = "RF")
  calibration_observed <- data.frame(value = dataset[, DepIndex], data_holdout_calibration = "calibration", type = "observed", method = "RF")
  RF_calibration_holdout <- rbind(holdout_prediction, calibration_prediction, holdout_observed, calibration_observed)

  holdout_calibration_data <- rbind(MLR_calibration_holdout,
                             BRNN_calibration_holdout,
                             MT_calibration_holdout,
                             RF_calibration_holdout)

  residuals_holdout_RF <- data.frame(residuals = holdout_prediction$value - holdout_observed$value,
                                     predicted = holdout_prediction$value, method = "RF")

  # Here i define the residual dataframe for the holdout
  residuals_holdout = rbind(residuals_holdout_MLR,
                            residuals_holdout_BRNN,
                            residuals_holdout_MT,
                            residuals_holdout_RF)

  residuals_holdout <- filter(residuals_holdout, method %in% methods)


  holdout_calibration_data$DE = 0

  DE_data_MLR <- data.frame(value = DE_predicted, data_holdout_calibration = "holdout", type = "predicted", method = "MLR", DE = 1)
  DE_data_BRNN <- data.frame(value = DE_predicted, data_holdout_calibration = "holdout", type = "predicted", method = "BRNN", DE = 1)
  DE_data_MT <- data.frame(value = DE_predicted, data_holdout_calibration = "holdout", type = "predicted", method = "MT", DE = 1)
  DE_data_RF <- data.frame(value = DE_predicted, data_holdout_calibration = "holdout", type = "predicted", method = "RF", DE = 1)

  DE_data <- rbind(DE_data_MLR,DE_data_BRNN,DE_data_MT,DE_data_RF)

  holdout_calibration_data <- rbind(holdout_calibration_data, DE_data)

  holdout_calibration_data <- dplyr::filter(holdout_calibration_data, method %in% methods)

  holdout_results <- holdout_calibration_data %>%
    group_by(method) %>%
    summarise(r_calibration = cor(value[data_holdout_calibration == "calibration" & type == "predicted"  & DE == 0 ],
                              value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ]),
              r_holdout = cor(value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ],
                           value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ]),

              RMSE_calibration = MLmetrics::RMSE(value[data_holdout_calibration == "calibration" & type == "predicted" & DE == 0 ],
                                             value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ]),
              RMSE_holdout =  MLmetrics::RMSE(value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ],
                                           value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ]),

              RRSE_calibration = MLmetrics::RRSE(value[data_holdout_calibration == "calibration" & type == "predicted" & DE == 0 ],
                                             value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ]),
              RRSE_holdout = MLmetrics::RRSE(value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ],
                                          value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ]),
              d_calibration  = 1 - (sum((value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ] -
                                       value[data_holdout_calibration == "calibration" & type == "predicted" & DE == 0 ]) ^ 2)) /
                sum((abs(value[data_holdout_calibration == "calibration" & type == "predicted" & DE == 0 ] -
                           mean(value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ])) +
                       abs(value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ] -
                             mean(value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ]))) ^ 2),
              d_holdout = 1 - (sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                   value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ]) ^ 2)) /
                sum((abs(value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ] -
                           mean(value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ])) +
                       abs(value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                             mean(value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ]))) ^ 2),

              RE_holdout = 1 - (sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                    value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ]) ^ 2) /
                               sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                      mean(value[data_holdout_calibration == "calibration" & type == "observed" & DE == 0 ])) ^ 2)),

              CE_holdout  = 1 - (sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                     value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ]) ^ 2) /
                                sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                       mean(value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ])) ^ 2)),
              DE_holdout = 1 - (sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                    value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 0 ]) ^ 2) /
                               sum((value[data_holdout_calibration == "holdout" & type == "observed" & DE == 0 ] -
                                      value[data_holdout_calibration == "holdout" & type == "predicted" & DE == 1 ]) ^ 2))


    )


  holdout_results_t <- t(holdout_results[,-1])
  colnames(holdout_results_t) <- holdout_results$method
  holdout_results_t <- round(holdout_results_t, digits)
  holdout_results_t <- data.frame(holdout_results_t, Period =  c("calibration", "holdout",
                                                           "calibration", "holdout",
                                                           "calibration", "holdout",
                                                           "calibration", "holdout",
                                                           "holdout","holdout",
                                                           "holdout"
  ),
  Metric = c("r", "r", "RMSE", "RMSE", "RRSE", "RRSE",
             "d", "d", "RE", "CE", "DE"))
  holdout_results_t <- select(holdout_results_t, Metric, Period, methods)
  row.names(holdout_results_t) <- NULL


} else {
  holdout_results_t <- "No holdout data was defined."
}

###################### Residual diagnostic ###################################################################
residuals_df <- data.frame(residuals = c(), predicted = c(), method = c())

data_observed <- dataset[,DepIndex]

#MLR model
if ("MLR" %in% methods){
MLR_model <- lm(formula, data = dataset)
MLR_predicted <- predict(MLR_model, dataset)
residuals_MLR <- MLR_predicted - data_observed
LOESS_MLR <- loess.smooth(MLR_predicted, residuals_MLR)
residuals_df_MLR <- data.frame(residuals = residuals_MLR,
                               predicted = MLR_predicted,
                               method = "MLR")
residuals_df <- rbind(residuals_df, residuals_df_MLR)

}

# BRNN model
if ("BRNN" %in% methods){
capture.output(BRNN_model <- brnn(formula, data = dataset, BRNN_neurons = BRNN_neurons, verbose = FALSE))
BRNN_predicted <- predict(BRNN_model, dataset)
residuals_BRNN <- BRNN_predicted - data_observed
LOESS_BRNN <- loess.smooth(BRNN_predicted, residuals_BRNN)
residuals_df_BRNN <- data.frame(residuals = residuals_BRNN,
                                predicted = BRNN_predicted,
                                method = "BRNN")
residuals_df <- rbind(residuals_df, residuals_df_BRNN)

}

# MT model
if ("MF" %in% methods){
dataset_ind <- data.frame(dataset[, -DepIndex])
colnames(dataset_ind) <- indep_names
dataset_dep <- dataset[, DepIndex]
MT_model <- cubist(x = dataset_ind, y = dataset_dep, committees = MT_committees,
                   neighbors = MT_neighbors, cubistControl(rules = MT_rules,
                   unbiased = MT_unbiased, extrapolation = MT_extrapolation, sample = MT_sample))
MT_predicted <- predict(MT_model, dataset_ind)
residuals_MT <- MT_predicted - data_observed
LOESS_MT <- loess.smooth(MT_predicted, residuals_MT)
residuals_df_MT <- data.frame(residuals = residuals_MT,
                              predicted = MT_predicted,
                              method = "MT")
residuals_df <- rbind(residuals_df, residuals_df_MT)

}

# RF model
if ("RF" %in% methods){
  RF_model <- randomForest(formula = formula, data = dataset, ntree = RF_ntree,
                           maxnodes = RF_maxnodes, mtry  = RF_mtry,
                           nodesize = RF_nodesize)
  RF_predicted <- predict(RF_model, dataset)
  residuals_RF <- RF_predicted - data_observed
  LOESS_RF <- loess.smooth(RF_predicted, residuals_RF)
  residuals_df_RF <- data.frame(residuals = residuals_RF,
                                predicted = RF_predicted,
                                method = "RF")
  residuals_df <- rbind(residuals_df, residuals_df_RF)
}

# Here we create plots for calibration

Normal_QQ_cal1 <- ggplot(residuals_df, aes(sample = residuals)) +
                    facet_grid(method~.)+ stat_qq() + stat_qq_line() +
                    xlab("Theoretical quantiles") +
                    ylab("Residuals") +
                    ggtitle("Normal Q-Q plot, Calibration data") +
                    theme_bw() +
                    theme(text = element_text(size = 15))

Residuals_vs_fitted_cal1 <- ggplot(residuals_df, aes(y = residuals, x = predicted)) +
  geom_point() +
  stat_smooth(method="loess", se = FALSE)+geom_hline(yintercept=0, col="red", linetype="dashed") +
  facet_grid(method ~ .) +
  xlab("Fitted values") +
  ylab("Residuals") +
  ggtitle("Residual vs Fitted Plot, Calibration data") +
  theme_bw() +
  theme(text = element_text(size = 15), axis.title.y = element_blank())

# Residual plots for the holdout data
if (!is.null(holdout)) {
  Normal_QQ_holdout1 <- ggplot(residuals_holdout, aes(sample = residuals)) +
  facet_grid(method~.)+ stat_qq() + stat_qq_line() +
  xlab("Theoretical quantiles") +
  ylab("Residuals") +
  ggtitle("Normal Q-Q plot, Holdout data") +
  theme_bw() +
  theme(text = element_text(size = 15))

  Residuals_vs_fitted_holdout1 <- ggplot(residuals_holdout, aes(y = residuals, x = predicted)) +
  geom_point() +
  stat_smooth(method="loess", se = FALSE)+geom_hline(yintercept=0, col="red", linetype="dashed") +
  facet_grid(method ~ .) +
  xlab("Fitted values") +
  ylab("Residuals") +
  ggtitle("Residual vs Fitted Plot, Holdout data") +
  theme_bw() +
  theme(text = element_text(size = 15), axis.title.y = element_blank())

} else {
  Normal_QQ_holdout1 <- "No holdout data was defined."
  Residuals_vs_fitted_holdout1 <- "No holdout data was defined."
}

# Residual plots for the edge data
if (numIND < 2 & edge_share > 0){
  Normal_QQ_edge1 <- ggplot(residuals_edge, aes(sample = residuals)) +
    facet_grid(method~.)+ stat_qq() + stat_qq_line() +
    xlab("Theoretical quantiles") +
    ylab("Residuals") +
    ggtitle("Normal Q-Q plot, edge data") +
    theme_bw() +
    theme(text = element_text(size = 15))

  Residuals_vs_fitted_edge1 <- ggplot(residuals_edge, aes(y = residuals, x = predicted)) +
    geom_point() +
    stat_smooth(method="loess", se = FALSE)+geom_hline(yintercept=0, col="red", linetype="dashed") +
    facet_grid(method ~ .) +
    xlab("Fitted values") +
    ylab("Residuals") +
    ggtitle("Residual vs Fitted Plot, edge data") +
    theme_bw() +
    theme(text = element_text(size = 15), axis.title.y = element_blank())
} else {
  Residuals_vs_fitted_edge1 <- "No edge data is available for regression problems with 2 or more independent variables."
  Normal_QQ_edge1 <- "No edge data is available for regression problems with 2 or more independent variables."
}


  #####################################################################

# If Calibration and Validation data should be returned, then this is our final results
  final_list <- list(mean_std = Results_mean_std, ranks = Results_ranks,
                     edge_results = edge_results_t,
                     holdout_results = holdout_results_t,
                     bias_cal = suppressMessages(gg_object_cal),
                     bias_val = suppressMessages(gg_object_val),
                     transfer_functions = suppressMessages(plot_1),
                     transfer_functions_together = suppressMessages(plot_2),
                     parameter_values = parameters,
                     PCA_output = PCA_result,
                     reconstructions = suppressMessages(plot_4),
                     reconstructions_together = suppressMessages(plot_3),
                     normal_QQ_cal = suppressMessages(Normal_QQ_cal1),
                     normal_QQ_holdout = suppressMessages(Normal_QQ_holdout1),
                     normal_QQ_edge = suppressMessages(Normal_QQ_edge1),
                     residuals_vs_fitted_cal = suppressMessages(Residuals_vs_fitted_cal1),
                     residuals_vs_fitted_holdout = suppressMessages(Residuals_vs_fitted_holdout1),
                     residuals_vs_fitted_edge = suppressMessages(Residuals_vs_fitted_edge1)
                     )

return(final_list) # Return the final list

}
