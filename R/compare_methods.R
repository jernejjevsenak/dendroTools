#' compare_methods
#'
#' Calculates performance metrics for train and test data of different
#' regression methods: multiple linear regression (MLR), ridge and lasso
#' regression, artificial neural networks with Bayesian regularization
#' training algorithm (ANN), M5P model trees (MT), model trees with bagging
#' (BMT), polynomial regression (POLY) and random forest of regression
#' trees (RF). Calculated performance metrics are correlation coefficient,
#' root mean squared error (RMSE), root relative squared error (RSSE), index
#' of agreement (d), reduction of error (RE), coefficient of efficiency
#' (CE) and mean bias.
#'
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param dataset a data frame with dependent and independent variables as
#' columns and (optional) years as row names.
#' @param k number of folds for cross-validation
#' @param repeats number of cross-validation repeats. Should be equal or more
#' than 2.
#' @param use_caret if set to TRUE, the package caret will be used to tune parameters
#' for regression methods
#' @param ANN_neurons positive integer that indicates the number of neurons used
#'  for brnn method
#' @param MT_M minimum number of instances used by model trees
#' @param MT_N unpruned (argument for model trees)
#' @param MT_U unsmoothed (argument for model trees)
#' @param MT_R use regression trees (argument for model trees)
#' @param BMT_P bagSizePercent (argument for bagging of model trees)
#' @param BMT_I number of iterations (argument for bagging of model trees)
#' @param BMT_M minimum number of instances used by model trees
#' @param BMT_N unpruned (argument for bagging of model trees)
#' @param BMT_U unsmoothed (argument for bagging of model trees)
#' @param BMT_R use regression trees (argument for bagging of model trees)
#' @param RF_mtry Number of variables randomly sampled as candidates at each
#' split (argument for random forest)
#' @param RF_maxnodes Maximum number of terminal nodes trees in the forest can
#' have (argument for random forest)
#' @param RF_ntree Number of trees to grow (argument for random forest)
#' @param RIDGE_lambda lambda argument for ridge regression
#' @param LASSO_lambda lambda argument for lasso regression
#' @param seed_factor an intiger that will be used to change the seed options
#' for different repeats. set.seed(seed_factor*5)
#' @param returns A character vector that specifies, whether a calibration and/ or
#' validation results should be returned.
#' @param digits intiger of number of digits to be displayed in the final
#' result tables
#' @param blocked_CV default is FALSE, if changed to TRUE, blocked cross-validation
#' will be used to compare regression methods.
#' @param PCA_transformation if set to TRUE, all independet variables will be
#' transformed using PCA transformation.
#' @param log_preprocess if set to TRUE, variables will be transformed with
#' logarithmic transformation before used in PCA
#' @param components_selection character string specifying how to select the Principal
#' Components used as predictors.
#' There are three options: "automatic", "manual" and "plot_selection". If
#' parameter is set to automatic, all scores with eigenvalues above 1 will be
#' selected. This threshold could be changed by changing the
#' eigenvalues_threhold argument. If parameter is set to "manual", user should
#' set the number of components with N_components argument. If component
#' selection is se to "plot_selection", Scree plot will be shown and user must
#' manually enter the number of components used as predictors.
#' @param eigenvalues_threhold threshold for automatic selection of Principal Components
#' @param N_components number of Principal Components used as predictors
#' @param polynomial_formula a symbolic description of polinomial model to be fitted
#' @param round_bias_cal number of digits for bias in calibration period. Effects
#' the outlook of the final ggplot  of mean bias for calibration data (element 3 of
#' the output list)
#' @param round_bias_val number of digits for bias in validation period. Effects
#' the outlook of the final ggplot of mean bias for validation data (element 4 of
#' the output list)
#'
#' @return a list with five elements:
#'          $mean_std,  data frame with calculated metrics for five regression methods.
#'           For each regression method and each calculated metric, mean and standard
#'           deviation are given.
#'          $ranks, data frame with ranks of calculated metrics: average rank and %rank_1
#'           are given.
#'          $bias_cal, ggplot object of mean bias for calibration data.
#'          $bias_val, ggplot object of mean bias for validation data. If returns argument
#'           is set to return only "Calibration" or "Validation" results, only the three
#'           relevant elements will be returned in the list.
#'          $parameters, a data frame with specifications of parameters used for different
#'           regression methods.
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
#' data(example_dataset_1)
#'
#' # An example with default settings of machine learning algorithms
#' experiment_1 <- compare_methods(formula = MVA~.,
#' dataset = example_dataset_1, k = 5, repeats = 100,
#' returns = c("Calibration", "Validation"), blocked_CV = TRUE, PCA_transformation = FALSE,
#' components_selection = "plot_selection", use_caret = TRUE,
#' polynomial_formula = "MVA ~ T_APR + T_aug_sep + T_APR^2")
#' experiment_1[[1]] # See a data frame results of mean and standard deviation
#' # for different methods
#' experiment_1[[2]] # See a data frame results of average rank and share of
#' # rank 1 for different methods
#' experiment_1[[3]] # See a ggplot of mean bias for calibration data
#' experiment_1[[4]] # See a ggplot of mean bias for validation data
#' experiment_1[[5]] # Data frame with parameters used for regression methods
#'
#' experiment_2 <- compare_methods(formula = MVA ~ .,
#' dataset = example_dataset_1, k = 5, repeats = 100, ANN_neurons = 1,
#' MT_M = 4, MT_N = FALSE, MT_U = FALSE, MT_R = FALSE, BMT_P = 100,
#' BMT_I = 100, BMT_M = 4, BMT_N = FALSE, BMT_U = FALSE, BMT_R = FALSE,
#' RF_mtry = 0, RF_maxnodes = 4, RF_ntree = 200, seed_factor = 5,
#' returns = c("Calibration"))
#' experiment_2[[1]]
#' experiment_2[[2]]
#' experiment_2[[3]]
#'
#' experiment_3 <- compare_methods(formula = MVA~.,
#' dataset = example_dataset_1, k = 5, repeats = 10,
#' use_caret = TRUE, returns = c("Validation"))
#' experiment_3[[1]]
#' experiment_3[[2]]
#' experiment_3[[3]]

compare_methods <- function(formula, dataset, k = 3, repeats = 2,
                            use_caret = TRUE,
                            ANN_neurons = 1, MT_M = 4, MT_N = F, MT_U = F,
                            MT_R = F, BMT_P = 100, BMT_I = 100, BMT_M = 4,
                            BMT_N = F, BMT_U = F, BMT_R = F, RF_mtry = 0,
                            RF_maxnodes = 4, RF_ntree = 200, RIDGE_lambda = 0.1,
                            LASSO_lambda = 0.1, polynomial_formula = "",
                            seed_factor = 5,
                            returns = c("Calibration", "Validation"),
                            digits = 3, blocked_CV = FALSE,
                            PCA_transformation = FALSE, log_preprocess = TRUE,
                            components_selection = 'automatic',
                            eigenvalues_threhold = 1, N_components = 2,
                            round_bias_cal = 15, round_bias_val = 4) {

dataset <- data.frame(dataset) # dataset needs to be of class data.frame!

# This function is used to calculate metrics r, RMSE, RRSE, d, RE, CE and bias
# for train and test data

#############################################################################

# Here, empty lists are defined, where calculations will be stored. Empty lists
# for bias are defined separately, since bias should not be averaged. It is
# later given as density plots
list_MLR <- list()
list_ANN <- list()
list_MT <- list()
list_BMT <- list()
list_RF <- list()
list_RIDGE <- list()
list_LASSO <- list()
list_POLY <- list()

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
    subset_vector <- PCA_result$sdev > eigenvalues_threhold
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
}

# Here we fit a lm model, just to get information about the number of independet variables;
# when formula is used in the form of: y~., we don't know the number of independet variables
# this information is used later

quazi_mod <- lm(formula, data = dataset)
numIND <- length(quazi_mod[[1]]) - 1
indep_names <- colnames(quazi_mod[[7]][[1]])[-1]

# Warnings
if (polynomial_formula == ""){
  warning("No polynomial formula specified. Polynomial calculations will not be performed!")
}

if (numIND < 2){
  warning("Only one independet variable is used. RIDGE and LASSO regression will not be used!")
}

# Here we use caret package to tune parameters of different methods

if (use_caret == TRUE){

  model = NULL

  # Optimization for ANN
  capture.output(model <- train(formula, data = dataset, method = "brnn"))
  ANN_neurons = as.numeric(model[[6]][1])

  # Optimization for MT
  capture.output(model <- train(formula, data = dataset, method = "M5"))
  if (as.matrix(model[[6]][1][1]) == 'Yes'){
    MT_N = FALSE
  } else {
    MT_N = TRUE
  }

  if (as.matrix(model[[6]][2]) == 'Yes'){
    MT_U = FALSE
  } else {
    MT_U = TRUE
  }

  if (as.matrix(model[[6]][3]) == 'No'){
    MT_R = FALSE
  } else {
    MT_R = TRUE
  }

  # Optimization for BMT, just take MT rezults
  BMT_N <- MT_N
  BMT_U <- MT_U
  BMT_R <- MT_R

  # Optimization for Random Forest
  suppressWarnings(capture.output(model <- train(formula, data = dataset, method = "rf")))
  RF_mtry = as.numeric(model[[6]][1])

  # Optimizacija za ridge in lasso regression
  if (numIND < 2){
    RIDGE_lambda <- NA
    LASSO_lambda <- NA
  } else {
  x = model.matrix(formula, dataset)[,-DepIndex]
  y = as.matrix(dataset[,DepIndex])

  ridge.mod <- glmnet(y = y, x = x,  alpha = 0)
  cv.out = cv.glmnet(x, y, alpha = 0)
  RIDGE_lambda <- cv.out$lambda.min

  lasso.mod <- glmnet(y = y, x = x,  alpha = 1)
  cv.out = cv.glmnet(x, y, alpha = 1)
  LASSO_lambda <- cv.out$lambda.min
  }
}

##################################################################################
##################################################################################
##################################################################################
# Normal cross-validation with repeats.

if (blocked_CV == FALSE){

# create progress bar
pb <- txtProgressBar(min = 0, max = repeats, style = 3)

b = 0 # place holder for saving rezults

for (m in 1:repeats){

foldi <- seq(1:k)
foldi <- paste("fold_", foldi)

#Randomly shuffle the data
set.seed(seed_factor * 5)
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
  MLR <- lm(formula, data = train)
  train_predicted <- predict(MLR, train)
  test_predicted <- predict(MLR, test)
  train_observed <- train[, DepIndex]
  test_observed <- test[, DepIndex]
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                         train_observed, test_observed, digits = 15)
  list_MLR[[b]] <- calculations

  #ANN Model
  capture.output(ANN <- brnn(formula, data = train, ANN_neurons = ANN_neurons, verbose = FALSE))
  train_predicted <- predict(ANN, train)
  test_predicted <- predict(ANN, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15)
  list_ANN[[b]] <- calculations

  # Model Trees
  MT_model <- M5P(formula, data = train,
                  control = Weka_control(M = MT_M, N =  MT_N, U = MT_U,
                                         R = MT_R))
  train_predicted <- predict(MT_model, train)
  test_predicted <- predict(MT_model, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15)
  list_MT[[b]] <- calculations



  #M5 Model with bagging
  BMT_model <- Bagging(formula,
                       data = train,
                       control = Weka_control(P = BMT_P, I = BMT_I,
                                              W = list("weka.classifiers.trees.M5P",
                                                       M = BMT_M, N = BMT_N,
                                                       U = BMT_U, R = BMT_R)))
  train_predicted <- predict(BMT_model, train)
  test_predicted <- predict(BMT_model, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15)
  list_BMT[[b]] <- calculations

  # Random Forest
  RF_model <- randomForest(formula = formula, data = train, RF_mtry = RF_mtry,
                               RF_maxnodes = RF_maxnodes, RF_ntree = RF_ntree)
  train_predicted <- predict(RF_model, train)
  test_predicted <- predict(RF_model, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15)
  list_RF[[b]] <- calculations



  # Ridge Regression!
  x = model.matrix(formula, train)[,-DepIndex]
  y = as.matrix(train[,DepIndex])

  if (numIND < 2){
    ridge.mod <- lm(formula, data = train)
    train_predicted <- predict(ridge.mod, train)
    test_predicted <- predict(ridge.mod, test)
  } else {
    ridge.mod <- glmnet(y = y, x = x,  alpha = 0)
    train_predicted <- predict(ridge.mod, s = RIDGE_lambda, newx = x)
    test_predicted <- predict(ridge.mod, s = RIDGE_lambda, newx = as.matrix(test[,-1]))
  }



  calculations <- calculate_metrics(train_predicted, test_predicted,
                                    train_observed, test_observed, digits = 15)
  list_RIDGE[[b]] <- calculations

  # Lasso Regression
  x = model.matrix(formula, train)[,-DepIndex]
  y = as.matrix(train[,DepIndex])

  if (numIND < 2){
    lasso.mod <- lm(formula, data = train)
    train_predicted <- predict(lasso.mod, train)
    test_predicted <- predict(lasso.mod, test)
  } else {
  lasso.mod <- glmnet(y = y, x = x,  alpha = 1)
  train_predicted <- predict(lasso.mod, s = LASSO_lambda, newx = x)
  test_predicted <- predict(lasso.mod, s = LASSO_lambda, newx = as.matrix(test[,-1]))
  }

  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15)
  list_LASSO[[b]] <- calculations

  # Polynomial regression
  if (polynomial_formula == ""){
    poly_model <- lm(formula, data = train)
  } else {
    poly_model <- lm(polynomial_formula, data = train)
  }

  train_predicted <- predict(poly_model, train)
  test_predicted <- predict(poly_model, test)
  calculations <- calculate_metrics(train_predicted, test_predicted,
                                     train_observed, test_observed, digits = 15)
  list_POLY[[b]] <- calculations


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

  b = 0 # place holder for saving rezults

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
    MLR <- lm(formula, data = train)
    train_predicted <- predict(MLR, train)
    test_predicted <- predict(MLR, test)
    train_observed <- train[, DepIndex]
    test_observed <- test[, DepIndex]
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_MLR[[b]] <- calculations

    #ANN Model
    capture.output(ANN <- brnn(formula, data = train, ANN_neurons = ANN_neurons, verbose = FALSE))
    train_predicted <- predict(ANN, train)
    test_predicted <- predict(ANN, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_ANN[[b]] <- calculations

    #M5 Model tree
    MT_model <- M5P(formula, data = train,
                    control = Weka_control(M = MT_M, N =  MT_N, U = MT_U,
                                           R = MT_R))
    train_predicted <- predict(MT_model, train)
    test_predicted <- predict(MT_model, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_MT[[b]] <- calculations

    #M5 Model with bagging
    BMT_model <- Bagging(formula,
                         data = train,
                         control = Weka_control(P = BMT_P, I = BMT_I,
                                                W = list("weka.classifiers.trees.M5P",
                                                         M = BMT_M, N = BMT_N,
                                                         U = BMT_U, R = BMT_R)))
    train_predicted <- predict(BMT_model, train)
    test_predicted <- predict(BMT_model, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_BMT[[b]] <- calculations

    ##Random Forest
    RegTree_Weka <- randomForest(formula = formula, data = train, RF_mtry = RF_mtry,
                                 RF_maxnodes = RF_maxnodes, RF_ntree = RF_ntree)
    train_predicted <- predict(RegTree_Weka, train)
    test_predicted <- predict(RegTree_Weka, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_RF[[b]] <- calculations

    # Ridge Regression!
    x = model.matrix(formula, train)[,-DepIndex]
    y = as.matrix(train[,DepIndex])

    if (numIND < 2){
      ridge.mod <- lm(formula, data = train)
      train_predicted <- predict(ridge.mod, train)
      test_predicted <- predict(ridge.mod, test)
    } else {
      ridge.mod <- glmnet(y = y, x = x,  alpha = 0)
      train_predicted <- predict(ridge.mod, s = RIDGE_lambda, newx = x)
      test_predicted <- predict(ridge.mod, s = RIDGE_lambda, newx = as.matrix(test[,-1]))
    }



    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_RIDGE[[b]] <- calculations

    # Lasso Regression
    x = model.matrix(formula, train)[,-DepIndex]
    y = as.matrix(train[,DepIndex])

    if (numIND < 2){
      lasso.mod <- lm(formula, data = train)
      train_predicted <- predict(lasso.mod, train)
      test_predicted <- predict(lasso.mod, test)
    } else {
      lasso.mod <- glmnet(y = y, x = x,  alpha = 1)
      train_predicted <- predict(lasso.mod, s = LASSO_lambda, newx = x)
      test_predicted <- predict(lasso.mod, s = LASSO_lambda, newx = as.matrix(test[,-1]))
    }

    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_LASSO[[b]] <- calculations

    # Polynomial regression
    if (polynomial_formula == ""){
      poly_model <- lm(formula, data = train)
    } else {
      poly_model <- lm(polynomial_formula, data = train)
    }
    train_predicted <- predict(poly_model, train)
    test_predicted <- predict(poly_model, test)
    calculations <- calculate_metrics(train_predicted, test_predicted,
                                       train_observed, test_observed, digits = 15)
    list_POLY[[b]] <- calculations

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
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_MLR <- data.frame(m)
df_MLR_bias <- df_MLR[c(13, 14), c(1: position)]
df_MLR_rank <- df_MLR[-c(13, 14), c(1: position)]
df_MLR_avg <- df_MLR[-c(13, 14), c(position + 1, position + 2)]
rownames(df_MLR_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")

listVec <- lapply(list_ANN, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_ANN <- data.frame(m)
df_ANN_bias <- df_ANN[c(13, 14), c(1: position)]
df_ANN_rank <- df_ANN[-c(13, 14), c(1: position)]
df_ANN_avg <- df_ANN[-c(13, 14), c(position + 1, position + 2)]
rownames(df_ANN_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")

listVec <- lapply(list_MT, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_MT <- data.frame(m)
df_MT_bias <- df_MT[c(13, 14), c(1: position)]
df_MT_rank <- df_MT[-c(13, 14), c(1: position)]
df_MT_avg <- df_MT[-c(13, 14), c(position + 1, position + 2)]
rownames(df_MT_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")

listVec <- lapply(list_BMT, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_BMT <- data.frame(m)
df_BMT_bias <- df_BMT[c(13, 14), c(1: position)]
df_BMT_rank <- df_BMT[-c(13, 14), c(1: position)]
df_BMT_avg <- df_BMT[-c(13, 14), c(position + 1, position + 2)]
rownames(df_BMT_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")

listVec <- lapply(list_RF, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_RF <- data.frame(m)
df_RF_bias <- df_RF[c(13, 14), c(1: position)]
df_RF_rank <- df_RF[-c(13, 14), c(1: position)]
df_RF_avg <- df_RF[-c(13, 14), c(position + 1, position + 2)]
rownames(df_RF_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")

listVec <- lapply(list_RIDGE, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_RIDGE <- data.frame(m)
df_RIDGE_bias <- df_RIDGE[c(13, 14), c(1: position)]
df_RIDGE_rank <- df_RIDGE[-c(13, 14), c(1: position)]
df_RIDGE_avg <- df_RIDGE[-c(13, 14), c(position + 1, position + 2)]
rownames(df_RIDGE_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                         "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                         "CE_cal", "CE_val")

if (numIND < 2){
  df_RIDGE_bias[df_RIDGE_bias<10000] <- NA
  df_RIDGE_rank[df_RIDGE_rank<10000] <- NA
  df_RIDGE_avg[df_RIDGE_avg<10000] <- NA
}

listVec <- lapply(list_LASSO, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_LASSO <- data.frame(m)
df_LASSO_bias <- df_LASSO[c(13, 14), c(1: position)]
df_LASSO_rank <- df_LASSO[-c(13, 14), c(1: position)]
df_LASSO_avg <- df_LASSO[-c(13, 14), c(position + 1, position + 2)]
rownames(df_LASSO_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                            "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                            "CE_cal", "CE_val")

if (numIND < 2){
  df_LASSO_bias[df_LASSO_bias<10000] <- NA
  df_LASSO_rank[df_LASSO_rank<10000] <- NA
  df_LASSO_avg[df_LASSO_avg<10000] <- NA
}

listVec <- lapply(list_POLY, c, recursive = TRUE)
m <- do.call(cbind, listVec)
averages <- apply(m, 1, mean)
std <- apply(m, 1, sd)
m <- cbind(m, averages, std)
df_POLY <- data.frame(m)
df_POLY_bias <- df_POLY[c(13, 14), c(1: position)]
df_POLY_rank <- df_POLY[-c(13, 14), c(1: position)]
df_POLY_avg <- df_POLY[-c(13, 14), c(position + 1, position + 2)]
rownames(df_POLY_avg) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                            "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                            "CE_cal", "CE_val")

if (polynomial_formula == ""){
  df_POLY_bias[df_POLY_bias<10000] <- NA
  df_POLY_rank[df_POLY_rank<10000] <- NA
  df_POLY_avg[df_POLY_avg<10000] <- NA
}


# Here, all data frames are binded together
df_all_avg <- round(cbind(df_MLR_avg, df_ANN_avg, df_MT_avg, df_BMT_avg, df_RF_avg, df_RIDGE_avg, df_LASSO_avg, df_POLY_avg), 8)

############################################################################################
# Calculation of ranks
df_all <- round(rbind(df_MLR_rank, df_ANN_rank, df_MT_rank, df_BMT_rank, df_RF_rank, df_RIDGE_rank, df_LASSO_rank, df_POLY_rank), 8)

# Now, all metrics (except bias) are extracted for calibration and validation
# data.
r_cal <- df_all[c(seq(1, 96, by = 12)), ]
r_val <- df_all[c(seq(2, 96, by = 12)), ]

RMSE_cal <- df_all[c(seq(3, 96, by = 12)), ]
RMSE_val <- df_all[c(seq(4, 96, by = 12)), ]

RSSE_cal <- df_all[c(seq(5, 96, by = 12)), ]
RSSE_val <- df_all[c(seq(6, 96, by = 12)), ]

d_cal <- df_all[c(seq(7, 96, by = 12)), ]
d_val <- df_all[c(seq(8, 96, by = 12)), ]

RE_cal <- df_all[c(seq(9, 96, by = 12)), ]
RE_val <- df_all[c(seq(10, 96, by = 12)), ]

CE_cal <- df_all[c(seq(11, 96, by = 12)), ]
CE_val <- df_all[c(seq(12, 96, by = 12)), ]


# Average rank and share of rank 1 is calculated
AVG_rank <- data.frame(rowMeans(apply(-r_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-r_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
r_cal_ranks <- cbind(AVG_rank, shareOne)
names(r_cal_ranks) <- c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-r_val, 2, rank, ties.method =  "min")))
shareOne <- data.frame(apply(apply(-r_val, 2, rank, ties.method =  "min"), 1,
                             count_ones) /  position)
r_val_ranks <- cbind(AVG_rank, shareOne)
names(r_val_ranks) <-  c("Average Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RMSE_cal, 2, rank,
                                      ties.method =  "min")))
shareOne <- data.frame(apply(apply(RMSE_cal, 2, rank, ties.method =  "min"),
                             1, count_ones) /  position)
RMSE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RMSE_cal_ranks) <-  c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RMSE_val, 2, rank,
                                      ties.method = "min")))
shareOne <- data.frame(apply(apply(RMSE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RMSE_val_ranks <- cbind(AVG_rank, shareOne)
names(RMSE_val_ranks) <-  c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RSSE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(RSSE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RSSE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RSSE_cal_ranks) <-  c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(RSSE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(RSSE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RSSE_val_ranks <- cbind(AVG_rank, shareOne)
names(RSSE_val_ranks) <-  c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-d_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-d_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
d_cal_ranks <- cbind(AVG_rank, shareOne)
names(d_cal_ranks) <-  c("Average Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-d_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-d_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
d_val_ranks <- cbind(AVG_rank, shareOne)
names(d_val_ranks) <-  c("Average Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-RE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-RE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
RE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RE_cal_ranks) <-  c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-RE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-RE_val, 2, rank, ties.method = "min"),
                             1, count_ones) /  position)
RE_val_ranks <- cbind(AVG_rank, shareOne)
names(RE_val_ranks) <-  c("Average Rank", "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-CE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-CE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
CE_cal_ranks <- cbind(AVG_rank, shareOne)
names(CE_cal_ranks) <- c("Average Rank",  "%rank_1")

AVG_rank <- data.frame(rowMeans(apply(-CE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-CE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) /  position)
CE_val_ranks <- cbind(AVG_rank, shareOne)
names(CE_val_ranks) <-  c("Average Rank",  "%rank_1")

# Results are rbinded together
ranks_together <- rbind(r_cal_ranks, r_val_ranks,
                       RMSE_cal_ranks, RMSE_val_ranks,
                       RSSE_cal_ranks, RSSE_val_ranks,
                       d_cal_ranks, d_val_ranks,
                       RE_cal_ranks, RE_val_ranks,
                       CE_cal_ranks, CE_val_ranks)

# Those variables have to be defined, solution suggest on Stackoverflow.com
ANN <- NULL
ANN_AR <- NULL
ANN_M <- NULL
ANN_S1 <- NULL
ANN_SD <- NULL
BMT <- NULL
BMT_AR <- NULL
BMT_S1 <- NULL
BMT_SD <- NULL
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
RIDGE_AR <- NULL
RIDGE_M <- NULL
RIDGE_S1 <- NULL
RIDGE_SD <- NULL
RIDGE <- NULL
LASSO_AR <- NULL
LASSO_M <- NULL
LASSO_S1 <- NULL
LASSO_SD <- NULL
LASSO <- NULL
POLY_AR <- NULL
POLY_M <- NULL
POLY_S1 <- NULL
POLY_SD <- NULL
POLY <- NULL
bias <- NULL
Method <- NULL
value <- NULL

ranks_together$Method <- c("MLR", "ANN", "MT", "BMT", "RF", "RIDGE", "LASSO", "POLY")
ranks_together$Period <- c(rep("cal", 8), rep("val", 8))
ranks_together$Metric <- c(rep("r", 16),
                            rep("RMSE", 16),
                            rep("RRSE", 16),
                            rep("d", 16),
                            rep("RE", 16),
                            rep("CE", 16))

colnames(ranks_together)[1] <- "Avg_rank"
togeter_AVG_rank <- reshape::cast(ranks_together,
                                  formula = Metric + Period ~ Method,
                                  value = c("Avg_rank"))
togeter_AVG_rank$Metric  <- factor(togeter_AVG_rank$Metric,
                                    levels = c("r", "RMSE", "RRSE", "d",
                                               "RE", "CE"))
togeter_AVG_rank <- togeter_AVG_rank[order(togeter_AVG_rank$Metric), ]
togeter_AVG_rank <- dplyr::select(togeter_AVG_rank, Metric, Period, MLR, RIDGE, LASSO, POLY,
                                  ANN, MT, BMT, RF)

colnames(ranks_together)[2] <- "Share_rank1"
together_share1 <- reshape::cast(ranks_together,
                                 formula = Metric + Period ~ Method,
                                 value = c("Share_rank1"))

together_share1$Metric  <- factor(together_share1$Metric,
                                   levels = c("r", "RMSE", "RRSE", "d",
                                              "RE", "CE"))

together_share1 <- together_share1[order(together_share1$Metric), ]
together_share1 <- dplyr::select(together_share1, Metric, Period, MLR, RIDGE, LASSO, POLY, ANN,
                                 MT, BMT, RF)

###############################################################################

colnames(df_all_avg) <- c("Mean MLR", "Std MLR", "Mean ANN", "Std ANN", "Mean MT",
                          "Std MT", "Mean BMT", "Std BMT", "Mean RF", "Std RF", "Mean RIDGE", "Std RIDGE",
                          "Mean LASSO", "Std LASSO", "Mean POLY", "Std POLY")
df_all_avg$Period <- c("cal", "val")
df_all_avg$Metric <- c("r", "r", "RMSE", "RMSE", "RRSE", "RRSE",
                            "d", "d", "RE", "RE", "CE", "CE")
row.names(df_all_avg) <- NULL

Rezults_mean_std <- dplyr::select(df_all_avg, Metric, Period, "Mean MLR", "Std MLR",
                                  "Mean RIDGE", "Std RIDGE", "Mean LASSO", "Std LASSO","Mean POLY", "Std POLY",
                                  "Mean ANN", "Std ANN",
                                  "Mean MT", "Std MT", "Mean BMT", "Std BMT", "Mean RF", "Std RF")

together_share1 <- together_share1[, -c(1,2)]
colnames(togeter_AVG_rank) <- c("Metric", "Period", "Avg rank MLR","Avg rank RIDGE", "Avg rank LASSO", "Avg rank POLY",
                                "Avg rank ANN", "Avg rank MT", "Avg rank BMT",
"Avg rank RF")
colnames(together_share1) <- c("%rank_1 MLR", "%rank_1 RIDGE","%rank_1 LASSO","%rank_1 POLY","%rank_1 ANN",
                                "%rank_1 MT", "%rank_1 BMT", "%rank_1 RF")
ranks <- cbind(togeter_AVG_rank, together_share1)
Rezults_ranks <- dplyr::select(ranks, Metric, Period,
                               "Avg rank MLR", "%rank_1 MLR",
                               "Avg rank RIDGE", "%rank_1 RIDGE",
                               "Avg rank LASSO", "%rank_1 LASSO",
                               "Avg rank POLY", "%rank_1 POLY",
                       "Avg rank ANN", "%rank_1 ANN",
                       "Avg rank MT", "%rank_1 MT",
                       "Avg rank BMT", "%rank_1 BMT",
                       "Avg rank RF", "%rank_1 RF")

##################################################################
# Here is a function to calculate bias
df_MLR_bias$Period <- c("Calibration", "Validation")
df_MLR_bias$Method <- "MLR"

df_RIDGE_bias$Period <- c("Calibration", "Validation")
df_RIDGE_bias$Method <- "RIDGE"

df_LASSO_bias$Period <- c("Calibration", "Validation")
df_LASSO_bias$Method <- "LASSO"

df_POLY_bias$Period <- c("Calibration", "Validation")
df_POLY_bias$Method <- "POLY"

df_ANN_bias$Period <- c("Calibration", "Validation")
df_ANN_bias$Method <- "ANN"

df_MT_bias$Period <- c("Calibration", "Validation")
df_MT_bias$Method <- "MT"

df_BMT_bias$Period <- c("Calibration", "Validation")
df_BMT_bias$Method <- "BMT"

df_RF_bias$Period <- c("Calibration", "Validation")
df_RF_bias$Method <- "RF"

bias_together <- rbind(df_MLR_bias,
                       df_ANN_bias,
                       df_MT_bias,
                       df_BMT_bias,
                       df_RF_bias,
                       df_RIDGE_bias,
                       df_LASSO_bias,
                       df_POLY_bias)


bias_together <- melt(bias_together, id.vars = c("Period", "Method"))

bias_together_calibration <- dplyr::filter(bias_together, Period == "Calibration")
bias_together_validation<- dplyr::filter(bias_together, Period == "Validation")

bias_together_calibration$Method <- factor(bias_together_calibration$Method,
                                 levels = c("MLR","RIDGE", "LASSO","POLY","ANN", "MT", "BMT", "RF"),
                                 ordered = TRUE)
bias_together_validation$Method <- factor(bias_together_validation$Method,
                                            levels = c("MLR", "RIDGE", "LASSO","POLY", "ANN", "MT", "BMT", "RF"),
                                           ordered = TRUE)

if (polynomial_formula == ""){
  bias_together_validation <- dplyr::filter(bias_together_validation, Method != "POLY")
  bias_together_calibration <- dplyr::filter(bias_together_calibration, Method != "POLY")
}

if (numIND < 2){
  bias_together_validation <- dplyr::filter(bias_together_validation, Method != "RIDGE" & Method != "LASSO")
  bias_together_calibration<- dplyr::filter(bias_together_calibration, Method != "RIDGE" & Method != "LASSO")
}

bias_together_calibration$value <- round(bias_together_calibration$value, round_bias_cal)

gg_object_cal <- ggplot(bias_together_calibration, aes(value)) +
  geom_density(aes(group = Method)) +
  geom_vline(xintercept = 0) +
  facet_grid(Method ~ ., scales = "free") +
  theme_bw() +
  theme(legend.position = "NONE", legend.title = element_blank(),
        text = element_text(size = 15))

bias_together_validation$value <- round(bias_together_validation$value, round_bias_val)

gg_object_val <- ggplot(bias_together_validation, aes(value)) +
  geom_density(aes(group = Method)) +
  geom_vline(xintercept = 0) +
  facet_grid(Method ~ .) +
  theme_bw() +
  theme(legend.position = "NONE", legend.title = element_blank(),
        text = element_text(size = 15))

##### Here both data frames are subset with round_df function #############

Rezults_mean_std <- round_df(Rezults_mean_std, digits = digits)
Rezults_ranks <- round_df(Rezults_ranks, digits = digits)

if (numIND < 2){
  Rezults_mean_std <- dplyr::select(Rezults_mean_std, -ends_with("RIDGE"), -ends_with("LASSO"))
  Rezults_ranks <- dplyr::select(Rezults_ranks, -ends_with("RIDGE"), -ends_with("LASSO"))
}

if (polynomial_formula == ""){
  Rezults_mean_std <- dplyr::select(Rezults_mean_std, -ends_with("POLY"))
  Rezults_ranks <- dplyr::select(Rezults_ranks, -ends_with("POLY"))
}

# Here, Calibration Validation subset is
a <- 0
b <- 0

if ("Calibration" %in% returns){
  a <- 1
}

if ("Validation" %in% returns){
  b <- 3
}

c <- a + b

# Here, all optimized parameters are saved in a data frame, which will be saved as
# a fifth elemnt of the final_list
parameters <- data.frame(
  Method = c("RIDGE", "LASSO", "ANN", "MT", "MT", "MT", "MT", "BMT", "BMT", "BMT", "BMT", "BMT", "BMT",
             "RF", "RF", "RF"),
  Parameter = c("RIDGE_lambda", "LASSO_lambda","ANN_neurons", "MT_M", "MT_N", "MT_U", "MT_R", "BMT_P", "BMT_I", "BMT_M",
                "BMT_N", "BMT_U", "BMT_R", "RF_mtry", "RF_maxnodes", "RF_ntree"),
  Value = c(round(RIDGE_lambda,2), round(LASSO_lambda,2), ANN_neurons, MT_M,
            ifelse(MT_N == 1, as.character("TRUE"), as.character("FALSE")),
            ifelse(MT_U == 1, as.character("TRUE"), as.character("FALSE")),
            ifelse(MT_R == 1, as.character("TRUE"), as.character("FALSE")), BMT_P, BMT_I, BMT_M,
            ifelse(BMT_N == 1, as.character("TRUE"), as.character("FALSE")),
            ifelse(BMT_U == 1, as.character("TRUE"), as.character("FALSE")),
            ifelse(BMT_R == 1, as.character("TRUE"), as.character("FALSE")),
            RF_mtry, RF_maxnodes, RF_ntree))

# If Calibration and Validation data should be returned, then this is our final results
if (c == 4){
  final_list <- list(mean_std <- Rezults_mean_std, ranks <- Rezults_ranks, bias_cal <- gg_object_cal,
                     bias_val <- gg_object_val, parameter_values <- parameters)
}

if (c == 1){
  Rezults_mean_std <- dplyr::filter(Rezults_mean_std, Period == "cal")
  Rezults_ranks <- dplyr::filter(Rezults_ranks, Period == "cal")
  final_list <- list(Rezults_mean_std, Rezults_ranks, gg_object_cal, parameters)
  final_list <- list(mean_std <- Rezults_mean_std, ranks <- Rezults_ranks, bias_cal <- gg_object_cal,
                     bias_val <- gg_object_val)
}

if (c == 3){
  Rezults_mean_std <- dplyr::filter(Rezults_mean_std, Period == "val")
  Rezults_ranks <- dplyr::filter(Rezults_ranks, Period == "val")
  final_list <- list(mean_std <- Rezults_mean_std, ranks <- Rezults_ranks,
                     bias_val <- gg_object_val, parameter_values <- parameters)
}

final_list

}
