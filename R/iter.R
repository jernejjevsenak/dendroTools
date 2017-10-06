#' iter
#'
#' splits data into desired number of folds and performs cross-validation
#' with measure calculations
#'
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param dataset  a data frame with dependent and independent variables as
#' columns and (optional) years as row names
#' @param k number of folds for cross-validation
#' @param neurons positive integer that indicates the number of neurons used
#'  for brnn method
#' @param multiply an intiger that will be used to change the seed options
#' for different repeats. set.seed(multiply*5)
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
#' @param RF_P bagSizePercent (argument for random forest)
#' @param RF_I number of iterations (argument for random forest)
#' @param RF_depth maxDepth (argument for random forest)
#'
#' @return list with data frames of calculated measures for different methods


iter <- function(formula, dataset, k, neurons, MT_M, MT_N, MT_U, MT_R,
                BMT_P, BMT_I, BMT_M, BMT_N, BMT_U, BMT_R, RF_P, RF_I, RF_depth,
                multiply) {

  # Here, idex of dependent variable is extracted and later used to locate the
  # observed values
  DepIndex <- grep(as.character(formula[[2]]), colnames(dataset))

  list_MLR <- list()
  list_ANN <- list()
  list_MT <- list()
  list_BMT <- list()
  list_RF <- list()

  foldi <- seq(1:k)
  foldi <- paste("fold_", foldi)

  #Randomly shuffle the data
  set.seed(multiply * 5)
  dataset <- dataset[sample(nrow(dataset)), ]

  #Create 10 equally size folds
  folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)


  #Perform k fold cross validation

  for (j in 1:k){
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
    calculations <- calculate_measures(train_predicted, test_predicted,
                                       train_observed, test_observed)
    list_MLR[[j]] <- calculations

    #ANN Model
    capture.output(ANN <- brnn(formula, data = train, neurons = neurons, verbose = FALSE))
    train_predicted <- predict(ANN, train)
    test_predicted <- predict(ANN, test)
    calculations <- calculate_measures(train_predicted, test_predicted,
                                       train_observed, test_observed)
    list_ANN[[j]] <- calculations

    #M5 Model tree
    MT_model <- M5P(formula, data = train,
                   control = Weka_control(M = MT_M, N =  MT_N, U = MT_U,
                                          R = MT_R))
    train_predicted <- predict(MT_model, train)
    test_predicted <- predict(MT_model, test)
    calculations <- calculate_measures(train_predicted, test_predicted,
                                       train_observed, test_observed)
    list_MT[[j]] <- calculations

    #M5 Model with bagging
    BMT_model <- Bagging(formula,
                         data = train,
                         control = Weka_control(P = BMT_P, I = BMT_I,
                                   W = list("weka.classifiers.trees.M5P",
                                   M = BMT_M, N = BMT_N,
                                   U = BMT_U, R = BMT_R)))
    train_predicted <- predict(BMT_model, train)
    test_predicted <- predict(BMT_model, test)
    calculations <- calculate_measures(train_predicted, test_predicted,
                                       train_observed, test_observed)
    list_BMT[[j]] <- calculations

    ##Regression Tree with random forest, WEKA
    # RF <- make_Weka_classifier("weka/classifiers/trees/RandomForest")
    #RegTree_Weka <- RF(formula, data = train,
     #                  control = Weka_control(P = RF_P, I = RF_I,
      #                                        depth = RF_depth))
    RegTree_Weka <- randomForest(formula = formula, data = dataset, mtry = mtry, maxnodes = 4, ntree = 200)
    train_predicted <- predict(RegTree_Weka, train)
    test_predicted <- predict(RegTree_Weka, test)
    calculations <- calculate_measures(train_predicted, test_predicted,
                                       train_observed, test_observed)
    list_RF[[j]] <- calculations

  }

  listVec <- lapply(list_MLR, c, recursive = TRUE)
  m <- do.call(cbind, listVec)
  vmesne <- apply(m, 1, mean)
  m <- cbind(m, vmesne)
  df_MLR <- data.frame(m)
  colnames(df_MLR) <- c(foldi, "Average")
  rownames(df_MLR) <- c("r_cal", "r_val", "RSSE_cal", "RMSE_cal", "RMSE_val",
                     "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                     "CE_cal", "CE_val", "bias_cal", "bias_val")

  listVec <- lapply(list_ANN, c, recursive = TRUE)
  m <- do.call(cbind, listVec)
  vmesne <- apply(m, 1, mean)
  m <- cbind(m, vmesne)
  df_ANN <- data.frame(m)
  colnames(df_ANN) <- c(foldi, "Average")
  rownames(df_ANN) <- c("r_cal", "r_val", "RSSE_cal", "RMSE_cal", "RMSE_val",
                        "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                        "CE_cal", "CE_val", "bias_cal", "bias_val")

  listVec <- lapply(list_MT, c, recursive = TRUE)
  m <- do.call(cbind, listVec)
  vmesne <- apply(m, 1, mean)
  m <- cbind(m, vmesne)
  df_MT <- data.frame(m)
  colnames(df_MT) <- c(foldi, "Average")
  rownames(df_MT) <- c("r_cal", "r_val", "RSSE_cal", "RMSE_cal", "RMSE_val",
                    "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                    "CE_cal", "CE_val", "bias_cal", "bias_val")

  listVec <- lapply(list_BMT, c, recursive = TRUE)
  m <- do.call(cbind, listVec)
  vmesne <- apply(m, 1, mean)
  m <- cbind(m, vmesne)
  df_BMT <- data.frame(m)
  colnames(df_BMT) <- c(foldi, "Average")
  rownames(df_BMT) <- c("r_cal", "r_val", "RSSE_cal", "RMSE_cal", "RMSE_val",
                       "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                       "CE_cal", "CE_val", "bias_cal", "bias_val")

  listVec <- lapply(list_RF, c, recursive = TRUE)
  m <- do.call(cbind, listVec)
  vmesne <- apply(m, 1, mean)
  m <- cbind(m, vmesne)
  df_RF <- data.frame(m)
  colnames(df_RF) <- c(foldi, "Average")
  rownames(df_RF) <- c("r_cal", "r_val", "RSSE_cal", "RMSE_cal", "RMSE_val",
                       "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                                 "CE_cal", "CE_val", "bias_cal", "bias_val")

  final <- list(df_MLR, df_ANN, df_MT, df_BMT, df_RF)

  final
}
