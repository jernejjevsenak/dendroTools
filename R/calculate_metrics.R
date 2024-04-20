#' calculate_metrics
#'
#' Calculates performance metrics for train and test data. Calculated
#' performance metrics are correlation coefficient (r), root mean squared
#' error (RMSE), root relative squared error (RRSE), index of agreement
#' (d), reduction of error (RE), coefficient of efficiency (CE), detrended
#' efficiency (DE) and bias.
#'
#' @param train_predicted a vector indicating predicted data for training set
#' @param test_predicted a vector indicating predicted data for testing set
#' @param train_observed a vector indicating observed data for training set
#' @param test_observed a vector indicating observed data for training set
#' @param digits integer of number of digits to be displayed
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. This
#' additional argument is needed to calculate DE metrics.
#' @param test data frame with test data.
#' @return a data frame of calculated test and train metrics
#' @export
#'
#' @references
#' Briffa, K.R., Jones, P.D., Pilcher, J.R., Hughes, M.K., 1988. Reconstructing
#' summer temperatures in northern Fennoscandinavia back to A.D.1700 using tree
#' ring data from Scots Pine. Arct. Alp. Res. 20, 385-394.
#'
#' Fritts, H.C., 1976. Tree Rings and Climate. Academic Press, London 567 pp.
#'
#' Lorenz, E.N., 1956. Empirical Orthogonal Functions and Statistical Weather
#' Prediction. Massachusetts Institute of Technology, Department of Meteorology.
#'
#' Willmott, C.J., 1981. On the validation of models. Phys. Geogr. 2, 184-194.
#'
#' Witten, I.H., Frank, E., Hall, M.A., 2011. Data Mining: Practical Machine
#' Learning Tools and Techniques, 3rd ed. Morgan Kaufmann Publishers, Burlington
#' 629 pp.
#'
#' @examples
#' \donttest{
#' data(example_dataset_1)
#' test_data <- example_dataset_1[1:30, ]
#' train_data <- example_dataset_1[31:55, ]
#' lin_mod <- lm(MVA ~., data = train_data)
#' train_predicted <- predict(lin_mod, train_data)
#' test_predicted <- predict(lin_mod, test_data)
#' train_observed <- train_data[, 1]
#' test_observed <- test_data[, 1]
#' calculate_metrics(train_predicted, test_predicted, train_observed,
#' test_observed, test = test_data, formula = MVA ~.)
#'
#' test_data <- example_dataset_1[1:20, ]
#' train_data <- example_dataset_1[21:55, ]
#' library(brnn)
#' lin_mod <- brnn(MVA ~., data = train_data)
#' train_predicted <- predict(lin_mod, train_data)
#' test_predicted <- predict(lin_mod, test_data)
#' train_observed <- train_data[, 1]
#' test_observed <- test_data[, 1]
#' calculate_metrics(train_predicted, test_predicted, train_observed,
#' test_observed, test = test_data, formula = MVA ~.)
#' }

calculate_metrics <- function(train_predicted, test_predicted,
                               train_observed, test_observed, digits = 4, formula, test){

  # Calculating metrics for train (calibration data)
  train_cor <- cor(train_predicted, train_observed)
  train_RMSE <- MLmetrics::RMSE(train_predicted, train_observed)
  train_RRSE <- MLmetrics::RRSE(train_predicted, train_observed)
  train_d <- 1 - (sum((train_observed - train_predicted) ^ 2)) /
    sum((abs(train_predicted - mean(train_observed)) +
             abs(train_observed - mean(train_observed))) ^ 2)
  bias_train <- mean(train_observed) - mean(train_predicted)

  train_CE <- NA
  train_RE <- NA
  train_DE <- NA

  train_metrics <- data.frame(cbind(train_cor, train_RMSE, train_RRSE,
                                     train_d, train_RE, train_CE, train_DE, bias_train))
  row.names(train_metrics) <- c(deparse(substitute(train_predicted)))
  colnames(train_metrics) <- c("cor", "RMSE", "RRSE", "d", "RE", "CE", "DE", "bias")

  #Calculations for test (validation) data
  test_cor <- cor(test_observed, test_predicted)
  test_RMSE <- MLmetrics::RMSE(test_predicted, test_observed)
  test_d <- 1 - (sum((test_observed - test_predicted) ^ 2)) /
    sum((abs(test_predicted - mean(test_observed)) +
             abs(test_observed - mean(test_observed))) ^ 2)

  test_RRSE <- MLmetrics::RRSE(test_predicted, test_observed)
  bias_test <- mean(test_observed) - mean(test_predicted)


  test_RE <- 1 - (sum((test_observed - test_predicted) ^ 2) /
                           sum((test_observed - mean(train_observed)) ^ 2))

  test_CE <- 1 - (sum((test_observed - test_predicted) ^ 2) /
                           sum((test_observed - mean(test_observed)) ^ 2))

  ###
  DepIndex <- grep(as.character(formula[[2]]), colnames(test))
  DepName <- colnames(test)[DepIndex]

  test_1 <- data.frame(test[, DepIndex])
  colnames(test_1) <- "DE_TRICK"
  test_1$sequence <- seq(1:nrow(test))

  MLR_DE <- lm(formula = DE_TRICK~., data = test_1)

  DE_predicted <- predict(MLR_DE, test_1)

  test_DE <- 1 - (sum((test_observed - test_predicted) ^ 2) /
                    sum((test_observed - DE_predicted) ^ 2))
  ##

  test_metrics <- data.frame(cbind(test_cor, test_RMSE, test_RRSE, test_d,
                                    test_RE, test_CE, test_DE, bias_test))
  row.names(test_metrics) <- c(deparse(substitute(test_predicted)))
  colnames(test_metrics) <- c("cor", "RMSE", "RRSE", "d", "RE", "CE", "DE", "bias")

  metrics <- round(rbind(train_metrics, test_metrics),  digits)
  row.names(metrics) <- c("train", "test")

  metrics
}
