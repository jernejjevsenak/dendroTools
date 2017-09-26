#' calculate_measures
#'
#' Calculates performance measures for train and test data. Calculated
#' performance measures are correlation coefficient (r), root mean squared
#' error (RMSE), root relative squared error (RSSE), index of agreement
#' (d), reduction of error (RE) coefficient of efficiency (CE) and bias.
#'
#' @param train_predicted a vector indicating predicted data for training set
#' @param test_predicted a vector indicating predicted data for testing set
#' @param train_observed a vector indicating observed data for training set
#' @param test_observed a vector indicating observed data for training set
#'
#' @return a data frame of calculated test and train measures
#' @export
#'
#' @examples
#' data(example_dataset_1)
#' test_data <- example_dataset_1[1:30, ]
#' train_data <- example_dataset_1[31:79, ]
#' lin_mod <- lm(MVA ~., data = train_data)
#' train_predicted <- predict(lin_mod, train_data)
#' test_predicted <- predict(lin_mod, test_data)
#' train_observed <- train_data[, 1]
#' test_observed <- test_data[, 1]
#' calculate_measures(train_predicted, test_predicted, train_observed,
#' test_observed)

calculate_measures <- function(train_predicted, test_predicted,
                               train_observed, test_observed){

  # Calculating measures for train (calibration data)
  train_cor <- cor(train_predicted, train_observed)
  a <- dcv::test.RE(train_observed, train_predicted)
  train_RMSE <- unname(a[[3]])
  train_RRSE <- MLmetrics::RRSE(train_predicted, train_observed)
  train_d <- 1 - (sum((train_observed - train_predicted) ^ 2)) /
    sum((abs(train_predicted - mean(train_observed)) +
             abs(train_observed - mean(train_observed))) ^ 2)
  bias_train <- mean(train_observed) - mean(train_predicted)
  train_CE <- 1 - (sum((train_observed - train_predicted) ^ 2) /
                     sum((train_observed - mean(train_observed)) ^ 2))
  train_RE <- 1 - (sum((train_observed - train_predicted) ^ 2) /
                     sum((train_observed - mean(test_observed)) ^ 2))

  train_measures <- data.frame(cbind(train_cor, train_RMSE, train_RRSE,
                                     train_d, train_RE, train_CE, bias_train))
  row.names(train_measures) <- c(deparse(substitute(train_predicted)))
  colnames(train_measures) <- c("cor", "RMSE", "RRSE", "d", "RE", "CE", "bias")

  #Calculations for test (validation) data
  test_cor <- cor(test_observed, test_predicted)
  a <- dcv::test.RE(test_observed, test_predicted)
  test_RMSE <- unname(a[[3]])
  test_d <- 1 - (sum((test_observed - test_predicted) ^ 2)) /
    sum((abs(test_predicted - mean(test_observed)) +
             abs(test_observed - mean(test_observed))) ^ 2)

  test_RRSE <- MLmetrics::RRSE(test_predicted, test_observed)
  bias_test <- mean(test_observed) - mean(test_predicted)
  test_CE <- 1 - (sum((test_observed - test_predicted) ^ 2) /
                    sum((test_observed - mean(test_observed)) ^ 2))
  test_RE <- 1 - (sum((test_observed - test_predicted) ^ 2) /
                    sum((test_observed - mean(train_observed)) ^ 2))

  test_measures <- data.frame(cbind(test_cor, test_RMSE, test_RRSE, test_d,
                                    test_RE, test_CE, bias_test))
  row.names(test_measures) <- c(deparse(substitute(test_predicted)))
  colnames(test_measures) <- c("cor", "RMSE", "RRSE", "d", "RE", "CE", "bias")

  measures <- round(rbind(train_measures, test_measures),  4)
  print(measures)
}
