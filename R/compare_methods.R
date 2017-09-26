#' compare_methods
#'
#' Calculates performance measures for train and test data of different
#' regression methods: multiple linear regression (MLR), artificial neural
#' networks with Bayesian regularization training algorithm (ANN), M5P model
#' trees (MT), model trees with bagging (BMT) and random forest of regression
#' trees (RF). Calculated performance measures are correlation coefficient,
#' root mean squared error (RMSE), root relative squared error (RSSE), index
#' of agreement (d), reduction of error (RE), coefficient of efficiency
#' (CE) nad mean bias.
#'
#' @param formula an object of class "formula" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted.
#' @param dataset a data frame with dependent and independent variables as
#' columns and (optional) years as row names.
#' @param k number of folds for cross-validation
#' @param repeats number of cross-vadlidation repeats
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
#' @return a list with two elements. Element one is a data frame with
#' calculated measures for five regression methods. Element two is a
#' ggplot object of bias for validation data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_dataset_1)
#'
#' # An example with default settings of machine learning algorithms
#' experiment_1 <- compare_methods(formula = MVA~.,
#' dataset = example_dataset_1, k = 3, repeats = 2)
#' experiment_1[[1]] # See a data frame results
#' experiment_1[[2]] # See a ggplot of mean bias for validation data
#'
#' experiment_2 <- compare_methods(formula = MVA~.,
#' dataset = example_dataset_1, k = 3, repeats = 2, neurons = 1,
#' MT_M = 4, MT_N = FALSE, MT_U = FALSE, MT_R = FALSE, BMT_P = 100,
#' BMT_I = 100, BMT_M = 4, BMT_N = FALSE, BMT_U = FALSE, BMT_R = FALSE,
#' RF_P = 100, RF_I = 100, RF_depth= 0, multiply = 5)
#' experiment_2[[1]] # See a data frame results
#' experiment_2[[2]] # See a ggplot of mean bias for validation data
#' }

compare_methods <- function(formula, dataset, k = 3, repeats = 2, neurons = 1,
                           MT_M = 4, MT_N = F, MT_U = F, MT_R = F, BMT_P = 100,
                           BMT_I = 100, BMT_M = 4, BMT_N = F, BMT_U = F,
                           BMT_R = F, RF_P = 100, RF_I = 100, RF_depth = 0,
                           multiply = 5) {

# This function is used to calculate measures r, RMSE, RRSE, d, RE, CE and bias
# for train and test data

#############################################################################
# Iter function is now used

# Here, empty lists are defined, where calculations will be stored. Empty lists
# for bias are defined separately, since bias should not be averaged. It is
# later given as density plots
list_MLR <- list()
list_ANN <- list()
list_MT <- list()
list_BMT <- list()
list_RF <- list()

# Empty lists for bias
list_MLR_bias <- list()
list_ANN_bias <- list()
list_MT_bias <- list()
list_BMT_bias <- list()
list_RF_bias <- list()

# Now, a for loop is used to calculate statistical measures with iter().
# Results are stored in a temporary_df.
for (m in 1:repeats){
  temporary_df <- iter(formula = formula, dataset = dataset, k = k,
                       neurons = neurons, MT_M = MT_M, MT_N = MT_N,
                       MT_U = MT_U, MT_R = MT_R, BMT_P = BMT_P,
                       BMT_I = BMT_I, BMT_M = BMT_M, BMT_N = BMT_N,
                       BMT_U = BMT_U, BMT_R = BMT_R, RF_P = RF_P,
                       RF_I = RF_I, RF_depth = RF_depth, multiply = m)

  # temporary_df is called and results are stored in a pre-defined lists
  # This is repeated two times, because bias goes to seperate lists
  list_MLR[[m]] <- temporary_df[[1]][k + 1]
  list_MLR_bias[[m]] <- temporary_df[[1]][k + 1]

  list_ANN[[m]] <- temporary_df[[2]][k + 1]
  list_ANN_bias[[m]] <- temporary_df[[2]][k + 1]

  list_MT[[m]] <- temporary_df[[3]][k + 1]
  list_MT_bias[[m]] <- temporary_df[[3]][k + 1]

  list_BMT[[m]] <- temporary_df[[4]][k + 1]
  list_BMT_bias[[m]] <- temporary_df[[4]][k + 1]

  list_RF[[m]] <- temporary_df[[5]][k + 1]
  list_RF_bias[[m]] <- temporary_df[[5]][k + 1]
}

# Here, lists are rearranged and measures are extracted
listVec <- lapply(list_MLR, c, recursive = TRUE)
m <- do.call(cbind, listVec)
vmesne <- apply(m, 1, mean)
m <- cbind(m, vmesne)
df_MLR <- data.frame(m)
df_MLR <- df_MLR[-c(13, 14), ]
rownames(df_MLR) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                    "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                    "CE_cal", "CE_val")
foldi <- paste("CV_", seq(1, repeats), sep = " ")
colnames(df_MLR) <- c(foldi, "Mean")

listVec <- lapply(list_ANN, c, recursive = TRUE)
m <- do.call(cbind, listVec)
vmesne <- apply(m, 1, mean)
m <- cbind(m, vmesne)
df_ANN <- data.frame(m)
df_ANN <- df_ANN[-c(13, 14), ]
rownames(df_ANN) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")
colnames(df_ANN) <- c(foldi, "Mean")


listVec <- lapply(list_MT, c, recursive = TRUE)
m <- do.call(cbind, listVec)
vmesne <- apply(m, 1, mean)
m <- cbind(m, vmesne)
df_MT <- data.frame(m)
df_MT <- df_MT[-c(13, 14), ]
rownames(df_MT) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                     "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                     "CE_cal", "CE_val")
colnames(df_MT) <- c(foldi, "Mean")


listVec <- lapply(list_BMT, c, recursive = TRUE)
m <- do.call(cbind, listVec)
vmesne <- apply(m, 1, mean)
m <- cbind(m, vmesne)
df_BMT <- data.frame(m)
df_BMT <- df_BMT[-c(13, 14), ]
rownames(df_BMT) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                      "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                      "CE_cal", "CE_val")
colnames(df_BMT) <- c(foldi, "Mean")


listVec <- lapply(list_RF, c, recursive = TRUE)
m <- do.call(cbind, listVec)
vmesne <- apply(m, 1, mean)
m <- cbind(m, vmesne)
df_RF <- data.frame(m)
df_RF <- df_RF[-c(13, 14), ]
rownames(df_RF) <- c("r_cal", "r_val", "RMSE_cal", "RMSE_val", "RSSE_cal",
                     "RSSE_val", "d_cal", "d_val", "RE_cal", "RE_val",
                     "CE_cal", "CE_val")
colnames(df_RF) <- c(foldi, "Mean")

# Here we make an empty row, so when we bind together all measures, that we
# have a clearly separated methods
temprow <- matrix(c(rep.int(NA, length(df_MLR))), nrow = 1,
                  ncol = length(df_MLR))
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(df_MLR)

# Here, all data frames are dinded together
df_all <- round(rbind(df_MLR, newrow, df_ANN, newrow, df_MT, newrow, df_BMT,
                      newrow, df_RF), 4)

# Standard deviation is calculated
df_all$sd <- apply(df_all[, c(1:repeats)], 1, sd)


# Now, all measures (except bias) are extracted for calibration and validation
# data.
r_cal <- df_all[c(seq(1, 64, by = 13)), c(1:repeats)]
r_val <- df_all[c(seq(2, 64, by = 13)), c(1:repeats)]

RMSE_cal <- df_all[c(seq(3, 64, by = 13)), c(1:repeats)]
RMSE_val <- df_all[c(seq(4, 64, by = 13)), c(1:repeats)]

RSSE_cal <- df_all[c(seq(5, 64, by = 13)), c(1:repeats)]
RSSE_val <- df_all[c(seq(6, 64, by = 13)), c(1:repeats)]

d_cal <- df_all[c(seq(7, 64, by = 13)), c(1:repeats)]
d_val <- df_all[c(seq(8, 64, by = 13)), c(1:repeats)]

RE_cal <- df_all[c(seq(9, 64, by = 13)), c(1:repeats)]
RE_val <- df_all[c(seq(10, 64, by = 13)), c(1:repeats)]

CE_cal <- df_all[c(seq(11, 64, by = 13)), c(1:repeats)]
CE_val <- df_all[c(seq(12, 64, by = 13)), c(1:repeats)]


# Average rank and share of rank 1 is calculated
AVG_rank <- data.frame(rowMeans(apply(-r_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-r_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
r_cal_ranks <- cbind(AVG_rank, shareOne)
names(r_cal_ranks) <- c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-r_val, 2, rank, ties.method =  "min")))
shareOne <- data.frame(apply(apply(-r_val, 2, rank, ties.method =  "min"), 1,
                             count_ones) / repeats)
r_val_ranks <- cbind(AVG_rank, shareOne)
names(r_val_ranks) <-  c("Average Rank",  "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(RMSE_cal, 2, rank,
                                      ties.method =  "min")))
shareOne <- data.frame(apply(apply(RMSE_cal, 2, rank, ties.method =  "min"),
                             1, count_ones) / repeats)
RMSE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RMSE_cal_ranks) <-  c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(RMSE_val, 2, rank,
                                      ties.method = "min")))
shareOne <- data.frame(apply(apply(RMSE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
RMSE_val_ranks <- cbind(AVG_rank, shareOne)
names(RMSE_val_ranks) <-  c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(RSSE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(RSSE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
RSSE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RSSE_cal_ranks) <-  c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(RSSE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(RSSE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
RSSE_val_ranks <- cbind(AVG_rank, shareOne)
names(RSSE_val_ranks) <-  c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-d_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-d_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
d_cal_ranks <- cbind(AVG_rank, shareOne)
names(d_cal_ranks) <-  c("Average Rank",  "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-d_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-d_val, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
d_val_ranks <- cbind(AVG_rank, shareOne)
names(d_val_ranks) <-  c("Average Rank",  "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-RE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-RE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
RE_cal_ranks <- cbind(AVG_rank, shareOne)
names(RE_cal_ranks) <-  c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-RE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-RE_val, 2, rank, ties.method = "min"),
                             1, count_ones) / repeats)
RE_val_ranks <- cbind(AVG_rank, shareOne)
names(RE_val_ranks) <-  c("Average Rank", "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-CE_cal, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-CE_cal, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
CE_cal_ranks <- cbind(AVG_rank, shareOne)
names(CE_cal_ranks) <- c("Average Rank",  "Share of Rank 1")

AVG_rank <- data.frame(rowMeans(apply(-CE_val, 2, rank, ties.method = "min")))
shareOne <- data.frame(apply(apply(-CE_val, 2, rank, ties.method = "min"), 1,
                             count_ones) / repeats)
CE_val_ranks <- cbind(AVG_rank, shareOne)
names(CE_val_ranks) <-  c("Average Rank",  "Share of Rank 1")

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
Measure <- NULL
Period <- NULL
RF <- NULL
RF_AR <- NULL
RF_M <- NULL
RF_S1 <- NULL
RF_SD <- NULL

bias <- NULL
method <- NULL

ranks_together$Method <- c("MLR", "ANN", "MT", "BMT", "RF")
ranks_together$Period <- c("cal", "cal", "cal", "cal", "cal", "val", "val",
                           "val", "val", "val")
ranks_together$Measure <- c("r", "r", "r", "r", "r", "r", "r", "r", "r", "r",
                           "RMSE", "RMSE", "RMSE", "RMSE", "RMSE", "RMSE",
                           "RMSE", "RMSE", "RMSE", "RMSE", "RSSE", "RSSE",
                           "RSSE", "RSSE", "RSSE", "RSSE", "RSSE", "RSSE",
                           "RSSE", "RSSE", "d", "d", "d", "d", "d", "d", "d",
                           "d", "d", "d", "RE", "RE", "RE", "RE", "RE", "RE",
                           "RE", "RE", "RE", "RE", "CE", "CE", "CE", "CE",
                           "CE", "CE", "CE", "CE", "CE", "CE")

colnames(ranks_together)[1] <- "Avg_rank"
togeter_AVG_rank <- reshape::cast(ranks_together,
                                  formula = Measure + Period ~ Method,
                                  value = c("Avg_rank"))
togeter_AVG_rank$Measure  <- factor(togeter_AVG_rank$Measure,
                                    levels = c("r", "RMSE", "RSSE", "d",
                                               "RE", "CE"))
togeter_AVG_rank <- togeter_AVG_rank[order(togeter_AVG_rank$Measure), ]
togeter_AVG_rank <- dplyr::select(togeter_AVG_rank, Measure, Period, MLR, ANN,
                                  MT, BMT, RF)

colnames(ranks_together)[2] <- "Share_rank1"
together_share1 <- reshape::cast(ranks_together,
                                 formula = Measure + Period ~ Method,
                                 value = c("Share_rank1"))

together_share1$Measure  <- factor(together_share1$Measure,
                                   levels = c("r", "RMSE", "RSSE", "d",
                                              "RE", "CE"))

together_share1 <- together_share1[order(together_share1$Measure), ]
together_share1 <- dplyr::select(together_share1, Measure, Period, MLR, ANN,
                                 MT, BMT, RF)

###############################################################################

df_means_sd <- rbind(df_MLR, df_ANN, df_MT, df_BMT, df_RF)
df_means_sd$sd <- apply(df_means_sd[, c(1:repeats)], 1, sd)
df_means_sd$Method <- c("MLR", "MLR", "MLR", "MLR", "MLR", "MLR", "MLR", "MLR",
                        "MLR", "MLR", "MLR", "MLR", "ANN", "ANN", "ANN", "ANN",
                        "ANN", "ANN", "ANN", "ANN", "ANN", "ANN", "ANN", "ANN",
                        "MT", "MT", "MT", "MT", "MT", "MT", "MT", "MT", "MT",
                        "MT", "MT", "MT", "BMT", "BMT", "BMT", "BMT", "BMT",
                        "BMT", "BMT", "BMT", "BMT", "BMT", "BMT", "BMT", "RF",
                        "RF", "RF", "RF", "RF", "RF", "RF", "RF", "RF", "RF",
                        "RF", "RF")
df_means_sd$Period <- c("cal", "val")
df_means_sd$Measure <- c("r", "r", "RMSE", "RMSE", "RSSE", "RSSE", "d", "d",
                         "RE", "RE", "CE", "CE")

together_means_sd <- reshape::cast(df_means_sd,
                                   formula = Measure + Period ~ Method,
                                   value = c("Mean"))

together_means_sd$Measure  <- factor(together_means_sd$Measure,
                                     levels = c("r", "RMSE", "RSSE", "d",
                                                "RE", "CE"))

together_means_sd <- together_means_sd[order(together_means_sd$Measure), ]
together_means <- dplyr::select(together_means_sd, Measure, Period, MLR,
                                ANN, MT, BMT, RF)

together_means_sd <- reshape::cast(df_means_sd,
                                   formula = Measure + Period ~ Method,
                                   value = c("sd"))

together_means_sd$Measure  <- factor(together_means_sd$Measure,
                                     levels = c("r", "RMSE", "RSSE", "d",
                                                "RE", "CE"))

together_means_sd <- together_means_sd[order(together_means_sd$Measure), ]
together_sd <- dplyr::select(together_means_sd, Measure, Period, MLR, ANN,
                             MT, BMT, RF)

colnames(together_means) <- c("Measure", "Period", "MLR_M", "ANN_M", "MT_M",
                              "BMT_M", "RF_M")

colnames(together_sd) <- c("Measure_SD", "Period_SD", "MLR_SD", "ANN_SD",
                           "MT_SD", "BMT_SD", "RF_SD")

colnames(togeter_AVG_rank) <- c("Measure_AR", "Period_AR", "MLR_AR", "ANN_AR",
                                "MT_AR", "BMT_AR", "RF_AR")

colnames(together_share1) <- c("Measure_S1", "Period_S1", "MLR_S1", "ANN_S1",
                               "MT_S1", "BMT_S1", "RF_S1")

TOGETHER <- cbind(together_means, together_sd, togeter_AVG_rank,
                  together_share1)
TOGETHER_MEAN_SD <- dplyr::select(TOGETHER, Measure, Period, MLR_M, MLR_SD,
                                 ANN_M, ANN_SD, MT_M, MT_SD, BMT_M, BMT_SD,
                                 RF_M, RF_SD)

TOGETHER_MEAN_SD[, -c(1, 2)] <- round(TOGETHER_MEAN_SD[, -c(1, 2)], 3)

TOGETHER_RANKS <- dplyr::select(TOGETHER, Measure, Period, MLR_AR, MLR_S1,
                               ANN_AR, ANN_S1, MT_AR, MT_S1, BMT_AR, BMT_S1,
                               RF_AR, RF_S1)

TOGETHER_RANKS[, -c(1, 2)] <- round(TOGETHER_RANKS[, -c(1, 2)], 3)

DF <- as.data.frame(cbind(TOGETHER_MEAN_SD, TOGETHER_RANKS[, -c(1, 2)]))

final_resuts <- dplyr::select(DF, Measure, Period,
       MLR_M, MLR_SD, MLR_AR, MLR_S1,
       ANN_M, ANN_SD, ANN_AR, ANN_S1,
       MT_M, MT_SD, MT_AR, MT_S1,
       BMT_M, BMT_SD, BMT_AR, BMT_S1,
       RF_M, RF_SD, RF_AR, RF_S1)

##################################################################
# Here is a function to calculate bias

listVec <- lapply(list_MLR_bias, c, recursive = TRUE)
m <- do.call(cbind, listVec)
MLR_bias <- m[14, ]
MLR_bias <- data.frame(bias = MLR_bias, method = "MLR")

listVec <- lapply(list_ANN_bias, c, recursive = TRUE)
m <- do.call(cbind, listVec)
ANN_bias <- m[14, ]
ANN_bias <- data.frame(bias = ANN_bias, method = "ANN")

listVec <- lapply(list_MT_bias, c, recursive = TRUE)
m <- do.call(cbind, listVec)
MT_bias <- m[14, ]
MT_bias <- data.frame(bias = MT_bias, method = "MT")

listVec <- lapply(list_BMT_bias, c, recursive = TRUE)
m <- do.call(cbind, listVec)
BMT_bias <- m[14, ]
BMT_bias <- data.frame(bias = BMT_bias, method = "BMT")

listVec <- lapply(list_RF_bias, c, recursive = TRUE)
m <- do.call(cbind, listVec)
RF_bias <- m[14, ]
RF_bias <- data.frame(bias = RF_bias, method = "RF")

bias_together <- rbind(MLR_bias, ANN_bias, MT_bias, BMT_bias, RF_bias)

gg_object <- ggplot(bias_together, aes(bias)) +
  geom_density(aes(group = method)) +
  geom_vline(xintercept = 0) +
  facet_grid(method ~.) +
  theme_bw() +
  theme(legend.position = "NONE", legend.title = element_blank(),
        text = element_text(size = 15),
        strip.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = pretty(bias_together$bias, n = 4)) +
  ggtitle("bias for validation data")

final_list <- list(final_resuts, gg_object)

final_list

}
