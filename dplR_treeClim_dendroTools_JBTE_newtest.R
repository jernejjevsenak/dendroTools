library(dplR)
library(reshape2)
library(dplyr)
library(dendroTools)
library(ggplot2)
library(reshape2)
library(treeclim)
library(lubridate)

data <- read.rwl("JB_crossdated_TER.rwl")

data_detrended <- detrend(data, method = "Spline", f = 0.5, nyrs = 30, make.plot = TRUE)

# Build master chronology, try with biweight = TRUE/FALSE
chronolongy <- chron(data_detrended, biweight = TRUE, prewhiten = TRUE) #biweight=true ir? usar mediana e =false ir? usar m?dias

JB_res <- chronolongy[c(35:nrow(chronolongy)),c(2,3),F] # para cortar no ano 1964 - 31
JB_std <- chronolongy[c(35:nrow(chronolongy)),c(1,3),F] # para cortar no ano 1964 - 31
cor(chronolongy$xxxstd, chronolongy$xxxres, use = "pairwise.complete")


#daily_temp <- read.csv("ERA5_TE_t.csv")
#daily_temp <- select(daily_temp, 1,3)
 monthly_temp <- read.csv("tmax_TE.csv")
 monthly_temp$ID <- NULL
 row.names(monthly_temp) <- monthly_temp$year
 monthly_temp$year <- NULL

 monthly_prec <- read.csv("prec_TE.csv")
 monthly_prec$ID <- NULL
 row.names(monthly_prec) <- monthly_prec$year
 monthly_prec$year <- NULL

 # monthly_temp <- data_transform(daily_temp, format = "monthly")
 # write.csv(monthly_temp, "monthly_temp_TER.csv")
 # colnames(daily_temp)[1]<- "date"
 # daily_temp <- data_transform(daily_temp, format = "daily", monthly_aggregate_function = "sum") # organizar os dados para di?rios - 366 dias e 118 anos


#daily_prec <- read.table("Psum_TER_new.txt") # CRU JRA
#daily_prec <- read.csv("ERA5_TE_p.csv")
#daily_prec <- select(daily_prec, 1,3)
#monthly_prec <- data_transform(daily_prec, format = "monthly", monthly_aggregate_function = "sum")
#write.csv(monthly_prec, "monthly_prec_TER.csv")
#colnames(daily_prec)[1] <- "date"
#daily_prec$year <- year(daily_prec$date)
#daily_prec <- mutate(daily_prec, date = ymd(date))

#daily_prec <- filter(daily_prec, date < "2020-01-01")
#daily_prec$doy <- yday(daily_prec$date)
#summary(daily_prec)

#daily_prec <- daily_prec[!duplicated(daily_prec[,1]),]

#daily_prec$date <- NULL
#daily_matrix <- dcast(year ~ doy, data = daily_prec, value.var = "p_sum")

#row.names(daily_matrix) <- daily_matrix$year
# daily_matrix$year <- NULL

#daily_prec <- daily_matrix

###############################################
############ Work with monthly data ###########
############ reshape climate data #############

#monthly_prec$year <- as.numeric(row.names(monthly_prec))
#monthly_prec <- melt(monthly_prec, id.vars = "year")
#colnames(monthly_prec) <- c("year", "month", "prec")

#monthly_temp$year <- as.numeric(row.names(monthly_temp))
#monthly_temp <- melt(monthly_temp, id.vars = "year")
#colnames(monthly_temp) <- c("year", "month", "temp")

#treeClim_data <- merge(monthly_temp, monthly_prec, by = c("year", "month"))

###########################################################################
######################### Apply treeclim ##################################

#dcc_response <- dcc(JB_res, treeClim_data)
#plot(dcc_response)

###########################################################################
#####################################################################################################
# Partial correlations
ds_T_res_seascorr <- monthly_response_seascorr(response = JB_res,
                                             env_data_primary = monthly_temp,
                                             env_data_control = monthly_prec,
                           row_names_subset = TRUE, boot = FALSE,
                           remove_insignificant = TRUE, alpha = 0.05)

plot(ds_T_res_seascorr, type = 1)
plot(ds_T_res_seascorr, type = 2)

df_temp <- data.frame(ds_T_res_seascorr$calculations)
df_temp$var <- "Temperature"

ds_P_res_seascorr <- monthly_response_seascorr(response = JB_res,
                                             env_data_primary = monthly_prec,
                                             env_data_control = monthly_temp,
                                             row_names_subset = TRUE, remove_insignificant = FALSE,
                                             #lower_limit = 21, upper_limit = 90, boot = FALSE,
                                             remove_insignificant = TRUE, alpha = 0.05)
df_prec <- data.frame(ds_P_res_seascorr$calculations)
df_prec$var <- "Precipitation"


heat_TE_parcial <- rbind(df_temp, df_prec)
heat_TE_parcial$temp_row_names <- seq(1, 12, by = 1)
heat_TE_parcial <- melt(heat_TE_parcial, id.vars = c("temp_row_names", "var"))

heat_TE_parcial <- mutate(heat_TE_parcial, var = factor(var, levels = c("Temperature", "Precipitation")))
heat_TE_parcial$Site <- "Terceira"
