library(dendroTools)
library(testthat)
library(dplyr)

data(LJ_daily_temperatures)
data(example_proxies_1)

MVA_parameter <- select(example_proxies_1, MVA)
TRW_parameter <- select(example_proxies_1, TRW)

# daily_response function should return a list with matrix and two characters
test3 <- daily_response(response = example_proxies_1,
  env_data = LJ_daily_temperatures, method = "lm",
  metric = "adj.r.squared", lower = 250, upper = 251, previous_year = FALSE,
  row_names_subset = TRUE)

expect_is(test3, "dmrs")
expect_is(test3[[1]], "matrix")
expect_is(test3[[2]], "character")
expect_is(test3[[2]], "character")

# stop functions were included to prevent wrong results
expect_error(daily_response(response = TRW_parameter,
                            env_data = LJ_daily_temperatures,
                            lower = 200, upper = 270, fixed_width = -368))

expect_error(daily_response(response = example_proxies_1,
  env_data = LJ_daily_temperatures, method = "cor", lower = 250,
  upper = 270, previous_year = FALSE))

expect_error(daily_response(response = example_proxies_1,
  env_data = LJ_daily_temperatures, method = "cor", lower = 280,
  upper = 270, previous_year = FALSE))

# If row.names of env_data and response do not match, error should be given
example_proxies_1_temporal <- example_proxies_1
row.names(example_proxies_1_temporal)[1] <- "1520" # random year is assigned

# to row.name of the firest row
expect_error(daily_response(response = example_proxies_1_temporal,
                            env_data = LJ_daily_temperatures,
                            method = "lm", lower = 260, upper = 270,
                            previous_year = FALSE, remove_insignificant = FALSE))

# The order of data is unimportant, if row_names_subset is set to TRUE and
# row.names are years. In this case, data is merged based on matching
# row names.
# will be ordered data
data(example_proxies_1)
MVA_parameter <- dplyr::select(example_proxies_1, MVA)
