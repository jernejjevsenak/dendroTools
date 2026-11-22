library(dendroTools)
library(testthat)

data(data_MVA)
data(LJ_monthly_temperatures)
data(LJ_monthly_precipitation)
data(LJ_daily_temperatures)

prec_wide <- reshape2::dcast(LJ_monthly_precipitation, Year ~ Month,
                             value.var = "Precipitation")
prec_wide <- prec_wide[stats::complete.cases(prec_wide), ]
row.names(prec_wide) <- prec_wide$Year
prec_wide$Year <- NULL

test_that("with boot = TRUE the coefficients are those of the observed data", {

  set.seed(3)

  a <- monthly_response(data_MVA, LJ_monthly_temperatures,
                        row_names_subset = TRUE, remove_insignificant = FALSE)

  # vectorised bootstrap
  b <- monthly_response(data_MVA, LJ_monthly_temperatures,
                        row_names_subset = TRUE, remove_insignificant = FALSE,
                        boot = TRUE, boot_n = 50)
  expect_equal(b$calculations, a$calculations)
  expect_true(any(is.finite(b$boot_lower)))

  # boot::boot()
  b <- suppressWarnings(
    monthly_response(data_MVA, LJ_monthly_temperatures,
                     row_names_subset = TRUE, remove_insignificant = FALSE,
                     boot = TRUE, boot_n = 50, boot_ci_type = "bca"))
  expect_equal(b$calculations, a$calculations)

  # R squared
  a <- monthly_response(data_MVA, LJ_monthly_temperatures, method = "lm",
                        row_names_subset = TRUE, remove_insignificant = FALSE)
  b <- monthly_response(data_MVA, LJ_monthly_temperatures, method = "lm",
                        row_names_subset = TRUE, remove_insignificant = FALSE,
                        boot = TRUE, boot_n = 50)
  expect_equal(b$calculations, a$calculations)

  # partial correlations
  a <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE)
  b <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE,
                                 boot = TRUE, boot_n = 50)
  expect_equal(b$calculations, a$calculations)

  b <- suppressWarnings(
    monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                              row_names_subset = TRUE,
                              remove_insignificant = FALSE,
                              boot = TRUE, boot_n = 50, boot_ci_type = "bca"))
  expect_equal(b$calculations, a$calculations)

  # daily windows of several widths
  a <- daily_response(data_MVA, LJ_daily_temperatures, row_names_subset = TRUE,
                      lower_limit = 50, upper_limit = 70,
                      skip_window_length = 10, skip_window_position = 30,
                      remove_insignificant = FALSE)
  b <- daily_response(data_MVA, LJ_daily_temperatures, row_names_subset = TRUE,
                      lower_limit = 50, upper_limit = 70,
                      skip_window_length = 10, skip_window_position = 30,
                      remove_insignificant = FALSE, boot = TRUE, boot_n = 30)
  expect_equal(b$calculations, a$calculations)
})

test_that("temporal stability of monthly_response() uses asymptotic p values", {

  # cor.test() was called with a misspelled exact = FALSE, so Kendall and
  # Spearman p values were exact whenever the data had no ties
  set.seed(7)
  years <- 1981:2010
  climate <- as.data.frame(matrix(stats::rnorm(30 * 12), 30, 12))
  row.names(climate) <- years
  response <- data.frame(proxy = climate[, 6] + stats::rnorm(30, sd = 2),
                         row.names = years)

  r <- monthly_response(response, climate, row_names_subset = TRUE,
                        cor_method = "kendall", k = 2,
                        temporal_stability_check = "progressive",
                        remove_insignificant = FALSE)

  o <- r$optimized_return
  asymptotic <- stats::cor.test(o[, 1], o[, 2], method = "kendall",
                                exact = FALSE)$p.value
  exact <- stats::cor.test(o[, 1], o[, 2], method = "kendall",
                           exact = TRUE)$p.value

  # the last step of the progressive check covers the whole period
  expect_equal(as.numeric(r$temporal_stability[2, 3]), round(asymptotic, 4))
  expect_false(isTRUE(all.equal(round(exact, 4), round(asymptotic, 4))))
})

test_that("bca with fewer replicates than years gives a warning", {

  collect_warnings <- function(expr) {
    messages <- character(0)
    withCallingHandlers(expr, warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    messages
  }

  # data_MVA and the Ljubljana temperatures share 73 years
  w <- collect_warnings(
    monthly_response(data_MVA, LJ_monthly_temperatures,
                     row_names_subset = TRUE, remove_insignificant = FALSE,
                     fixed_width = 3, boot = TRUE, boot_n = 20,
                     boot_ci_type = "bca"))
  expect_true(any(grepl("Set boot_n to at least 73", w)))

  w <- collect_warnings(
    monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                              row_names_subset = TRUE,
                              remove_insignificant = FALSE, fixed_width = 3,
                              boot = TRUE, boot_n = 20,
                              boot_ci_type = "bca"))
  expect_true(any(grepl("Set boot_n to at least", w)))

  # enough replicates
  w <- collect_warnings(
    monthly_response(data_MVA, LJ_monthly_temperatures,
                     row_names_subset = TRUE, remove_insignificant = FALSE,
                     fixed_width = 3, boot = TRUE, boot_n = 80,
                     boot_ci_type = "bca"))
  expect_false(any(grepl("boot_ci_type = 'bca' needs", w)))
})
