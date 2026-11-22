library(dendroTools)
library(testthat)

data(data_MVA)
data(LJ_monthly_temperatures)
data(LJ_monthly_precipitation)

prec_wide <- reshape2::dcast(LJ_monthly_precipitation, Year ~ Month,
                             value.var = "Precipitation")
prec_wide <- prec_wide[stats::complete.cases(prec_wide), ]
row.names(prec_wide) <- prec_wide$Year
prec_wide$Year <- NULL

gap <- LJ_monthly_temperatures
gap["1975", ] <- NA

test_that("sidak_window_matrix() leaves out incomplete years, not windows", {

  sidak_window_matrix <- dendroTools:::sidak_window_matrix

  set.seed(5)
  x <- matrix(stats::rnorm(60), 20, 3)
  expect_identical(sidak_window_matrix(x), x)

  # leading missing values, as after dc_method = "Ar"
  y <- x
  y[1:2, 1] <- NA
  y[1, 2] <- NA
  expect_identical(sidak_window_matrix(y), x[3:20, ])

  # a window without any value is dropped
  expect_identical(sidak_window_matrix(cbind(y, NA)), x[3:20, ])

  # a missing control value counts as missing
  ctrl <- matrix(0, 20, 3)
  ctrl[4, 3] <- NA
  expect_identical(sidak_window_matrix(y, x_control = ctrl),
                   x[c(3, 5:20), ])

  # with fewer than 10 complete years, the incomplete windows are dropped,
  # as in earlier versions
  w <- x
  w[1:15, 1] <- NA
  expect_identical(sidak_window_matrix(w), x[, 2:3])
})

test_that("sidak uses every window after Ar and with a missing year", {

  r <- suppressWarnings(
    monthly_response(data_MVA, LJ_monthly_temperatures,
                     row_names_subset = TRUE,
                     remove_insignificant_method = "sidak",
                     dc_method = "Ar", dc_difference = TRUE))
  expect_identical(r$significance_filter$n_windows_searched, 78L)

  r <- suppressWarnings(
    monthly_response(data_MVA, gap, row_names_subset = TRUE,
                     cor_na_use = "pairwise.complete.obs",
                     remove_insignificant_method = "sidak"))
  expect_true(r$significance_filter$applied)
  expect_identical(r$significance_filter$n_windows_searched, 78L)
})

test_that("cross-validation uses the years that have climate values", {

  r <- suppressWarnings(
    monthly_response(data_MVA, LJ_monthly_temperatures,
                     row_names_subset = TRUE, remove_insignificant = FALSE,
                     k = 4, dc_method = "Ar", dc_difference = TRUE))
  cv <- r$cross_validation
  expect_false(anyNA(cv[, c("cor", "RMSE", "RRSE", "d")]))
  expect_false(anyNA(cv[cv$Period == "Validation", c("RE", "CE", "DE")]))

  # the warning about the missing year is expected
  r <- suppressWarnings(
    monthly_response(data_MVA, gap, row_names_subset = TRUE,
                     remove_insignificant = FALSE, fixed_width = 3,
                     cor_na_use = "pairwise.complete.obs"))
  expect_false(anyNA(r$cross_validation[, c("cor", "RMSE", "RRSE", "d")]))
})

test_that("seascorr p values for temporal stability spend a degree of freedom on the control", {

  r <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE, k = 3)

  o <- r$optimized_return
  folds <- cut(seq_len(nrow(o)), breaks = 3, labels = FALSE)
  manual <- vapply(1:3, function(m) {
    d <- o[folds == m, ]
    pr <- as.numeric(psych::partial.r(d, c(1, 2), 3)[2])
    n <- nrow(d)
    t_stat <- pr * sqrt((n - 3) / (1 - pr ^ 2))
    2 * stats::pt(-abs(t_stat), n - 3)
  }, numeric(1))

  expect_equal(as.numeric(r$temporal_stability[[3]]), round(manual, 4))
})

test_that("missing-value warnings report only gaps as long as a window", {

  collect_warnings <- function(expr) {
    messages <- character(0)
    withCallingHandlers(expr, warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    messages
  }

  one_month <- LJ_monthly_temperatures
  one_month["1985", 5] <- NA
  three_months <- LJ_monthly_temperatures
  three_months["1985", 4:6] <- NA

  # a single missing month is absorbed by windows of three months
  w <- collect_warnings(
    monthly_response(data_MVA, one_month, row_names_subset = TRUE,
                     fixed_width = 3, remove_insignificant = FALSE))
  expect_length(w, 0)

  w <- collect_warnings(
    monthly_response(data_MVA, three_months, row_names_subset = TRUE,
                     fixed_width = 3, remove_insignificant = FALSE))
  expect_true(any(grepl("env_data contains gaps.*1985", w)))

  # the seascorr functions check the primary and the control data separately
  prec_gap <- prec_wide
  prec_gap["1985", 4:6] <- NA

  w <- collect_warnings(
    monthly_response_seascorr(data_MVA, one_month, prec_gap,
                              row_names_subset = TRUE, fixed_width = 3,
                              remove_insignificant = FALSE))
  expect_true(any(grepl("env_data_control contains gaps.*1985", w)))
  expect_false(any(grepl("env_data_primary", w)))
})

test_that("a window without any value is NA, not 0 or Inf", {

  prec_gap <- prec_wide
  prec_gap["1985", 4:6] <- NA

  r <- suppressWarnings(
    monthly_response(data_MVA, prec_gap, row_names_subset = TRUE,
                     fixed_width = 3, aggregate_function = "sum",
                     cor_na_use = "pairwise.complete.obs",
                     remove_insignificant = FALSE))
  climate <- r$aggregated_climate

  expect_true(is.na(climate["1985", "4_6"]))
  expect_identical(sum(is.na(climate)), 1L)
  # a partly missing window adds the months that are available
  expect_equal(as.numeric(climate["1985", "3_5"]), prec_wide["1985", 3])

  data(LJ_daily_temperatures)
  temperature_gap <- LJ_daily_temperatures
  temperature_gap["1985", 101:140] <- NA

  r <- suppressWarnings(
    daily_response(data_MVA, temperature_gap, row_names_subset = TRUE,
                   fixed_width = 30, skip_window_position = 10,
                   aggregate_function = "min",
                   cor_na_use = "pairwise.complete.obs",
                   remove_insignificant = FALSE))
  climate <- r$aggregated_climate

  expect_false(any(is.infinite(climate)))
  expect_true(is.na(climate["1985", "101_130"]))
  expect_equal(as.numeric(climate["1985", "91_120"]),
               min(as.numeric(LJ_daily_temperatures["1985", 91:100])))
})

test_that("gaps outside the analysed years raise no warning", {

  late_gap <- LJ_monthly_temperatures
  late_gap["2015", 4:6] <- NA
  early_gap <- LJ_monthly_temperatures
  early_gap["1985", 4:6] <- NA

  # data_MVA ends before 2015
  expect_warning(monthly_response(data_MVA, late_gap, row_names_subset = TRUE,
                                  fixed_width = 3,
                                  remove_insignificant = FALSE), NA)

  # 1985 lies outside subset_years
  expect_warning(monthly_response(data_MVA, early_gap, row_names_subset = TRUE,
                                  fixed_width = 3, subset_years = c(1990, 2010),
                                  remove_insignificant = FALSE), NA)

  expect_warning(monthly_response_seascorr(data_MVA, late_gap, prec_wide,
                                           row_names_subset = TRUE,
                                           fixed_width = 3,
                                           remove_insignificant = FALSE), NA)
})
