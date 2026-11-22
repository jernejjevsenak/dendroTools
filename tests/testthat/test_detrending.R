library(dendroTools)
library(testthat)

data(data_MVA)
data(LJ_monthly_temperatures)
data(LJ_monthly_precipitation)
data(LJ_daily_temperatures)
data(LJ_daily_precipitation)

detrend_climate <- dendroTools:::detrend_climate

rev_rows <- function(x) x[rev(seq_len(nrow(x))), , drop = FALSE]

prec_wide <- reshape2::dcast(LJ_monthly_precipitation, Year ~ Month,
                             value.var = "Precipitation")
prec_wide <- prec_wide[stats::complete.cases(prec_wide), ]
row.names(prec_wide) <- prec_wide$Year
prec_wide$Year <- NULL

test_that("detrend_climate() agrees with dplR::detrend.series() for positive series", {

  y <- rowSums(prec_wide[, as.character(5:8)])

  for (dc in c("Spline", "Mean", "Friedman")) {
    for (difference in c(FALSE, TRUE)) {
      ours <- detrend_climate(y, dc_method = dc, dc_difference = difference)
      dplr <- suppressWarnings(dplR::detrend.series(y, method = dc,
                                                    difference = difference,
                                                    make.plot = FALSE))
      expect_equal(ours, as.numeric(dplr))
    }
  }
})

test_that("Ar detrending agrees with dplR 1.7.8 and later", {

  skip_if(utils::packageVersion("dplR") < "1.7.8")

  # an autocorrelated positive series, so that the AR model has an order of
  # at least 1; for order 0 the old and the new scaling are the same
  set.seed(11)
  y <- 20 + as.numeric(stats::arima.sim(list(ar = 0.6), 60))
  expect_gte(stats::ar(y)$order, 1)

  for (difference in c(FALSE, TRUE)) {
    ours <- detrend_climate(y, dc_method = "Ar", dc_difference = difference)
    dplr <- suppressWarnings(dplR::detrend.series(y, method = "Ar",
                                                  difference = difference,
                                                  make.plot = FALSE))
    expect_equal(ours, as.numeric(dplr))
  }
})

test_that("detrend_climate() fits the dplR methods in chronological order", {

  set.seed(11)
  y <- 20 + as.numeric(stats::arima.sim(list(ar = 0.6), 60))
  years <- 1951:2010
  forward <- detrend_climate(y, dc_method = "Ar", dc_difference = TRUE)

  # the same series stored in descending order
  backward <- detrend_climate(rev(y), dc_method = "Ar", dc_difference = TRUE,
                              years = rev(years))
  expect_equal(rev(backward), forward)

  # SLD keeps the row order, as in earlier versions
  expect_identical(detrend_climate(rev(y), dc_method = "SLD", dc_scale = TRUE,
                                   years = rev(years)),
                   detrend_climate(rev(y), dc_method = "SLD", dc_scale = TRUE))
})

test_that("without detrending the results carry dc_method = 'none'", {

  a <- monthly_response(data_MVA, LJ_monthly_temperatures,
                        row_names_subset = TRUE, remove_insignificant = FALSE)
  b <- monthly_response(data_MVA, LJ_monthly_temperatures,
                        row_names_subset = TRUE, remove_insignificant = FALSE,
                        dc_method = "none")

  expect_identical(a$calculations, b$calculations)
  expect_identical(a$detrending, list(dc_method = "none"))
})

test_that("SLD in monthly_response() is the residual of a linear fit, scaled", {

  r <- monthly_response(data_MVA, LJ_monthly_temperatures, fixed_width = 3,
                        row_names_subset = TRUE, remove_insignificant = FALSE,
                        dc_method = "SLD")

  climate <- r$aggregated_climate
  window <- as.numeric(strsplit(colnames(climate)[4], "_")[[1]])
  x <- rowMeans(LJ_monthly_temperatures[row.names(climate),
                                        window[1]:window[2]])
  expected <- stats::residuals(stats::lm(x ~ seq_along(x)))
  expected <- as.numeric(expected / stats::sd(expected))

  expect_equal(as.numeric(climate[, 4]), expected)
  expect_equal(r$detrending, list(dc_method = "SLD", dc_scale = TRUE))
  expect_identical(row.names(r$optimized_return_all),
                   row.names(LJ_monthly_temperatures))
})

test_that("Mean detrending by differences does not change correlations", {

  a <- monthly_response(data_MVA, LJ_monthly_temperatures,
                        row_names_subset = TRUE, remove_insignificant = FALSE)
  b <- monthly_response(data_MVA, LJ_monthly_temperatures,
                        row_names_subset = TRUE, remove_insignificant = FALSE,
                        dc_method = "Mean", dc_difference = TRUE)
  expect_equal(a$calculations, b$calculations)

  a <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE)
  b <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE,
                                 dc_method = "Mean", dc_difference = TRUE)
  expect_equal(a$calculations, b$calculations)
})

test_that("the optimum matches $optimized_return after detrending", {

  r <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE,
                                 dc_method = "Spline", dc_difference = TRUE)

  o <- r$optimized_return
  recomputed <- psych::partial.r(data.frame(x = o[, 1], y = o[, 2], z = o[, 3]),
                                 c("x", "y"), "z")[2]

  expect_equal(r$calculations[which.max(abs(r$calculations))],
               as.numeric(recomputed))
})

test_that("seascorr results do not depend on the row order of the input", {

  a <- suppressWarnings(
    monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                              row_names_subset = TRUE,
                              remove_insignificant = FALSE,
                              dc_method = "Ar", dc_difference = TRUE))
  b <- suppressWarnings(
    monthly_response_seascorr(data_MVA, rev_rows(LJ_monthly_temperatures),
                              prec_wide, row_names_subset = TRUE,
                              remove_insignificant = FALSE,
                              dc_method = "Ar", dc_difference = TRUE))

  expect_equal(a$calculations, b$calculations)
  expect_equal(a$optimized_return_all, b$optimized_return_all)
})

test_that("ratio detrending stops for data with negative values", {

  expect_error(monthly_response(data_MVA, LJ_monthly_temperatures,
                                row_names_subset = TRUE, dc_method = "Spline"),
               "env_data contains negative values")

  expect_error(monthly_response_seascorr(data_MVA, prec_wide,
                                         LJ_monthly_temperatures,
                                         row_names_subset = TRUE,
                                         dc_method = "Friedman"),
               "env_data_control contains negative values")

  expect_error(suppressWarnings(
    daily_response_seascorr(data_MVA, LJ_daily_temperatures,
                            LJ_daily_precipitation,
                            tidy_env_data_control = TRUE,
                            row_names_subset = TRUE, fixed_width = 30,
                            dc_method = "Mean")),
    "env_data_primary contains negative values")

  expect_error(monthly_response(data_MVA, LJ_monthly_temperatures,
                                row_names_subset = TRUE,
                                dc_method = "ModNegExp"),
               "dc_method should be one of")
})

test_that("prewhitening switches to pairwise deletion with a warning", {

  expect_warning(
    r <- monthly_response(data_MVA, LJ_monthly_temperatures,
                          row_names_subset = TRUE,
                          remove_insignificant = FALSE,
                          dc_method = "Ar", dc_difference = TRUE),
    "pairwise.complete.obs")
  expect_true(any(is.finite(r$calculations)))

  # a window containing a dry month is NA under ratio detrending, and must
  # not stop the analysis
  r <- suppressWarnings(
    monthly_response(data_MVA, prec_wide, row_names_subset = TRUE,
                     remove_insignificant = FALSE, aggregate_function = "sum",
                     dc_method = "Ar"))
  expect_true(any(is.finite(r$calculations)))
})

test_that("daily_response_seascorr() reports the detrending settings", {

  r <- suppressWarnings(
    daily_response_seascorr(data_MVA, LJ_daily_temperatures,
                            LJ_daily_precipitation,
                            tidy_env_data_control = TRUE,
                            row_names_subset = TRUE, fixed_width = 40,
                            skip_window_position = 40,
                            remove_insignificant = FALSE,
                            dc_method = "Friedman", dc_difference = TRUE))

  expect_identical(r$detrending$primary$dc_method, "Friedman")
  expect_identical(r$detrending$primary$dc_scale, FALSE)
  expect_identical(r$detrending$control, r$detrending$primary)
  expect_true(any(is.finite(r$calculations)))
})

test_that("dc_env_data_control detrends the control windows separately", {

  r <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE,
                                 dc_method = "Spline", dc_difference = TRUE,
                                 dc_env_data_control =
                                   list(dc_difference = FALSE))

  expect_identical(r$detrending$primary$dc_difference, TRUE)
  expect_identical(r$detrending$control$dc_difference, FALSE)
  expect_identical(r$detrending$control$dc_method, "Spline")

  # one window recomputed by hand: May to July
  years <- row.names(r$aggregated_climate_primary)
  temperature <- rowMeans(LJ_monthly_temperatures[years, 5:7])
  precipitation <- rowMeans(prec_wide[years, 5:7])

  expect_equal(as.numeric(r$aggregated_climate_primary[, "5_7"]),
               detrend_climate(temperature, dc_method = "Spline",
                               dc_difference = TRUE, years = years))
  expect_equal(as.numeric(r$aggregated_climate_control[, "5_7"]),
               detrend_climate(precipitation, dc_method = "Spline",
                               dc_difference = FALSE, years = years))

  # the control windows can also be left as they are
  r <- monthly_response_seascorr(data_MVA, LJ_monthly_temperatures, prec_wide,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE,
                                 dc_method = "Spline", dc_difference = TRUE,
                                 dc_env_data_control = list(dc_method = "none"))

  expect_identical(r$detrending$control, list(dc_method = "none"))
  expect_equal(as.numeric(r$aggregated_climate_control[, "5_7"]),
               as.numeric(precipitation))
})

test_that("dc_env_data_control is checked", {

  # temperature as control: the default ratio detrending is not possible
  expect_error(monthly_response_seascorr(data_MVA, prec_wide,
                                         LJ_monthly_temperatures,
                                         row_names_subset = TRUE,
                                         dc_method = "Spline"),
               "dc_env_data_control = list\\(dc_difference = TRUE\\)")

  r <- monthly_response_seascorr(data_MVA, prec_wide, LJ_monthly_temperatures,
                                 row_names_subset = TRUE,
                                 remove_insignificant = FALSE,
                                 dc_method = "Spline",
                                 dc_env_data_control = list(dc_difference = TRUE))
  expect_identical(r$detrending$primary$dc_difference, FALSE)
  expect_identical(r$detrending$control$dc_difference, TRUE)

  expect_error(monthly_response_seascorr(data_MVA, LJ_monthly_temperatures,
                                         prec_wide, row_names_subset = TRUE,
                                         dc_env_data_control =
                                           list(dc_nyears = 20)),
               "unknown settings: dc_nyears")

  expect_error(monthly_response_seascorr(data_MVA, LJ_monthly_temperatures,
                                         prec_wide, row_names_subset = TRUE,
                                         dc_env_data_control = "none"),
               "must be NULL or a named list")

  expect_error(monthly_response_seascorr(data_MVA, LJ_monthly_temperatures,
                                         prec_wide, row_names_subset = TRUE,
                                         dc_method = "Spline",
                                         dc_difference = TRUE,
                                         dc_env_data_control = list(dc_f = 2)),
               "dc_env_data_control.*dc_f must be")
})

test_that("daily_response() fits Ar chronologically", {

  a <- suppressWarnings(
    daily_response(data_MVA, LJ_daily_temperatures, row_names_subset = TRUE,
                   fixed_width = 60, skip_window_position = 60,
                   dc_method = "Ar", dc_difference = TRUE))
  b <- suppressWarnings(
    daily_response(data_MVA, rev_rows(LJ_daily_temperatures),
                   row_names_subset = TRUE, fixed_width = 60,
                   skip_window_position = 60, dc_method = "Ar",
                   dc_difference = TRUE))

  expect_equal(a$calculations, b$calculations)
  expect_equal(a$optimized_return_all,
               b$optimized_return_all[row.names(a$optimized_return_all), ,
                                      drop = FALSE])
})

test_that("daily_response() does not stop on windows emptied by ratio detrending", {

  # many 3-day precipitation sums are zero, which ratio detrending turns into
  # an empty window
  r <- suppressWarnings(
    daily_response(data_MVA, LJ_daily_precipitation, tidy_env_data = TRUE,
                   row_names_subset = TRUE, fixed_width = 3,
                   aggregate_function = "sum", skip_window_position = 60,
                   dc_method = "Ar"))

  expect_is(r, "dmrs")
})

test_that("without detrending, integer climate data stay integer", {

  p <- reshape2::dcast(LJ_daily_precipitation, Year ~ DOY,
                       value.var = "Precipitation")
  row.names(p) <- p$Year
  p$Year <- NULL
  p[] <- lapply(p, function(z) as.integer(round(z)))

  r <- daily_response(data_MVA, p, row_names_subset = TRUE, fixed_width = 30,
                      aggregate_function = "sum", skip_window_position = 60)

  expect_true(is.integer(r$optimized_return[, 2]))
})
