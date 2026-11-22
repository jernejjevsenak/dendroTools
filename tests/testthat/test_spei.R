library(dendroTools)
library(testthat)

data(data_MVA)
data(LJ_monthly_temperatures)
data(LJ_monthly_precipitation)
data(LJ_daily_temperatures)
data(LJ_daily_precipitation)

# ---------------------------------------------------------------------------
# Reference implementation. It follows the code of SPEI::spei() 1.8.1 for the
# default settings (rectangular kernel, log-logistic distribution, ub-pwm fit)
# step by step: rolling sums of a continuous series, then one fit for each
# position in the yearly cycle, with the probability weighted moments of
# TLMoments::PWM(), the L-moments of lmomco::pwm2lmom(), the parameters of
# lmom::pelglo() and the probabilities of lmom::cdfglo().
# ---------------------------------------------------------------------------

ref_pwm <- function(x, r) {
  xs <- sort(x)
  n <- length(xs)
  w <- vapply(seq_len(n), function(j) prod(j - seq_len(r)), numeric(1))
  sum(w * xs) / prod(n:(n - r))
}

ref_lmom <- function(b) {
  L <- vapply(0:2, function(r) {
    sum(vapply(0:r, function(k) {
      (-1)^(r - k) * choose(r, k) * choose(r + k, k) * b[k + 1]
    }, numeric(1)))
  }, numeric(1))
  c(L[1], L[2], L[3] / L[2])
}

ref_pelglo <- function(lm) {
  g <- -lm[3]
  if (abs(g) <= 1e-6) return(c(lm[1], lm[2], 0))
  gg <- g * pi / sin(g * pi)
  a <- lm[2] / gg
  c(lm[1] - a * (1 - gg) / g, a, g)
}

ref_cdfglo <- function(x, para) {
  if (para[3] == 0) {
    y <- (x - para[1]) / para[2]
  } else {
    y <- -1 / para[3] * log(pmax(0, 1 - para[3] * (x - para[1]) / para[2]))
  }
  1 / (1 + exp(-y))
}

ref_fit <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  xm <- x[ok]
  if (length(xm) < 4 || is.na(sd(xm)) || sd(xm) == 0) return(out)
  lm <- ref_lmom(vapply(0:2, function(r) ref_pwm(xm, r), numeric(1)))
  if (lm[2] <= 0 || abs(lm[3]) >= 1) return(out)
  out[ok] <- qnorm(ref_cdfglo(xm, ref_pelglo(lm)))
  out
}

# SPEI of a continuous series that starts at the first position of the cycle
ref_spei <- function(x, scale, frequency) {
  n <- length(x)
  acu <- rep(NA_real_, n)
  for (t in seq(scale, n)) acu[t] <- sum(x[(t - scale + 1):t])
  cycle_pos <- (seq_len(n) - 1) %% frequency + 1
  out <- rep(NA_real_, n)
  for (cc in seq_len(frequency)) {
    idx <- which(cycle_pos == cc)
    out[idx] <- ref_fit(acu[idx])
  }
  out
}

# Expected aggregated_climate of monthly_response() for aggregate_function =
# "SPEI": every window is the SPEI of the continuous monthly series at its
# last month, with the window width as scale.
expected_monthly_spei <- function(wb, window_names, years, n_previous_years) {

  wb <- wb[order(as.integer(row.names(wb))), ]
  wb_years <- as.integer(row.names(wb))
  series <- as.vector(t(as.matrix(wb)))

  out <- vapply(window_names, function(nm) {
    cols <- as.integer(strsplit(nm, "_")[[1]])
    scale <- cols[2] - cols[1] + 1
    lag_end <- n_previous_years - (cols[2] - 1) %/% 12
    end_month <- (cols[2] - 1) %% 12 + 1
    spei <- ref_spei(series, scale, 12)
    pos <- (match(as.integer(years) - lag_end, wb_years) - 1) * 12 + end_month
    spei[pos]
  }, numeric(length(years)))

  matrix(out, nrow = length(years), dimnames = list(years, window_names))
}

leap <- function(y) (y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0

# Expected SPEI of one daily window, by brute force: the window sum for every
# year of the record, skipping DOY 366 of non-leap years, then one fit.
expected_daily_spei <- function(wb, window_name, years, n_previous_years) {

  wb_years <- as.integer(row.names(wb))
  m <- as.matrix(wb)
  cols <- as.integer(strsplit(window_name, "_")[[1]])
  last_block <- (cols[2] - 1) %/% 366
  lag_end <- n_previous_years - last_block

  sums <- vapply(wb_years, function(end_year) {
    total <- 0
    observed <- 0
    for (cc in cols[1]:cols[2]) {
      year <- end_year - (last_block - (cc - 1) %/% 366)
      doy <- (cc - 1) %% 366 + 1
      if (!(year %in% wb_years)) return(NA_real_)
      v <- m[match(year, wb_years), doy]
      if (is.na(v)) {
        if (doy == 366 && !leap(year)) next
        return(NA_real_)
      }
      total <- total + v
      observed <- observed + 1
    }
    if (observed == 0) NA_real_ else total
  }, numeric(1))

  ref_fit(sums)[match(as.integer(years) - lag_end, wb_years)]
}

# ---------------------------------------------------------------------------
# Water balance from the Ljubljana data, whose daily temperatures miss a few
# days (used to test missing values). Minimum and maximum temperatures are
# derived from the mean temperatures with a range of 10 degrees.
# ---------------------------------------------------------------------------

monthly_prec <- reshape2::dcast(LJ_monthly_precipitation, Year ~ Month,
                                value.var = "Precipitation")
monthly_prec <- years_to_rownames(monthly_prec, "Year")
monthly_years <- intersect(row.names(monthly_prec),
                           row.names(LJ_monthly_temperatures))
monthly_temp <- LJ_monthly_temperatures[monthly_years, ]
monthly_wb <- water_balance(prec = monthly_prec[monthly_years, ],
                            tmin = monthly_temp - 5, tmax = monthly_temp + 5,
                            lat = 46.07)

daily_prec <- reshape2::dcast(LJ_daily_precipitation, Year ~ DOY,
                              value.var = "Precipitation")
daily_prec <- years_to_rownames(daily_prec, "Year")
daily_years <- intersect(row.names(daily_prec), row.names(LJ_daily_temperatures))
daily_temp <- LJ_daily_temperatures[daily_years, ]
daily_wb <- water_balance(prec = daily_prec[daily_years, ],
                          tmin = daily_temp - 5, tmax = daily_temp + 5,
                          lat = 46.07, tmean = daily_temp)

# A copy without missing days other than DOY 366 of non-leap years
daily_wb_complete <- daily_wb
for (i in seq_len(nrow(daily_wb_complete))) {
  z <- unlist(daily_wb_complete[i, 1:365])
  z[is.na(z)] <- mean(z, na.rm = TRUE)
  daily_wb_complete[i, 1:365] <- z
}

# ---------------------------------------------------------------------------

test_that("the SPEI fit reproduces the reference implementation", {

  set.seed(3)
  for (skew in c(-2, 0, 2)) {
    x <- rnorm(60) + skew * rexp(60)
    expect_equal(dendroTools:::spei_standardize(x), ref_fit(x),
                 tolerance = 1e-12)
  }

  x <- c(rnorm(30), NA, NA, rnorm(10))
  expect_equal(dendroTools:::spei_standardize(x), ref_fit(x),
               tolerance = 1e-12)

  # too few values or constant values
  expect_true(all(is.na(dendroTools:::spei_standardize(c(1, 2, 3, NA)))))
  expect_true(all(is.na(dendroTools:::spei_standardize(rep(5, 10)))))
})

test_that("the SPEI fit recovers a known log-logistic distribution", {

  # quantiles of a generalized logistic distribution (xi = 10, alpha = 4,
  # k = -0.2) at the plotting positions of 2000 values
  p <- (seq_len(2000) - 0.5) / 2000
  x <- 10 + 4 / -0.2 * (1 - ((1 - p) / p)^-0.2)

  expect_equal(dendroTools:::spei_standardize(x), qnorm(p), tolerance = 0.02)
})

# SPEI::spei() of a monthly series. SPEI 1.8 prints its settings, which are
# captured; earlier versions have no verbose argument.
spei_package <- function(series, scale, start_year) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(SPEI::spei(
    ts(series, frequency = 12, start = c(start_year, 1)), scale = scale)))
  as.numeric(out$fitted)
}

test_that("SPEI agrees with SPEI::spei() when the SPEI package is installed", {

  skip_if_not_installed("SPEI")

  set.seed(11)
  wb_series <- rnorm(12 * 70, mean = 20, sd = 40)
  for (scale in c(1, 3, 6, 12)) {
    expect_equal(ref_spei(wb_series, scale, 12),
                 spei_package(wb_series, scale, 1950), tolerance = 1e-10)
  }

  result <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", remove_insignificant = FALSE))

  wb_sorted <- monthly_wb[order(as.integer(row.names(monthly_wb))), ]
  spei_6 <- spei_package(as.vector(t(as.matrix(wb_sorted))), 6,
                         min(as.integer(row.names(wb_sorted))))

  years <- row.names(result$aggregated_climate)
  pos <- (match(as.integer(years), as.integer(row.names(wb_sorted))) - 1) * 12 + 8
  expect_equal(unname(result$aggregated_climate[, "3_8"]), spei_6[pos],
               tolerance = 1e-10)
})

test_that("monthly SPEI windows are the SPEI of the continuous record", {

  for (n_prev in 0:2) {

    args <- list(response = data_MVA, env_data = monthly_wb,
                 row_names_subset = TRUE, aggregate_function = "SPEI",
                 remove_insignificant = FALSE)

    if (n_prev > 0) {
      args$previous_year <- TRUE
      args$number_previous_years <- n_prev
      args$lower_limit <- 1
      args$upper_limit <- 14
    }

    result <- suppressWarnings(do.call(monthly_response, args))

    climate <- result$aggregated_climate
    expected <- expected_monthly_spei(monthly_wb, colnames(climate),
                                      row.names(climate), n_prev)

    expect_equal(unname(as.matrix(climate)), unname(expected),
                 tolerance = 1e-10, info = paste("n_prev =", n_prev))

    expected_cor <- cor(data_MVA[row.names(climate), 1], expected[, 1])
    expect_equal(result$calculations[1, 1], expected_cor, tolerance = 1e-10)
  }
})

test_that("daily SPEI windows are the standardised DOY window sums", {

  result <- suppressWarnings(daily_response(
    response = data_MVA, env_data = daily_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", fixed_width = 45,
    remove_insignificant = FALSE))

  climate <- result$aggregated_climate
  windows <- c("1_45", "100_144", "300_344", "322_366")

  for (w in windows) {
    expect_equal(unname(climate[, w]),
                 expected_daily_spei(daily_wb, w, row.names(climate), 0),
                 tolerance = 1e-10, info = w)
  }

  # windows across the turn of the year, with the previous year included
  result_prev <- suppressWarnings(daily_response(
    response = data_MVA, env_data = daily_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", fixed_width = 45, previous_year = TRUE,
    remove_insignificant = FALSE))

  climate_prev <- result_prev$aggregated_climate
  windows_prev <- c("100_144", "340_384", "366_410", "500_544", "688_732")

  for (w in windows_prev) {
    expect_equal(unname(climate_prev[, w]),
                 expected_daily_spei(daily_wb, w, row.names(climate_prev), 1),
                 tolerance = 1e-10, info = w)
  }

  # a window in the previous year is the same window of the current year,
  # one year earlier
  current <- climate[, "100_144"]
  names(current) <- row.names(climate)
  previous <- climate_prev[, "100_144"]
  names(previous) <- as.integer(row.names(climate_prev)) - 1
  common <- intersect(names(current), names(previous))
  expect_equal(unname(previous[common]), unname(current[common]))
})

test_that("SPEI does not depend on the analysed years", {

  full <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", remove_insignificant = FALSE))

  part <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", remove_insignificant = FALSE,
    subset_years = c(1960, 1990)))

  years <- row.names(part$aggregated_climate)
  expect_equal(part$aggregated_climate,
               full$aggregated_climate[years, ])

  # the optimal window over all years of env_data
  expect_equal(nrow(full$optimized_return_all), nrow(monthly_wb))
})

test_that("SPEI windows with missing days are missing, with a warning", {

  expect_warning(
    result <- daily_response(response = data_MVA, env_data = daily_wb,
                             row_names_subset = TRUE,
                             aggregate_function = "SPEI", fixed_width = 30,
                             remove_insignificant = FALSE),
    "SPEI is not calculated for windows that include a missing value.*cor_na_use was changed")

  # LJ_daily_temperatures misses DOY 31 in 1942
  climate <- result$aggregated_climate
  expect_true(is.na(climate["1942", "31_60"]))
  expect_false(is.na(climate["1942", "32_61"]))
  expect_false(anyNA(result$calculations))

  # DOY 366 of non-leap years is skipped
  expect_false(is.na(climate["1941", "337_366"]))

  expect_silent(daily_response(response = data_MVA,
                               env_data = daily_wb_complete,
                               row_names_subset = TRUE,
                               aggregate_function = "SPEI", fixed_width = 30,
                               remove_insignificant = FALSE))
})

test_that("SPEI can only be detrended by differences", {

  expect_error(monthly_response(response = data_MVA, env_data = monthly_wb,
                                row_names_subset = TRUE,
                                aggregate_function = "SPEI",
                                dc_method = "Spline"),
               "only be detrended by differences")

  detrended <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", dc_method = "Spline", dc_difference = TRUE,
    remove_insignificant = FALSE))

  raw <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", remove_insignificant = FALSE))

  years <- row.names(raw$aggregated_climate)
  expect_equal(unname(detrended$aggregated_climate[, "4_6"]),
               dendroTools:::detrend_climate(raw$aggregated_climate[, "4_6"],
                                             dc_method = "Spline",
                                             dc_difference = TRUE,
                                             years = years))

  sld <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", dc_method = "SLD",
    remove_insignificant = FALSE))
  expect_equal(unname(sld$aggregated_climate[, "4_6"]),
               dendroTools:::detrend_climate(raw$aggregated_climate[, "4_6"],
                                             dc_method = "SLD",
                                             dc_scale = TRUE,
                                             years = years))
})

test_that("SPEI works for the primary and the control variable of seascorr", {

  result <- suppressWarnings(monthly_response_seascorr(
    response = data_MVA, env_data_primary = monthly_wb,
    env_data_control = LJ_monthly_temperatures,
    aggregate_function_env_data_primary = "SPEI", row_names_subset = TRUE,
    remove_insignificant = FALSE, pcor_control = TRUE))

  spei <- suppressWarnings(monthly_response(
    response = data_MVA, env_data = monthly_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", remove_insignificant = FALSE))

  years <- row.names(result$aggregated_climate_primary)
  expect_equal(result$aggregated_climate_primary[, "5_7"],
               spei$aggregated_climate[years, "5_7"])

  d <- data.frame(x = data_MVA[years, 1],
                  y = result$aggregated_climate_primary[, "5_7"],
                  z = result$aggregated_climate_control[, "5_7"])
  expect_equal(result$calculations["3", "5"],
               psych::partial.r(d, c("x", "y"), "z")[2])

  # SPEI as the control variable, and ratio detrending of the primary
  # temperature windows is refused by dc_check_settings() as before
  swapped <- suppressWarnings(monthly_response_seascorr(
    response = data_MVA, env_data_primary = LJ_monthly_temperatures,
    env_data_control = monthly_wb,
    aggregate_function_env_data_control = "SPEI", row_names_subset = TRUE,
    remove_insignificant = FALSE))
  expect_equal(swapped$calculations, result$calculations_control)

  expect_error(monthly_response_seascorr(
    response = data_MVA, env_data_primary = LJ_monthly_temperatures,
    env_data_control = monthly_wb,
    aggregate_function_env_data_control = "SPEI", row_names_subset = TRUE,
    dc_method = "Spline", dc_difference = TRUE,
    dc_env_data_control = list(dc_difference = FALSE)),
    "env_data_control is aggregated as SPEI")

  daily <- suppressWarnings(daily_response_seascorr(
    response = data_MVA, env_data_primary = daily_wb,
    env_data_control = LJ_daily_temperatures,
    aggregate_function_env_data_primary = "SPEI", row_names_subset = TRUE,
    fixed_width = 45, remove_insignificant = FALSE))

  daily_spei <- suppressWarnings(daily_response(
    response = data_MVA, env_data = daily_wb, row_names_subset = TRUE,
    aggregate_function = "SPEI", fixed_width = 45,
    remove_insignificant = FALSE))

  years <- row.names(daily$aggregated_climate_primary)
  expect_equal(daily$aggregated_climate_primary[, "120_164"],
               daily_spei$aggregated_climate[years, "120_164"])
})

test_that("water_balance() follows the Hargreaves equation", {

  # one day by hand: DOY 172 at 46.07 degrees north
  J <- 172
  latr <- 46.07 / 57.2957795
  delta <- 0.409 * sin(0.0172 * J - 1.39)
  dr <- 1 + 0.033 * cos(0.0172 * J)
  ws <- acos(-tan(latr) * tan(delta))
  Ra <- 37.6 * dr * (ws * sin(latr) * sin(delta) +
                       cos(latr) * cos(delta) * sin(ws))
  pet <- 0.0023 * 0.408 * Ra * (20 + 17.8) * sqrt(10)

  one_year <- function(value) {
    out <- data.frame(matrix(value, nrow = 1, ncol = 366))
    row.names(out) <- "2001"
    out
  }

  wb <- water_balance(prec = one_year(0), tmin = one_year(15),
                      tmax = one_year(25), lat = 46.07)
  expect_equal(-wb[1, 172], pet)
  expect_true(pet > 4 && pet < 5.5)

  # a mean temperature, if given, replaces (tmin + tmax) / 2
  wb_mean <- water_balance(prec = one_year(0), tmin = one_year(15),
                           tmax = one_year(25), lat = 46.07,
                           tmean = one_year(22))
  expect_equal(-wb_mean[1, 172], pet * (22 + 17.8) / (20 + 17.8))

  # monthly: middle day of the month times its length, February in a leap
  # year has 29 days
  month_data <- function(value, years) {
    out <- data.frame(matrix(value, nrow = length(years), ncol = 12))
    row.names(out) <- years
    out
  }
  wb_month <- water_balance(prec = month_data(0, c(2000, 2001)),
                            tmin = month_data(0, c(2000, 2001)),
                            tmax = month_data(10, c(2000, 2001)),
                            lat = 46.07)
  pet_day <- function(J) {
    delta <- 0.409 * sin(0.0172 * J - 1.39)
    dr <- 1 + 0.033 * cos(0.0172 * J)
    ws <- acos(-tan(latr) * tan(delta))
    Ra <- 37.6 * dr * (ws * sin(latr) * sin(delta) +
                         cos(latr) * cos(delta) * sin(ws))
    0.0023 * 0.408 * Ra * (5 + 17.8) * sqrt(10)
  }
  expect_equal(-wb_month["2000", 2], 29 * pet_day(32 + 14))
  expect_equal(-wb_month["2001", 2], 28 * pet_day(32 + 13))
  expect_equal(-wb_month["2001", 7], 31 * pet_day(182 + 14))

  # polar night and polar day
  wb_polar <- water_balance(prec = one_year(0), tmin = one_year(0),
                            tmax = one_year(10), lat = 75)
  expect_equal(wb_polar[1, 355], 0)
  expect_true(wb_polar[1, 172] < 0)

  expect_error(water_balance(prec = one_year(0), tmin = month_data(0, 2001),
                             tmax = one_year(10), lat = 46),
               "same dimensions")
  expect_error(water_balance(prec = one_year(0), tmin = one_year(0),
                             tmax = one_year(10), lat = 100),
               "lat must be")
})

test_that("water_balance() agrees with SPEI::hargreaves() for monthly data", {

  skip_if_not_installed("SPEI")

  years <- 1990:1999
  set.seed(5)
  tmin <- data.frame(matrix(rnorm(120, 5, 6), nrow = 10), row.names = years)
  tmax <- tmin + runif(120, 4, 14)
  prec <- data.frame(matrix(rexp(120, 1 / 80), nrow = 10), row.names = years)

  pet <- NULL
  utils::capture.output(pet <- suppressWarnings(SPEI::hargreaves(
    Tmin = ts(as.vector(t(as.matrix(tmin))), frequency = 12, start = c(1990, 1)),
    Tmax = ts(as.vector(t(as.matrix(tmax))), frequency = 12, start = c(1990, 1)),
    lat = 46.07)))
  pet <- as.numeric(pet)

  wb <- water_balance(prec = prec, tmin = tmin, tmax = tmax, lat = 46.07)
  expect_equal(as.vector(t(as.matrix(prec - wb))), pet, tolerance = 1e-10)
})

test_that("invalid aggregate functions are still rejected", {

  expect_error(monthly_response(response = data_MVA,
                                env_data = LJ_monthly_temperatures,
                                row_names_subset = TRUE,
                                aggregate_function = "spei"),
               "SPEI")
})
