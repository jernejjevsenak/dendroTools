library(dendroTools)
library(testthat)

data(SEZ_rwl)
data(SEZ_daily_tmin)
data(SEZ_daily_tavg)
data(SEZ_daily_tmax)
data(SEZ_daily_psum)

test_that("the Sezana datasets have the documented structure", {

  expect_s3_class(SEZ_rwl, "rwl")
  expect_equal(dim(SEZ_rwl), c(140L, 41L))
  expect_equal(range(as.integer(row.names(SEZ_rwl))), c(1882L, 2021L))

  climate <- list(SEZ_daily_tmin, SEZ_daily_tavg, SEZ_daily_tmax,
                  SEZ_daily_psum)

  for (x in climate) {
    expect_equal(nrow(x), 25933L)
    expect_s3_class(x$date, "Date")
    expect_true(all(diff(x$date) == 1))
    expect_false(anyNA(x[[2]]))
  }

  expect_equal(names(SEZ_daily_psum), c("date", "p_sum"))
  expect_true(all(SEZ_daily_psum$p_sum >= 0))
})

test_that("the Sezana climate gives a complete daily water balance", {

  prec <- data_transform(SEZ_daily_psum, format = "daily")
  tmin <- data_transform(SEZ_daily_tmin, format = "daily")
  tmax <- data_transform(SEZ_daily_tmax, format = "daily")

  expect_equal(dim(prec), c(71L, 366L))

  wb <- water_balance(prec = prec, tmin = tmin, tmax = tmax, lat = 45.71)

  # only DOY 366 of non-leap years is missing
  expect_false(anyNA(as.matrix(wb)[, 1:365]))
})
