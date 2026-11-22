library(dendroTools)
library(testthat)

data(data_MVA)
data(LJ_monthly_temperatures)
data(LJ_monthly_precipitation)
data(LJ_daily_temperatures)
data(LJ_daily_precipitation)

# Temperature as primary and precipitation as control, or the other way round
# with every paired setting swapped.
monthly_seascorr <- function(swap = FALSE, ...) {
  if (!swap) {
    monthly_response_seascorr(response = data_MVA,
                              env_data_primary = LJ_monthly_temperatures,
                              env_data_control = LJ_monthly_precipitation,
                              tidy_env_data_primary = FALSE,
                              tidy_env_data_control = TRUE,
                              aggregate_function_env_data_control = "sum",
                              row_names_subset = TRUE, ...)
  } else {
    monthly_response_seascorr(response = data_MVA,
                              env_data_primary = LJ_monthly_precipitation,
                              env_data_control = LJ_monthly_temperatures,
                              tidy_env_data_primary = TRUE,
                              tidy_env_data_control = FALSE,
                              aggregate_function_env_data_primary = "sum",
                              row_names_subset = TRUE, ...)
  }
}

daily_seascorr <- function(swap = FALSE, ...) {
  if (!swap) {
    daily_response_seascorr(response = data_MVA,
                            env_data_primary = LJ_daily_temperatures,
                            env_data_control = LJ_daily_precipitation,
                            tidy_env_data_primary = FALSE,
                            tidy_env_data_control = TRUE,
                            aggregate_function_env_data_control = "sum",
                            row_names_subset = TRUE, ...)
  } else {
    daily_response_seascorr(response = data_MVA,
                            env_data_primary = LJ_daily_precipitation,
                            env_data_control = LJ_daily_temperatures,
                            tidy_env_data_primary = TRUE,
                            tidy_env_data_control = FALSE,
                            aggregate_function_env_data_primary = "sum",
                            row_names_subset = TRUE, ...)
  }
}

# The control surface is the primary surface of the swapped run, and adding
# it leaves the primary results unchanged.
expect_control_as_swapped <- function(fun, ...) {
  with_control <- suppressWarnings(fun(pcor_control = TRUE, seed = 7, ...))
  swapped <- suppressWarnings(fun(swap = TRUE, seed = 7, ...))
  without <- suppressWarnings(fun(seed = 7, ...))

  expect_equal(with_control$calculations_control, swapped$calculations)
  expect_equal(with_control$boot_lower_control, swapped$boot_lower)
  expect_equal(with_control$boot_upper_control, swapped$boot_upper)
  expect_equal(with_control$significance_filter_control,
               swapped$significance_filter)

  for (element in c("calculations", "boot_lower", "boot_upper",
                    "significance_filter", "optimized_return",
                    "temporal_stability")) {
    expect_identical(with_control[[element]], without[[element]])
  }

  expect_null(without$calculations_control)
  invisible(with_control)
}

test_that("pcor_control gives the swapped run in monthly_response_seascorr()", {

  for (method in c("alpha", "fdr_BY", "fdr_BH", "sidak", "maxstat",
                   "cluster")) {
    expect_control_as_swapped(monthly_seascorr,
                              remove_insignificant_method = method,
                              perm_n = 99)
  }

  # bootstrap intervals, vectorised and with boot::boot()
  expect_control_as_swapped(monthly_seascorr, fixed_width = 3, boot = TRUE,
                            boot_n = 100,
                            remove_insignificant_method = "boot")
  expect_control_as_swapped(monthly_seascorr, fixed_width = 3, boot = TRUE,
                            boot_n = 100, boot_ci_type = "bca",
                            remove_insignificant = FALSE)
})

test_that("pcor_control gives the swapped run in daily_response_seascorr()", {

  expect_control_as_swapped(daily_seascorr, fixed_width = 30,
                            skip_window_position = 10)
  expect_control_as_swapped(daily_seascorr, lower_limit = 30, upper_limit = 40,
                            skip_window_length = 5, skip_window_position = 10,
                            remove_insignificant_method = "maxstat",
                            perm_n = 99)
  expect_control_as_swapped(daily_seascorr, fixed_width = 30,
                            skip_window_position = 20, boot = TRUE,
                            boot_n = 30,
                            remove_insignificant_method = "boot")
})

test_that("variance_partitioning splits the explained variance", {

  x <- monthly_seascorr(remove_insignificant = FALSE)
  vp <- stats::setNames(x$variance_partitioning$Value,
                        x$variance_partitioning$Variable)

  d <- x$optimized_return
  colnames(d) <- c("y", "p", "c")
  r2 <- function(formula) summary(stats::lm(formula, data = d))$r.squared

  expect_equal(vp[["r_response_primary"]], stats::cor(d$y, d$p))
  expect_equal(vp[["r_squared"]], r2(y ~ p + c))
  expect_equal(vp[["unique_primary"]], r2(y ~ p + c) - r2(y ~ c))
  expect_equal(vp[["unique_control"]], r2(y ~ p + c) - r2(y ~ p))
  expect_equal(vp[["shared"]],
               r2(y ~ p + c) - vp[["unique_primary"]] - vp[["unique_control"]])

  # the primary partial correlation is the optimum of the surface
  optimum <- x$calculations[which.max(abs(x$calculations))]
  expect_equal(vp[["pcor_primary"]], optimum)

  # Spearman: the split refers to ranks
  x <- monthly_seascorr(remove_insignificant = FALSE, pcor_method = "spearman")
  vp <- stats::setNames(x$variance_partitioning$Value,
                        x$variance_partitioning$Variable)
  d <- as.data.frame(lapply(as.data.frame(x$optimized_return), rank))
  colnames(d) <- c("y", "p", "c")
  expect_equal(vp[["r_squared"]], r2(y ~ p + c))

  # Kendall: correlations only
  x <- monthly_seascorr(remove_insignificant = FALSE, pcor_method = "kendall")
  vp <- stats::setNames(x$variance_partitioning$Value,
                        x$variance_partitioning$Variable)
  expect_true(all(is.finite(vp[c("r_response_primary", "pcor_primary",
                                 "pcor_control")])))
  expect_true(all(is.na(vp[c("r_squared", "unique_primary",
                             "unique_control", "shared")])))
})

test_that("summary() and plot() report the control partial correlations", {

  x <- monthly_seascorr(remove_insignificant = FALSE, pcor_control = TRUE)
  swapped <- monthly_seascorr(swap = TRUE, remove_insignificant = FALSE)

  s <- summary(x)
  s_swapped <- summary(swapped)
  value <- function(table, variable) table$Value[table$Variable == variable]

  expect_identical(value(s, "maximal_calculated_metric_control"),
                   value(s_swapped, "maximal_calculated_metric"))
  expect_identical(value(s, "optimal_time_window_control"),
                   value(s_swapped, "optimal_time_window"))

  vp <- x$variance_partitioning
  expect_equal(as.numeric(value(s, "control_metric_in_optimal_time_window")),
               round(vp$Value[vp$Variable == "pcor_control"], 3))

  # the primary rows are those without pcor_control
  s_without <- summary(monthly_seascorr(remove_insignificant = FALSE))
  expect_identical(s[seq_len(nrow(s_without)), ], s_without)

  expect_s3_class(plot(x, type = 3), "ggplot")
  expect_error(plot(swapped, type = 3), "pcor_control = TRUE")
  expect_error(monthly_seascorr(pcor_control = "yes"),
               "pcor_control must be either TRUE or FALSE")
})
