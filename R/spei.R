#' Standardised Precipitation-Evapotranspiration Index of one window
#'
#' Internal helper. Turns the accumulated climatic water balance of one window,
#' one value per year, into SPEI values. A three-parameter log-logistic
#' distribution is fitted to the values by unbiased probability weighted
#' moments, and the fitted probabilities are transformed to standard normal
#' quantiles. These are the defaults of SPEI::spei() (distribution =
#' 'log-Logistic', fit = 'ub-pwm'): the moments are those of
#' TLMoments::PWM(), the parameters those of lmom::pelglo() and the
#' probabilities those of lmom::cdfglo(). Values outside the support of the
#' fitted distribution get a probability of 0 or 1, which is moved inside by
#' .Machine$double.eps so that their SPEI is finite.
#'
#' @param x numeric vector of window sums, NA for incomplete windows
#'
#' @return numeric vector of SPEI values, NA where x is NA or where the
#' distribution cannot be fitted (fewer than four values, constant values or
#' invalid L-moments)
#' @keywords internal
#' @noRd

spei_standardize <- function(x) {

  out <- rep(NA_real_, length(x))
  ok <- which(!is.na(x))
  n <- length(ok)

  if (n < 4) {
    return(out)
  }

  v <- sort(x[ok])
  v_sd <- stats::sd(v)

  if (!is.finite(v_sd) || v_sd == 0) {
    return(out)
  }

  # Unbiased probability weighted moments b0, b1 and b2
  j <- seq_len(n)
  b0 <- sum(v) / n
  b1 <- sum((j - 1) * v) / (n * (n - 1))
  b2 <- sum((j - 1) * (j - 2) * v) / (n * (n - 1) * (n - 2))

  # L-moments
  l1 <- b0
  l2 <- 2 * b1 - b0
  l3 <- 6 * b2 - 6 * b1 + b0

  if (!is.finite(l1) || !is.finite(l2) || l2 <= 0) {
    return(out)
  }

  t3 <- l3 / l2

  if (!is.finite(t3) || abs(t3) >= 1) {
    return(out)
  }

  # Parameters of the generalized logistic distribution (Hosking 1997)
  k <- -t3

  if (abs(k) <= 1e-6) {
    k <- 0
    alpha <- l2
    xi <- l1
  } else {
    kk <- k * pi / sin(k * pi)
    alpha <- l2 / kk
    xi <- l1 - alpha * (1 - kk) / k
  }

  # Cumulative probabilities and standard normal quantiles
  y <- (x[ok] - xi) / alpha

  if (k != 0) {
    y <- -log(pmax(0, 1 - k * y)) / k
  }

  p <- 1 / (1 + exp(-y))
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)

  out[ok] <- stats::qnorm(p)
  out
}


#' Leap years
#'
#' @param years integer vector of calendar years
#' @return logical vector
#' @keywords internal
#' @noRd

spei_is_leap_year <- function(years) {
  (years %% 4 == 0 & years %% 100 != 0) | years %% 400 == 0
}


#' Reference data used to compute SPEI for any window
#'
#' Internal helper. SPEI is standardised against all years of the climate
#' data, before they are matched with the response. The SPEI of a year
#' therefore does not depend on which years are analysed, and a window is
#' standardised against the same years as in SPEI::spei() applied to the
#' continuous climate record.
#'
#' For each number of spanned years m = 0, ..., n_previous_years, a matrix
#' with m previous years is built from env, in the same way as the lagged
#' climate matrix of daily_response() and monthly_response(): rows are the
#' years for which the m previous years are available, columns run from the
#' oldest of these years to the current one. Cumulative sums over the columns
#' give every window sum in one subtraction.
#'
#' A window is incomplete, and its sum NA, if any of its values is missing.
#' The only exception is DOY 366 of a non-leap year, which does not exist and
#' is skipped. If the row names are not calendar years, a missing DOY 366 is
#' always skipped.
#'
#' @param env data frame of the climate water balance, years as rows and one
#' year of days (366) or months (12) as columns, all years of the record
#' @param n_previous_years number of previous years in the analysis layout
#' @param daily TRUE for daily data, FALSE for monthly data
#'
#' @return a list used by spei_window_values()
#' @keywords internal
#' @noRd

spei_reference <- function(env, n_previous_years = 0L, daily = TRUE) {

  env <- data.frame(env, check.names = FALSE)
  env_years_chr <- row.names(env)
  env <- as.matrix(env)
  storage.mode(env) <- "double"
  dimnames(env) <- NULL

  block_width <- ncol(env)
  env_years <- suppressWarnings(as.integer(env_years_chr))
  years_known <- !anyNA(env_years) && !any(duplicated(env_years))

  blocks <- vector("list", n_previous_years + 1L)

  for (m in seq(0L, n_previous_years)) {

    if (m == 0L) {

      mat <- env
      row_years <- env_years_chr
      block_years <- if (years_known) matrix(env_years, ncol = 1) else NULL

    } else {

      current <- env_years[vapply(env_years, function(y) {
        all((y - m):y %in% env_years)
      }, logical(1))]
      current <- sort(current)

      mat <- do.call(cbind, lapply(seq(m, 0L), function(lag_i) {
        env[match(current - lag_i, env_years), , drop = FALSE]
      }))

      row_years <- as.character(current)
      block_years <- vapply(seq(m, 0L), function(lag_i) current - lag_i,
                            numeric(length(current)))
      block_years <- matrix(block_years, nrow = length(current))
    }

    missing <- is.na(mat)

    # DOY 366 of a non-leap year does not exist and is skipped
    if (daily && block_width == 366L && length(mat) > 0) {

      for (b in seq_len(m + 1L)) {

        doy366 <- b * block_width

        skip <- if (is.null(block_years)) {
          rep(TRUE, nrow(mat))
        } else {
          !spei_is_leap_year(block_years[, b])
        }

        missing[skip, doy366] <- FALSE
      }
    }

    values <- mat
    values[is.na(values)] <- 0

    cum_values <- matrix(0, nrow = nrow(mat), ncol = ncol(mat) + 1L)
    cum_missing <- matrix(0L, nrow = nrow(mat), ncol = ncol(mat) + 1L)
    cum_observed <- matrix(0L, nrow = nrow(mat), ncol = ncol(mat) + 1L)

    for (k in seq_len(ncol(mat))) {
      cum_values[, k + 1L] <- cum_values[, k] + values[, k]
      cum_missing[, k + 1L] <- cum_missing[, k] + missing[, k]
      cum_observed[, k + 1L] <- cum_observed[, k] + !is.na(mat[, k])
    }

    blocks[[m + 1L]] <- list(years = row_years,
                             cum_values = cum_values,
                             cum_missing = cum_missing,
                             cum_observed = cum_observed)
  }

  list(block_width = block_width,
       n_previous_years = as.integer(n_previous_years),
       blocks = blocks)
}


#' SPEI of one window for selected years
#'
#' Internal helper. The window is given by its columns in the analysis layout
#' of daily_response() or monthly_response(), which runs from the oldest
#' previous year to the current year. The window is standardised as the
#' calendar window it covers: a window that lies within a previous year is
#' standardised against all years of the record, like the same window in the
#' current year, and its value for year t is the value of year t - 1 (or
#' earlier).
#'
#' @param ref list returned by spei_reference()
#' @param cols contiguous column indices of the window
#' @param years row names (years) for which values are returned
#'
#' @return numeric vector of SPEI values, one per element of years
#' @keywords internal
#' @noRd

spei_window_values <- function(ref, cols, years) {

  first_col <- min(cols)
  last_col <- max(cols)

  block_width <- ref$block_width
  first_block <- (first_col - 1L) %/% block_width
  last_block <- (last_col - 1L) %/% block_width

  # number of previous years spanned by the window and lag of its last year
  m <- last_block - first_block
  lag_end <- ref$n_previous_years - last_block

  first_ref <- first_col - first_block * block_width
  last_ref <- last_col - first_block * block_width

  blk <- ref$blocks[[m + 1L]]

  sums <- blk$cum_values[, last_ref + 1L] - blk$cum_values[, first_ref]
  n_missing <- blk$cum_missing[, last_ref + 1L] - blk$cum_missing[, first_ref]
  n_observed <- blk$cum_observed[, last_ref + 1L] - blk$cum_observed[, first_ref]
  sums[n_missing > 0 | n_observed == 0] <- NA

  spei <- spei_standardize(sums)

  keys <- if (lag_end == 0L) {
    as.character(years)
  } else {
    as.character(as.integer(years) - lag_end)
  }

  out <- spei[match(keys, blk$years)]
  names(out) <- years
  out
}


#' Years with missing values for SPEI
#'
#' Internal helper. Unlike the other aggregate functions, SPEI does not skip
#' missing values: a window with a missing day or month has no SPEI in that
#' year. Reports the calendar years whose missing values enter the analysed
#' rows, other than DOY 366 of non-leap years. In the lagged layout, a row
#' also holds the previous years, so a gap in year t - 1 is reported as
#' t - 1, not as the analysed year t.
#'
#' @param env analysed climate data, in the lagged layout
#' @param block_width 366 for daily data, 12 for monthly data
#'
#' @return a character vector of years
#' @keywords internal
#' @noRd

spei_gap_years <- function(env, block_width) {

  if (!anyNA(env)) {
    return(character(0))
  }

  missing <- is.na(as.matrix(env))
  n_blocks <- ncol(missing) %/% block_width
  years <- suppressWarnings(as.integer(row.names(env)))

  # Without calendar years there is one year per row and DOY 366 is skipped
  if (anyNA(years)) {

    if (block_width == 366L) {
      missing[, seq(366L, ncol(missing), by = 366L)] <- FALSE
    }

    return(row.names(env)[rowSums(missing) > 0])
  }

  gap_years <- integer(0)

  for (b in seq_len(n_blocks)) {

    block <- missing[, (b - 1L) * block_width + seq_len(block_width),
                     drop = FALSE]
    block_years <- years - (n_blocks - b)

    if (block_width == 366L) {
      block[!spei_is_leap_year(block_years), 366L] <- FALSE
    }

    gap_years <- c(gap_years, block_years[rowSums(block) > 0])
  }

  as.character(sort(unique(gap_years)))
}


#' Warning text for incomplete SPEI windows
#'
#' @param env_name name of the climate data frame, used in the message
#' @param gap_years years returned by spei_gap_years()
#' @param unit "days" or "months"
#' @param switched TRUE if cor_na_use was changed to 'pairwise.complete.obs'
#' @param old_na_use the previous value of cor_na_use
#' @param detrending TRUE if climate detrending is applied
#' @param na_arg name of the argument that controls missing values
#'
#' @return a character string
#' @keywords internal
#' @noRd

spei_gap_message <- function(env_name, gap_years, unit = "days",
                             switched = FALSE, old_na_use = "everything",
                             detrending = FALSE, na_arg = "cor_na_use") {

  msg <- paste0(env_name, " has missing ", unit, " in the following years: ",
                format_year_list(gap_years), ". SPEI is not calculated for ",
                "windows that include a missing value, so such windows are ",
                "missing in the analysed years that use them.")

  if (detrending) {
    msg <- paste0(msg, " Climate detrending is fitted on the years with SPEI ",
                  "values and the missing values are reinserted afterwards.")
  }

  if (switched) {
    msg <- paste0(msg, " ", na_arg, " was changed from '", old_na_use,
                  "' to 'pairwise.complete.obs' so that correlations are ",
                  "calculated over the available years.")
  }

  msg
}


#' Check detrending settings for SPEI
#'
#' Internal helper. SPEI values are centred on zero, so ratio detrending is
#' not defined for them.
#'
#' @param aggregate_function the aggregate function of the climate data
#' @param dc_method,dc_difference the detrending settings that apply to it
#' @param env_name optional name of the climate data frame, used in the
#' message
#' @param fix sentence telling how to fix the settings
#'
#' @return NULL, invisibly; stops for ratio detrending of SPEI
#' @keywords internal
#' @noRd

spei_check_detrending <- function(aggregate_function, dc_method, dc_difference,
                                  env_name = NULL,
                                  fix = "Set dc_difference = TRUE.") {

  if (identical(aggregate_function, "SPEI") && !is.null(dc_method) &&
      !identical(dc_method, "none") && !identical(dc_method, "SLD") &&
      isFALSE(dc_difference)) {

    stop(paste0(if (!is.null(env_name)) paste0(env_name, " is aggregated as ",
                                               "SPEI. "),
                "SPEI values are standardised and centred on zero, so they ",
                "can only be detrended by differences (dc_difference = TRUE), ",
                "not by ratios. ", fix), call. = FALSE)
  }

  invisible(NULL)
}


#' Climatic water balance from temperature and precipitation
#'
#' Calculates the climatic water balance, i.e. precipitation minus potential
#' evapotranspiration (P - PET), from daily or monthly precipitation and
#' minimum and maximum temperatures. The result is the input for
#' \code{aggregate_function = "SPEI"} in \code{\link{daily_response}},
#' \code{\link{monthly_response}}, \code{\link{daily_response_seascorr}} and
#' \code{\link{monthly_response_seascorr}}. PET is estimated with the
#' Hargreaves equation (Hargreaves and Samani 1985; Allen et al. 1998), using
#' the equations of \code{SPEI::hargreaves()}.
#'
#' @param prec,tmin,tmax data frames of precipitation sums (mm) and of minimum
#' and maximum temperatures (degrees Celsius), with years as row names. Daily
#' data have 366 columns, one per day of year (DOY); monthly data have 12
#' columns, January to December. This is the layout of env_data in
#' \code{\link{daily_response}} and \code{\link{monthly_response}}; daily data
#' with a date column can be converted with \code{\link{data_transform}}. All
#' data frames must have the same dimensions and row names, and for monthly
#' data the row names must be calendar years.
#' @param lat latitude of the site in decimal degrees, negative in the southern
#' hemisphere
#' @param tmean optional data frame of mean temperatures in the same layout.
#' If NULL (default), the mean temperature is (tmin + tmax) / 2, as in
#' \code{SPEI::hargreaves()}.
#'
#' @details
#' Daily PET (mm per day) is
#' \deqn{PET = 0.0023 \cdot 0.408 \cdot R_a (T_{mean} + 17.8)
#' \sqrt{T_{max} - T_{min}}}{PET = 0.0023 * 0.408 * Ra * (Tmean + 17.8) *
#' sqrt(Tmax - Tmin)}
#' where \eqn{R_a}{Ra} is the extraterrestrial radiation of the day of year
#' (MJ m-2 day-1) and 0.408 converts it into mm of evaporation per day
#' (Allen et al. 1998, equations 20, 21 and 52). For monthly data, the same
#' equation is applied to the middle day of each month and multiplied by the
#' number of days of the month, which gives the values of
#' \code{SPEI::hargreaves()} for a monthly time series. Negative temperature
#' ranges and negative PET are set to 0. Beyond the polar circles, the sunset
#' hour angle is pi during the polar day and 0 during the polar night
#' (\code{SPEI::hargreaves()} uses the largest angle of the series for the
#' polar day). Missing values are kept missing.
#'
#' @return a data frame of the climatic water balance (mm), with the
#' dimensions, row names and column names of prec
#'
#' @references
#' Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998) Crop
#' evapotranspiration - Guidelines for computing crop water requirements. FAO
#' Irrigation and drainage paper 56. FAO, Rome.
#'
#' Hargreaves, G.H., Samani, Z.A. (1985) Reference crop evapotranspiration
#' from temperature. Applied Engineering in Agriculture, 1(2), 96-99.
#'
#' Begueria, S., Vicente-Serrano, S.M. (2023) SPEI: Calculation of the
#' Standardized Precipitation-Evapotranspiration Index. R package version
#' 1.8.1.
#'
#' @export
#'
#' @examples
#' # daily climatic water balance of the black pine site Sezana
#' prec <- data_transform(SEZ_daily_psum, format = "daily")
#' tmin <- data_transform(SEZ_daily_tmin, format = "daily")
#' tmax <- data_transform(SEZ_daily_tmax, format = "daily")
#' tmean <- data_transform(SEZ_daily_tavg, format = "daily")
#'
#' wb <- water_balance(prec = prec, tmin = tmin, tmax = tmax, lat = 45.71,
#'                     tmean = tmean)
#' wb[1:5, 1:5]
#'
#' # monthly climatic water balance
#' monthly_wb <- water_balance(
#'   prec = data_transform(SEZ_daily_psum, format = "monthly",
#'                         monthly_aggregate_function = "sum"),
#'   tmin = data_transform(SEZ_daily_tmin, format = "monthly",
#'                         monthly_aggregate_function = "mean"),
#'   tmax = data_transform(SEZ_daily_tmax, format = "monthly",
#'                         monthly_aggregate_function = "mean"),
#'   lat = 45.71)
#' round(colMeans(monthly_wb), 1)

water_balance <- function(prec, tmin, tmax, lat, tmean = NULL) {

  inputs <- list(prec = prec, tmin = tmin, tmax = tmax)
  if (!is.null(tmean)) {
    inputs$tmean <- tmean
  }

  for (nm in names(inputs)) {
    if (!is.data.frame(inputs[[nm]]) && !is.matrix(inputs[[nm]])) {
      stop(nm, " must be a data frame or a matrix with years as row names.",
           call. = FALSE)
    }
  }

  inputs <- lapply(inputs, function(z) data.frame(z, check.names = FALSE))

  ref_dim <- dim(inputs$prec)
  ref_rows <- row.names(inputs$prec)

  for (nm in names(inputs)) {

    if (!identical(dim(inputs[[nm]]), ref_dim)) {
      stop("prec, tmin, tmax and tmean must have the same dimensions, but ",
           nm, " has ", nrow(inputs[[nm]]), " rows and ", ncol(inputs[[nm]]),
           " columns, while prec has ", ref_dim[1], " rows and ", ref_dim[2],
           " columns.", call. = FALSE)
    }

    if (!identical(row.names(inputs[[nm]]), ref_rows)) {
      stop("prec, tmin, tmax and tmean must have the same row names (years) ",
           "in the same order, but ", nm, " differs from prec.", call. = FALSE)
    }

    if (!all(vapply(inputs[[nm]], is.numeric, logical(1)))) {
      stop(nm, " must contain numeric values only.", call. = FALSE)
    }
  }

  if (!is.numeric(lat) || length(lat) != 1 || is.na(lat) ||
      lat < -90 || lat > 90) {
    stop("lat must be a single latitude in decimal degrees, between -90 and 90.",
         call. = FALSE)
  }

  n_col <- ref_dim[2]

  if (!(n_col %in% c(365L, 366L, 12L))) {
    stop("The data must have 366 columns (daily data, one per day of year) ",
         "or 12 columns (monthly data), but they have ", n_col, ".",
         call. = FALSE)
  }

  daily <- n_col != 12L
  n_row <- ref_dim[1]

  p <- as.matrix(inputs$prec)
  t_min <- as.matrix(inputs$tmin)
  t_max <- as.matrix(inputs$tmax)
  t_mean <- if (is.null(tmean)) (t_min + t_max) / 2 else as.matrix(inputs$tmean)

  if (daily) {

    # day of year of every column
    J <- matrix(rep(seq_len(n_col), each = n_row), nrow = n_row)
    month_length <- NULL

  } else {

    years <- suppressWarnings(as.integer(ref_rows))

    if (anyNA(years)) {
      stop("For monthly data, the row names must be calendar years, which ",
           "are needed for the lengths of the months.", call. = FALSE)
    }

    leap <- spei_is_leap_year(years)
    month_length <- matrix(rep(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
                               each = n_row), nrow = n_row)
    month_length[leap, 2] <- 29

    # middle day of each month, as in SPEI::hargreaves()
    first_day <- t(apply(month_length, 1, function(z) cumsum(c(0, z))[1:12] + 1))
    if (n_row == 1) first_day <- matrix(first_day, nrow = 1)
    J <- first_day + round(month_length / 2 - 1)
  }

  # extraterrestrial radiation, MJ m-2 day-1
  delta <- 0.409 * sin(0.0172 * J - 1.39)
  dr <- 1 + 0.033 * cos(0.0172 * J)
  latr <- lat / 57.2957795
  sset <- -tan(latr) * tan(delta)
  omegas <- acos(pmin(pmax(sset, -1), 1))
  Ra <- 37.6 * dr * (omegas * sin(latr) * sin(delta) +
                       cos(latr) * cos(delta) * sin(omegas))
  Ra[Ra < 0] <- 0

  # Hargreaves PET, mm day-1
  t_range <- pmax(0, t_max - t_min)
  pet <- 0.0023 * 0.408 * Ra * (t_mean + 17.8) * sqrt(t_range)
  pet[!is.na(pet) & pet < 0] <- 0

  if (!daily) {
    pet <- pet * month_length
  }

  out <- data.frame(p - pet, check.names = FALSE)
  row.names(out) <- ref_rows
  colnames(out) <- colnames(inputs$prec)
  out
}
