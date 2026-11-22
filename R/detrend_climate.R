#' Detrend an aggregated climate series
#'
#' Internal helper used by daily_response(), monthly_response(),
#' daily_response_seascorr() and monthly_response_seascorr() to detrend the
#' aggregated climate series of a single moving window, i.e. a vector with one
#' value per year.
#'
#' The fitted curves are the same as in dplR::detrend(): a cubic smoothing
#' spline (dplR::caps), a horizontal mean, an autoregressive model
#' (stats::ar) and the Friedman super smoother (stats::supsmu). For strictly
#' positive series the returned values are therefore identical to
#' dplR::detrend() in dplR 1.7.8 and later, which scales the prewhitened
#' series by the mean of the original series. The ARSTAN-style safeguards of
#' dplR, which replace a
#' spline or a super smoother by the series mean and set negative
#' autoregressive fits to zero whenever the fitted curve is not strictly
#' positive, are deliberately not reproduced. They are specific to ring
#' widths and would silently change the detrending method for climate
#' variables that cross zero, such as temperature in degrees Celsius.
#'
#' The dplR methods depend on the direction of time, "Ar" most of all, because
#' prewhitening explains each year by the years before it. When years are
#' supplied and the series is not in ascending order, for example after the
#' descending sort used by the seascorr functions, the series is put in
#' chronological order for fitting and returned in its original order. "SLD"
#' always keeps the order of x, as in earlier dendroTools versions; its
#' residuals are the same whether time runs forwards or backwards.
#'
#' @param x a numeric vector, the aggregated climate series of one window
#' @param dc_method detrending method, see daily_response()
#' @param dc_nyrs spline rigidity, see daily_response()
#' @param dc_f spline frequency response, see daily_response()
#' @param dc_difference ratios (FALSE) or differences (TRUE)
#' @param dc_span span of the Friedman super smoother
#' @param dc_bass smoothness of the Friedman super smoother
#' @param dc_scale logical, divide the result by its standard deviation
#' @param years optional vector with the year of each value of x, for example
#' the row names of the climate data. Used only to fit the dplR methods in
#' chronological order. NULL (default) uses the order of x.
#'
#' @return a numeric vector of the same length as x
#' @keywords internal
#' @noRd

detrend_climate <- function(x,
                            dc_method = NULL,
                            dc_nyrs = NULL,
                            dc_f = 0.5,
                            dc_difference = FALSE,
                            dc_span = "cv",
                            dc_bass = 0,
                            dc_scale = FALSE,
                            years = NULL) {

  x <- as.numeric(x)

  if (is.null(dc_method) || identical(dc_method, "none")) {
    return(x)
  }

  # Fit the dplR methods in chronological order. A series that is already in
  # ascending order is processed exactly as before.
  if (!is.null(years) && !identical(dc_method, "SLD")) {

    time_key <- suppressWarnings(as.numeric(as.character(years)))

    if (length(time_key) == length(x) && !anyNA(time_key) &&
        is.unsorted(time_key)) {

      ord <- order(time_key)
      out <- rep(NA_real_, length(x))
      out[ord] <- detrend_climate(x[ord],
                                  dc_method = dc_method,
                                  dc_nyrs = dc_nyrs,
                                  dc_f = dc_f,
                                  dc_difference = dc_difference,
                                  dc_span = dc_span,
                                  dc_bass = dc_bass,
                                  dc_scale = dc_scale)
      return(out)
    }
  }

  n <- length(x)
  out <- rep(NA_real_, n)

  # Missing values are removed before fitting and reinserted afterwards, so
  # that interior gaps are handled in the same way as in dplR::detrend().
  good <- which(!is.na(x))

  min_n <- if (identical(dc_method, "SLD")) 3L else 5L

  if (length(good) < min_n) {
    return(out)
  }

  y <- x[good]
  ny <- length(y)
  tt <- seq_len(ny)

  # Ratio detrending is undefined when the series itself is not strictly
  # positive, for example a precipitation window in which one year is
  # completely dry. Such a window is returned as NA. dplR instead recodes
  # zeros to 0.001, which produces an enormous ratio for that year.
  if (!identical(dc_method, "SLD") && isFALSE(dc_difference) && any(y <= 0)) {
    return(out)
  }

  if (dc_method == "SLD") {

    # Simple linear detrending. Always residual based, so dc_difference is not
    # relevant here. The fitted values are calculated from the coefficients
    # rather than with predict(), so that the length is correct when x
    # contains missing values.
    cf <- stats::coef(stats::lm(y ~ tt))
    out[good] <- as.numeric(y - (cf[1] + cf[2] * tt))

  } else if (dc_method == "Ar") {

    # Prewhitening. stats::ar() returns one missing residual for each lag of
    # the fitted model, so the detrended series starts with NA values.
    ar_fit <- try(stats::ar(y), silent = TRUE)

    if (inherits(ar_fit, "try-error")) {
      return(out)
    }

    ar_series <- ar_fit$resid + ar_fit$x.mean

    # As in dplR 1.7.8 and later (OpenDendro/dplR issue 22), the prewhitened
    # series is centred or scaled by the mean of the original series rather
    # than by its own mean. Ratios are only reached for positive series, so
    # this mean is positive.
    y_mean <- mean(y)

    if (isTRUE(dc_difference)) {
      out[good] <- as.numeric(ar_series - y_mean)
    } else {
      out[good] <- as.numeric(ar_series / y_mean)
    }

  } else {

    if (dc_method == "Spline") {

      nyrs2 <- if (is.null(dc_nyrs)) floor(ny * 0.67) else dc_nyrs
      curve <- try(dplR::caps(y = y, nyrs = nyrs2, f = dc_f), silent = TRUE)

    } else if (dc_method == "Mean") {

      curve <- rep.int(mean(y), ny)

    } else if (dc_method == "Friedman") {

      curve <- try(stats::supsmu(x = tt, y = y, span = dc_span,
                                 periodic = FALSE, bass = dc_bass)[["y"]],
                   silent = TRUE)

    } else {

      stop(paste0("Unknown dc_method: ", dc_method))

    }

    if (inherits(curve, "try-error")) {
      return(out)
    }

    # A spline or a super smoother can undershoot below zero even when the
    # series is positive. dplR replaces the whole curve by the series mean in
    # that case, silently changing the method. The window is returned as NA
    # instead.
    if (isFALSE(dc_difference) && any(!is.finite(curve) | curve <= 0)) {
      return(out)
    }

    if (isTRUE(dc_difference)) {
      out[good] <- as.numeric(y - curve)
    } else {
      out[good] <- as.numeric(y / curve)
    }

  }

  if (isTRUE(dc_scale)) {

    sd_out <- stats::sd(out, na.rm = TRUE)

    if (is.finite(sd_out) && sd_out > 0) {
      out <- out / sd_out
    }
  }

  out
}


#' Check the climate detrending arguments
#'
#' Internal helper shared by daily_response(), monthly_response(),
#' daily_response_seascorr() and monthly_response_seascorr(), so that the dc_*
#' arguments are checked, and their defaults resolved, identically in all four
#' functions.
#'
#' @param dc_method,dc_nyrs,dc_f,dc_difference,dc_span,dc_bass,dc_scale the
#' detrending arguments, see daily_response()
#' @param env_list a named list of the climate data frames that will be
#' detrended. Ratio detrending (dc_difference = FALSE) stops if any of them
#' contains negative values, and the error names the data frame.
#' @param ratio_hint optional sentence appended to that error
#'
#' @return a list with four elements: dc_method (NULL when no detrending is
#' applied), dc_scale (the resolved logical, or NULL when no detrending is
#' applied), summary, the list returned to the user as $detrending, and args,
#' all resolved settings, for dc_make_apply()
#' @keywords internal
#' @noRd

dc_check_settings <- function(dc_method = NULL,
                              dc_nyrs = NULL,
                              dc_f = 0.5,
                              dc_difference = FALSE,
                              dc_span = "cv",
                              dc_bass = 0,
                              dc_scale = NULL,
                              env_list = list(),
                              ratio_hint = NULL) {

  dc_known_methods <- c("none", "SLD", "Spline", "Mean", "Ar", "Friedman")

  if (identical(dc_method, "none")) {
    dc_method <- NULL
  }

  if (!is.null(dc_method)) {

    if (!is.character(dc_method) || length(dc_method) != 1 ||
        is.na(dc_method)) {
      stop("dc_method must be a single character string or NULL.",
           call. = FALSE)
    }

    if (!(dc_method %in% dc_known_methods)) {
      stop(paste0("dc_method should be one of ",
                  paste(dc_known_methods, collapse = ", "),
                  " but instead it is: ", dc_method), call. = FALSE)
    }

    if (!is.logical(dc_difference) || length(dc_difference) != 1 ||
        is.na(dc_difference)) {
      stop("dc_difference must be either TRUE or FALSE.", call. = FALSE)
    }

    if (!is.null(dc_nyrs) &&
        (!is.numeric(dc_nyrs) || length(dc_nyrs) != 1 || is.na(dc_nyrs) ||
         dc_nyrs < 3)) {
      stop("dc_nyrs must be NULL or a single numeric value of at least 3.",
           call. = FALSE)
    }

    if (!is.numeric(dc_f) || length(dc_f) != 1 || is.na(dc_f) ||
        dc_f <= 0 || dc_f >= 1) {
      stop("dc_f must be a single numeric value between 0 and 1.",
           call. = FALSE)
    }

    if (!identical(dc_span, "cv") &&
        (!is.numeric(dc_span) || length(dc_span) != 1 || is.na(dc_span) ||
         dc_span <= 0 || dc_span > 1)) {
      stop("dc_span must be 'cv' or a single numeric value between 0 and 1.",
           call. = FALSE)
    }

    if (!is.numeric(dc_bass) || length(dc_bass) != 1 || is.na(dc_bass) ||
        dc_bass < 0 || dc_bass > 10) {
      stop("dc_bass must be a single numeric value between 0 and 10.",
           call. = FALSE)
    }

    # dc_scale = NULL means "use the convention of the selected method". SLD
    # has always been rescaled by the standard deviation of the residuals, so
    # that default is kept for backward compatibility. The dplR methods are
    # not rescaled, exactly as in dplR::detrend().
    if (is.null(dc_scale)) {
      dc_scale <- identical(dc_method, "SLD")
    }

    if (!is.logical(dc_scale) || length(dc_scale) != 1 || is.na(dc_scale)) {
      stop("dc_scale must be TRUE, FALSE or NULL.", call. = FALSE)
    }

    # Ratio detrending divides the climate series by the fitted curve and is
    # only defined for strictly positive values. Variables that cross zero,
    # such as temperature in degrees Celsius, are rejected up front, because
    # the ratio is meaningless for them and explodes near zero.
    # Non-negative variables with occasional zeros, such as precipitation, are
    # allowed here: aggregation over a window normally removes the zeros, and
    # the few windows where it does not are returned as NA by
    # detrend_climate().
    if (dc_method != "SLD" && isFALSE(dc_difference)) {

      for (env_name in names(env_list)) {

        env_min <- suppressWarnings(min(as.matrix(env_list[[env_name]]),
                                        na.rm = TRUE))

        if (!is.finite(env_min) || env_min < 0) {
          stop(paste0("dc_difference = FALSE performs ratio detrending, which ",
                      "is only defined for positive climate data, but ",
                      env_name, " contains negative values (minimum: ",
                      format(env_min), "). Set dc_difference = TRUE to ",
                      "detrend by differences, which is the appropriate ",
                      "choice for temperature and for any variable that can ",
                      "be negative.",
                      if (!is.null(ratio_hint)) paste0(" ", ratio_hint)),
               call. = FALSE)
        }
      }
    }
  }

  # Settings are reported back with the results.
  if (is.null(dc_method)) {

    detrending_summary <- list(dc_method = "none")

  } else if (dc_method == "SLD") {

    detrending_summary <- list(dc_method = dc_method,
                               dc_scale = dc_scale)

  } else {

    detrending_summary <- list(dc_method = dc_method,
                               dc_nyrs = dc_nyrs,
                               dc_f = dc_f,
                               dc_difference = dc_difference,
                               dc_span = dc_span,
                               dc_bass = dc_bass,
                               dc_scale = dc_scale)
  }

  if (is.null(dc_method)) {
    dc_scale <- NULL
  }

  list(dc_method = dc_method,
       dc_scale = dc_scale,
       summary = detrending_summary,
       args = list(dc_method = dc_method,
                   dc_nyrs = dc_nyrs,
                   dc_f = dc_f,
                   dc_difference = dc_difference,
                   dc_span = dc_span,
                   dc_bass = dc_bass,
                   dc_scale = dc_scale))
}


#' Detrending settings for the control windows of the seascorr functions
#'
#' Internal helper. Starts from the settings given for the primary windows and
#' replaces those listed in dc_env_data_control.
#'
#' @param shared a named list of the dc_* arguments as the user supplied them
#' @param overrides NULL or the named list given as dc_env_data_control
#'
#' @return a named list of dc_* arguments for the control windows
#' @keywords internal
#' @noRd

dc_control_settings <- function(shared, overrides) {

  if (is.null(overrides)) {
    return(shared)
  }

  if (!is.list(overrides) || is.data.frame(overrides) ||
      (length(overrides) > 0 &&
       (is.null(names(overrides)) || any(!nzchar(names(overrides)))))) {
    stop(paste0("dc_env_data_control must be NULL or a named list of ",
                "detrending settings, for example ",
                "list(dc_difference = FALSE)."), call. = FALSE)
  }

  unknown <- setdiff(names(overrides), names(shared))

  if (length(unknown) > 0) {
    stop(paste0("dc_env_data_control contains unknown settings: ",
                paste(unknown, collapse = ", "), ". Allowed are ",
                paste(names(shared), collapse = ", "), "."), call. = FALSE)
  }

  if (anyDuplicated(names(overrides))) {
    stop("dc_env_data_control names a setting more than once.",
         call. = FALSE)
  }

  # single-bracket assignment keeps settings that are set to NULL
  for (nm in names(overrides)) {
    shared[nm] <- list(overrides[[nm]])
  }

  shared
}


#' Detrending function for one climate variable
#'
#' Internal helper. Returns a function of the aggregated series of one window
#' and its years. Without detrending the series is returned untouched, so that
#' results are exactly those of earlier versions.
#'
#' @param dc_settings a result of dc_check_settings()
#'
#' @return a function(x, years = NULL)
#' @keywords internal
#' @noRd

dc_make_apply <- function(dc_settings) {

  args <- dc_settings$args

  if (is.null(args$dc_method)) {
    return(function(x, years = NULL) x)
  }

  function(x, years = NULL) {

    detrend_climate(x,
                    dc_method = args$dc_method,
                    dc_nyrs = args$dc_nyrs,
                    dc_f = args$dc_f,
                    dc_difference = args$dc_difference,
                    dc_span = args$dc_span,
                    dc_bass = args$dc_bass,
                    dc_scale = isTRUE(args$dc_scale),
                    years = years)
  }
}
