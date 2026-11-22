#' Years in which a climate window can be completely empty
#'
#' Internal helper. The window aggregation removes missing values, so a window
#' is missing for a given year only when every day or month it covers is
#' missing. A year is therefore reported only when it contains a run of
#' missing values at least as long as the shortest window that is evaluated.
#' Scattered gaps, such as DOY 366 in non-leap years, are not reported.
#'
#' @param env a data frame or matrix of climate data, years as rows
#' @param min_window_width the shortest window width that is evaluated
#'
#' @return a character vector of years (row names)
#' @keywords internal
#' @noRd

window_gap_years <- function(env, min_window_width) {

  if (!anyNA(env)) {
    return(character(0))
  }

  longest_na_run <- function(z) {

    if (!anyNA(z)) {
      return(0L)
    }

    runs <- rle(is.na(z))
    max(runs$lengths[runs$values])
  }

  na_runs <- apply(as.matrix(env), 1, longest_na_run)

  row.names(env)[na_runs >= min_window_width]
}


#' Compact listing of years for warnings
#'
#' Internal helper. Lists at most ten years and counts the rest.
#'
#' @param years a character or numeric vector of years
#'
#' @return a character string
#' @keywords internal
#' @noRd

format_year_list <- function(years) {

  if (length(years) > 10) {
    paste0(paste(years[1:10], collapse = ", "), " and ",
           length(years) - 10, " more")
  } else {
    paste(years, collapse = ", ")
  }
}


#' Warning text for climate windows that are empty in some years
#'
#' Internal helper. Builds the warning issued when climate detrending meets
#' windows that are completely missing in some years.
#'
#' @param env_name name of the climate data frame, used in the message
#' @param gap_years years returned by window_gap_years()
#' @param min_window_width the shortest window width that is evaluated
#' @param unit "days" or "months"
#'
#' @return a character string
#' @keywords internal
#' @noRd

dc_gap_message <- function(env_name, gap_years, min_window_width,
                           unit = "days") {

  paste0(env_name, " contains gaps of at least ", min_window_width,
         " consecutive ", unit, ", the shortest window width, in the ",
         "following years: ", format_year_list(gap_years), ". Climate ",
         "detrending is fitted on the non-missing years of the affected ",
         "windows and the missing values are reinserted afterwards.")
}


#' Warning text for windows that return NA because of missing values
#'
#' Internal helper. Replaces the earlier anyNA() test, which also fired for
#' scattered missing values that the window aggregation removes. Only gaps
#' long enough to empty a whole window make a correlation NA when missing
#' values are not excluded.
#'
#' @param env_name name of the climate data frame, used in the message
#' @param gap_years years returned by window_gap_years()
#' @param min_window_width the shortest window width that is evaluated
#' @param unit "days" or "months"
#' @param na_arg name of the argument that controls missing values, e.g.
#' "cor_na_use"
#' @param na_value its current value
#'
#' @return a character string
#' @keywords internal
#' @noRd

na_gap_message <- function(env_name, gap_years, min_window_width,
                           unit = "days", na_arg = "cor_na_use",
                           na_value = "everything") {

  paste0(env_name, " contains gaps of at least ", min_window_width,
         " consecutive ", unit, ", the shortest window width, in the ",
         "following years: ", format_year_list(gap_years), ". Windows that ",
         "lie within such a gap are missing in those years, and with ",
         na_arg, " = '", na_value, "' they return NA. Set ", na_arg,
         " = 'pairwise.complete.obs' to compute over the available years ",
         "instead.")
}
