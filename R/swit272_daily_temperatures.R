#' Daily temperatures for swit272 chronology
#'
#' Mean daily temperature in Celsius for the period 1950 - 2019.
#' This gridded E-OBS data on 0.1° regular grid, version 20e.
#' Extracted data is for the grid point with lon = 9.75 and lat = 46.45.
#'
#' @format A data frame with 25414 rows and 2 variables:
#' \describe{
#'   \item{date}{date as a factor, Year-Month-Day}
#'   \item{t_avg}{mean daily temperature in degrees Celsius}
#' }
#' @source \url{https://www.ecad.eu/download/ensembles/download.php}
#'
#' @details The E-OBS data are subject to the ECA&D data policy: they may only
#' be used for non-commercial research and non-commercial education.
#'
#' We acknowledge the E-OBS dataset from the EU-FP6 project UERRA (http://www.uerra.eu) and the Copernicus Climate Change Service, and the data providers in the ECA&D project (https://www.ecad.eu). Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: An Ensemble Version of the E-OBS Temperature and Precipitation Datasets, J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200
#'
#' @export
"swit272_daily_temperatures"
