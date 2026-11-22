#' monthly_response
#'
#' Function calculates all possible values of a selected statistical metric
#' between one or more response variables and monthly sequences of environmental
#' data. Calculations are based on moving window which slides through monthly
#' environmental data. All calculated metrics are stored in a matrix. The
#' location of stored calculated metric in the matrix is indicating a window
#' width (row names) and a location in a matrix of monthly sequences of
#' environmental data (column names).
#'
#' @param response a data frame with tree-ring proxy variables as columns and
#' (optional) years as row names. Row.names should be matched with those from a
#' env_data data frame. If not, set row_names_subset = TRUE.
#' @param env_data a data frame of monthly sequences of environmental
#' data as columns and years as row names. Each row represents a year and
#' each column represents a day of a year (or month). Row.names should be matched with
#' those from a response data frame. If not, set row_names_subset = TRUE.
#' Alternatively, env_data could be a tidy data with three columns,
#' i.e. Year, DOY (Month) and third column representing values of mean temperatures,
#' sum of precipitation etc. If tidy data is passed to the function, set the argument
#' tidy_env_data to TRUE. For aggregate_function = 'SPEI', env_data holds the
#' monthly climatic water balance (precipitation minus potential
#' evapotranspiration), for example from water_balance().
#' @param lower_limit lower limit of window width (i.e. number of consecutive months
#' to be used for calculations)
#' @param upper_limit upper limit of window width (i.e. number of consecutive months
#' to be used for calculations)
#' @param fixed_width fixed width used for calculations (i.e. number of consecutive
#' months to be used for calculations)
#' @param method a character string specifying which method to use. Current
#' possibilities are "cor" (default), "lm" and "brnn".
#' @param metric a character string specifying which metric to use. Current
#' possibilities are "r.squared" and "adj.r.squared". If method = "cor",
#' metric is not relevant.
#' @param cor_method a character string indicating which correlation
#' coefficient is to be computed. One of "pearson" (default), "kendall", or
#' "spearman".
#' @param previous_year logical. If TRUE, previous-year climate data are included
#' in the analysis. If FALSE, no previous-year climate data are included and
#' number_previous_years is ignored.
#' @param number_previous_years integer between 1 and 5 specifying how many
#' previous years should be included in the environmental matrix when
#' previous_year = TRUE. For example, number_previous_years = 2 uses monthly
#' climate data from years t - 2, t - 1 and t for response year t. If NULL and
#' previous_year = TRUE, one previous year is included.
#' @param reference_window character string, the reference_window argument describes,
#' how each calculation is referred. There are two different options: 'start'
#' (default) and 'end'. If the reference_window argument is set to 'start',
#' then each calculation is related to the starting month of window. If the
#' reference_window argument is set to 'end', then each calculation is related
#' to the ending day of window calculation.
#' @param neurons positive integer that indicates the number of neurons used
#'  for brnn method
#' @param brnn_smooth if set to TRUE, a smoothing algorithm is applied that
#' removes unrealistic calculations which are a result of neural net failure.
#' @param remove_insignificant if set to TRUE, removes all correlations bellow
#' the significant threshold level, based on a selected alpha. For "lm" and
#' "brnn" method, squared threshold is used, which corresponds to R squared
#' statistics.
#' @param remove_insignificant_method the criterion used when
#' remove_insignificant = TRUE. One of "alpha" (per-window significance
#' threshold; the default), "boot" (bootstrap confidence interval excludes
#' zero; requires boot = TRUE), "fdr_BY" (false discovery rate,
#' Benjamini-Yekutieli, valid under arbitrary dependence between windows),
#' "fdr_BH" (Benjamini-Hochberg, valid under positive dependence and less
#' conservative), "maxstat" (max-statistic permutation test), "sidak"
#' (Sidak correction on the effective number of tests) or "cluster"
#' (cluster-extent inference). The "alpha" and "boot" criteria are applied per
#' window and therefore do not account for the number of windows screened; the
#' remaining criteria do.
#'
#' "maxstat" builds a null distribution of the largest absolute correlation
#' obtained anywhere in the search by repeatedly shuffling the response, and
#' retains only windows exceeding the resulting threshold. It makes no
#' distributional assumptions but requires perm_n re-evaluations of the whole
#' search. It is available for method = "cor" and for method = "lm" with a
#' single predictor, but not for method = "brnn".
#'
#' "sidak" targets the same family-wise error rate as "maxstat" but estimates
#' it analytically. The effective number of independent tests among the
#' screened windows is obtained from the eigenvalues of their correlation
#' matrix following Li and Ji (2005), and a Sidak correction is applied to
#' that number rather than to the nominal count of windows. No permutations
#' are required, so it is effectively instantaneous and is available for all
#' methods including "brnn".
#' Years that are missing in some of the windows, such as the first years
#' removed by dc_method = "Ar", are left out of this estimate, rather than
#' the windows themselves.
#'
#' "sidak" is an approximation to "maxstat". Across 54 combinations of site
#' and climate variable the two thresholds agreed to within 0.016 on average,
#' with a maximum discrepancy of 0.038. The approximation is close to exact
#' for strongly autocorrelated variables such as temperature, and tends to be
#' slightly permissive for weakly autocorrelated variables such as daily
#' precipitation, because the eigenvalue approach does not represent temporal
#' autocorrelation. "maxstat" is therefore preferable for final results,
#' while "sidak" is well suited to exploratory work and to method = "brnn".
#'
#' "cluster" is cluster-extent inference, adapted from neuroimaging. The
#' surface is thresholded at cluster_threshold, connected regions of
#' supra-threshold windows are identified, and each region's mass, the sum of
#' absolute correlations inside it, is compared against the largest mass
#' arising in permuted surfaces. It controls the family-wise error rate over
#' CLUSTERS rather than over windows, so it establishes that a region contains
#' signal without identifying which windows within it carry that signal.
#'
#' In tree-ring response surfaces this method has little discriminating power.
#' A real climate signal forms one large connected region whose mass greatly
#' exceeds anything the permuted surfaces produce, so every cluster that forms
#' is retained and the permutation step rarely rejects. The number of windows
#' kept is then determined by cluster_threshold rather than by the data. It is
#' provided for completeness and comparison; "maxstat" or "sidak" are
#' preferable when the question is which windows carry signal. It is available
#' only for method = "cor".
#' @param cluster_threshold the absolute correlation at which clusters are
#' formed when remove_insignificant_method = "cluster". Either "nominal", the
#' default, which uses the critical correlation at the chosen alpha, or a
#' single number between 0 and 1. This is a free parameter and it largely
#' determines how many windows are retained, so the value used is reported in
#' $significance_filter. Ignored by all other criteria.
#' @param perm_n number of permutations used when
#' remove_insignificant_method = "maxstat". Defaults to 999. Ignored by all
#' other criteria.
#' @param alpha significance level used to remove insignificant calculations.
#' @param row_names_subset if set to TRUE, row.names are used to subset
#' env_data and response data frames. Only years from both data frames are
#' kept.
#' @param aggregate_function character string specifying how the monthly data should be
#' aggregated. The default is 'mean', the other options are 'median', 'sum',
#' 'quantile' and 'SPEI'.
#' Missing months are skipped, so 'sum' adds only the months that are
#' available; a window without any value is NA.
#' 'SPEI' calculates the Standardized Precipitation-Evapotranspiration Index
#' (Vicente-Serrano et al. 2010) of each window. env_data must then hold the
#' monthly climatic water balance, i.e. precipitation minus potential
#' evapotranspiration (P - PET), which can be calculated with
#' water_balance(). The water balance is summed over the window and
#' standardised with a log-logistic distribution fitted by unbiased
#' probability weighted moments, the default of SPEI::spei(). Each window is
#' standardised against all years of env_data, before they are matched with
#' the response. A window of K months ending in month m therefore has the
#' value of SPEI::spei(scale = K) in month m, the SPEI of a year does not
#' depend on the analysed years, and a window in a previous year is
#' standardised as the same window in the current year. Unlike the other
#' options, 'SPEI' does not skip missing months: a window that includes a
#' missing month is NA in that year. SPEI values can only be detrended by
#' differences (dc_difference = TRUE).
#' @param quantile_prob numeric value between 0 and 1 specifying the quantile
#' probability used when aggregate_function = 'quantile'. For example,
#' quantile_prob = 0.95 calculates the 95th percentile. The default is 0.5.
#' @param temporal_stability_check character string, specifying, how temporal stability
#' between the optimal selection and response variable(s) will be analysed. Current
#' possibilities are "sequential", "progressive" and "running_window". Sequential check
#' will split data into k splits and calculate selected metric for each split. Progressive
#' check will split data into k splits, calculate metric for the first split and then
#' progressively add 1 split at a time and calculate selected metric. For running window,
#' select the length of running window with the k_running_window argument.
#' @param k_running_window the length of running window for temporal stability check.
#' Applicable only if temporal_stability argument is set to running window.
#' @param k integer, number of breaks (splits) for temporal stability and cross validation
#' analysis.
#' @param cross_validation_type character string, specifying, how to perform cross validation
#' between the optimal selection and response variables. If the argument is set to "blocked",
#' years will not be shuffled. If the argument is set to "randomized", years will be shuffled.
#' @param subset_years a subset of years to be analyzed. Should be given in the form of
#' subset_years = c(1980, 2005)
#' @param ylimits limit of the y axes for plot_extreme. It should be given in
#' the form of: ylimits = c(0,1)
#' @param seed optional seed argument for reproducible results
#' @param tidy_env_data if set to TRUE, env_data should be inserted as a data frame with three
#' columns: "Year", "Month", "Precipitation/Temperature/etc."
#' @param boot logical, if TRUE, bootstrap confidence intervals are calculated
#' for the correlation coefficients, R squared or adjusted R squared metrics and
#' returned as $boot_lower and $boot_upper. The metrics themselves are those of
#' the observed data, as with boot = FALSE. Versions before 2.0.0 reported the
#' mean of the bootstrap replicates instead.
#' @param boot_n The number of bootstrap replicates
#' @param boot_ci_type A character string representing the type of bootstrap
#' intervals: "norm" (default), "basic", "perc" or "bca". "bca" needs boot_n of
#' at least the number of analysed years. With fewer replicates
#' boot::boot.ci() cannot estimate bca intervals, so they are all NA, and a
#' warning is given. Studentized intervals ("stud") are not available.
#' @param boot_conf_int A scalar or vector containing the confidence level(s) of
#' the required interval(s)
#' @param month_interval a vector of two values defining the interval of months
#' used for calculations. Positive values indicate months in the current year.
#' Negative values indicate months in the previous-year block and are only used
#' when previous_year = TRUE. If previous_year = FALSE and negative values are
#' supplied, month_interval is ignored and the analysis is performed for the
#' current year only using month_interval = c(1, 12). If number_previous_years > 1,
#' the previous-year block starts with the earliest included previous year. For
#' example, previous_year = TRUE, number_previous_years = 2 and
#' month_interval = c(-1, 12) analyses the full sequence from January of year
#' t - 2 to December of year t.
#' @param dc_method a character string specifying how the aggregated climate
#' series of each moving window should be detrended before the climate-growth
#' relationship is calculated. Detrending is applied to the aggregated series,
#' i.e. to a vector with one value per year, and not to the raw monthly data.
#' Possible values are "none" (or NULL, the default, no detrending), "SLD"
#' (Simple Linear Detrending, the historical dendroTools option), and the
#' smoothers "Spline", "Mean", "Ar" and "Friedman", which follow
#' dplR::detrend(). The smoothers are fitted in chronological order of the
#' years given by the row names. "Ar" (prewhitening) leaves missing values at
#' the start of each detrended window, so for method = "cor" cor_na_use is
#' changed to "pairwise.complete.obs", with a warning, unless it already
#' excludes missing values.
#' The growth-specific options of dplR::detrend() ("ModNegExp",
#' "ModHugershoff", "AgeDepSpline") are not available, because they assume
#' strictly positive series with a biological growth trend.
#' @param dc_nyrs a number giving the rigidity of the smoothing spline when
#' dc_method = "Spline". Defaults to NULL, i.e. 0.67 times the length of the
#' series, exactly as in dplR::detrend(). Ignored by the other methods.
#' @param dc_f a number between 0 and 1 giving the frequency response at the
#' wavelength given by dc_nyrs when dc_method = "Spline". Defaults to 0.5.
#' Ignored by the other methods.
#' @param dc_difference logical. If FALSE (default, as in dplR::detrend()),
#' the detrended series is calculated as a ratio between the observed series
#' and the fitted curve. If TRUE, it is calculated as a difference. Ratio
#' detrending is only defined for strictly positive data, such as precipitation
#' sums, so the function stops if env_data contains negative values. Windows
#' whose aggregated series still contains a zero, for example a completely dry
#' month, are returned as NA. For temperature and for any other variable that
#' can be negative, use dc_difference = TRUE. Ignored when dc_method = "SLD",
#' which is always residual based.
#' @param dc_span the span used by the Friedman super smoother when
#' dc_method = "Friedman". Either "cv" (default, cross-validation) or a
#' number between 0 and 1. See also stats::supsmu(). Ignored by the other
#' methods.
#' @param dc_bass a number between 0 and 10 controlling the smoothness of the
#' Friedman super smoother when dc_method = "Friedman". Defaults to 0. See also
#' stats::supsmu(). Ignored by the other methods.
#' @param dc_scale logical or NULL. If TRUE, the detrended series is divided by
#' the standard deviation of its own values. The default (NULL) uses the
#' convention of the selected method: TRUE for dc_method = "SLD", which
#' preserves the behaviour of earlier dendroTools versions, and FALSE for the
#' dplR methods, which are not rescaled in dplR::detrend(). Rescaling does not
#' affect correlations, but it does affect the scale of $optimized_return and
#' of the transfer function.
#' @param cor_na_use an optional character string giving a method for computing
#' covariances in the presence of missing values for correlation coefficients.
#' This must be (an abbreviation of) one of the strings "everything" (default),
#' "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". See
#' also the documentation for the base cor() function.
#'
#' @return a list with 20 elements:
#' \enumerate{
#'  \item $calculations - a matrix with calculated metrics
#'  \item $method - the character string of a method
#'  \item $metric - the character string indicating the metric used for calculations
#'  \item $analysed_period - the character string specifying the analysed period based on the information from row names. If there are no row names, this argument is given as NA
#'  \item $optimized_return - data frame with two columns, response variable and aggregated (averaged) monthly data that return the optimal results. This data.frame could be directly used to calibrate a model for climate reconstruction
#'  \item $optimized_return_all - a data frame with aggregated monthly data, that returned the optimal result for the entire env_data (and not only subset of analysed years)
#'  \item $transfer_function - a ggplot object: scatter plot of optimized return and a transfer line of the selected method
#'  \item $temporal_stability - a data frame with calculations of selected metric for different temporal subsets
#'  \item $cross_validation - a data frame with cross validation results
#'  \item $plot_heatmap - ggplot2 object: a heatmap of calculated metrics
#'  \item $plot_extreme - ggplot2 object: line or bar plot of a row with the highest value in a matrix of calculated metrics
#'  \item $type - the character string describing type of analysis: daily or monthly
#'  \item $reference_window - character string, which reference window was used for calculations
#'  \item $boot_lower - matrix with lower limit of confidence intervals of bootstrap calculations
#'  \item $boot_upper - matrix with upper limit of confidence intervals of bootstrap calculations
#'  \item $significance_filter - list describing which significance filter was
#'  applied, the criterion used, and how many matrix cells were retained or removed
#'  \item $aggregated_climate - matrix with all aggregated climate series
#'  \item $previous_year - logical indicating whether previous-year climate data were used
#'  \item $number_previous_years - integer indicating how many previous years were used
#'  \item $detrending - information about the climate detrending applied,
#'  i.e. the method and the settings that were used, or "none"
#'}
#'
#' @references
#' Vicente-Serrano, S.M., Begueria, S., Lopez-Moreno, J.I. (2010) A
#' multiscalar drought index sensitive to global warming: the Standardized
#' Precipitation Evapotranspiration Index. Journal of Climate, 23, 1696-1718.
#' \doi{10.1175/2009JCLI2909.1}
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' # The examples below are enclosed within donttest{} to minimize the execution
#' # time during R package checks.
#'
#' # Load the dendroTools R package
#' library(dendroTools)
#'
#' # Load data used for examples
#' data(data_MVA)
#' data(data_TRW)
#' data(data_TRW_1)
#' data(example_proxies_individual)
#' data(example_proxies_1)
#' data(LJ_monthly_temperatures)
#' data(LJ_monthly_precipitation)
#'
#' # 1 Example with tidy precipitation data
#' example_tidy_data <- monthly_response(response = data_MVA,
#'     lower_limit = 1, upper = 24, dc_method = "SLD",
#'     env_data = LJ_monthly_precipitation, fixed_width = 0,
#'     method = "cor", row_names_subset = TRUE,
#'     remove_insignificant = FALSE, previous_year = FALSE,
#'     reference_window = "end",
#'     alpha = 0.05, aggregate_function = 'sum', boot = FALSE,
#'     tidy_env_data = TRUE, boot_n = 100, month_interval = c(-5, 10))
#'
#' # summary(example_tidy_data)
#' # plot(example_tidy_data, type = 1)
#' # plot(example_tidy_data, type = 2)
#'
#' # 2 Example with split data for early and late
#' example_MVA_early <- monthly_response(response = data_MVA,
#'     env_data = LJ_monthly_temperatures,
#'     method = "cor", row_names_subset = TRUE, previous_year = TRUE,
#'     remove_insignificant = TRUE, alpha = 0.05,
#'     subset_years = c(1940, 1980), aggregate_function = 'mean')
#'
#' example_MVA_late <- monthly_response(response = data_MVA,
#'     env_data = LJ_monthly_temperatures,
#'     method = "cor", row_names_subset = TRUE, alpha = 0.05,
#'     previous_year = TRUE, remove_insignificant = TRUE,
#'     subset_years = c(1981, 2010), aggregate_function = 'mean')
#'
#' # summary(example_MVA_late)
#' # plot(example_MVA_early, type = 1)
#' # plot(example_MVA_late, type = 1)
#' # plot(example_MVA_early, type = 2)
#' # plot(example_MVA_late, type = 2)
#'
#' # 3 Example negative correlations
#' example_neg_cor <- monthly_response(response = data_TRW_1, alpha = 0.05,
#'    env_data = LJ_monthly_temperatures,
#'    method = "cor", row_names_subset = TRUE,
#'    remove_insignificant = TRUE, boot = FALSE)
#'
#' # summary(example_neg_cor)
#' # plot(example_neg_cor, type = 1)
#' # plot(example_neg_cor, type = 2)
#' # example_neg_cor$temporal_stability
#'
#' # 4 Example of multiproxy analysis
#' # summary(example_proxies_1)
#' # cor(example_proxies_1)
#'
#' example_multiproxy <- monthly_response(response = example_proxies_1,
#'    env_data = LJ_monthly_temperatures,
#'    method = "lm", metric = "adj.r.squared",
#'    row_names_subset = TRUE, previous_year = FALSE,
#'    remove_insignificant = TRUE, alpha = 0.05)
#'
#' # summary(example_multiproxy)
#' # plot(example_multiproxy, type = 1)
#'
#' # 5 Example to test the temporal stability
#' example_MVA_ts <- monthly_response(response = data_MVA,
#'    env_data = LJ_monthly_temperatures,
#'    method = "lm", metric = "adj.r.squared", row_names_subset = TRUE,
#'    remove_insignificant = TRUE, alpha = 0.05,
#'    temporal_stability_check = "running_window", k_running_window = 10)
#'
#' # summary(example_MVA_ts)
#' # example_MVA_ts$temporal_stability
#'
#' # 6 Example using quantiles for the aggregation
#' example_q95 <- monthly_response(
#'   response = data_MVA,
#'   env_data = LJ_monthly_temperatures,
#'   method = "cor",
#'   fixed_width = 3,
#'   row_names_subset = TRUE,
#'   previous_year = TRUE,
#'   aggregate_function = "quantile",
#'   quantile_prob = 0.95
#' )
#'
#' # summary(example_q95)
#' # example_q95$temporal_stability
#'
#' # Example using two previous years plus the current year
#' example_two_previous_years <- monthly_response(
#'   response = data_MVA,
#'   env_data = LJ_monthly_temperatures,
#'   method = "cor",
#'   fixed_width = 3,
#'   row_names_subset = TRUE,
#'   previous_year = TRUE,
#'   number_previous_years = 2,
#'   month_interval = c(-1, 10),
#'   aggregate_function = "mean"
#' )
#'
#' # Example with climate detrending, following dplR::detrend().
#' # LJ_monthly_temperatures are temperatures in degrees Celsius and therefore
#' # cross zero, so differences must be used instead of the default ratios.
#' example_detrended <- monthly_response(response = data_MVA,
#'   env_data = LJ_monthly_temperatures,
#'   method = "cor", row_names_subset = TRUE,
#'   remove_insignificant = FALSE,
#'   dc_method = "Spline", dc_difference = TRUE)
#'
#' # example_detrended$detrending
#'
#' # Precipitation sums are positive, so the default ratios can be used
#' example_detrended_ratio <- monthly_response(response = data_MVA,
#'   env_data = LJ_monthly_precipitation, tidy_env_data = TRUE,
#'   method = "cor", row_names_subset = TRUE,
#'   aggregate_function = "sum", remove_insignificant = FALSE,
#'   dc_method = "Friedman")
#'
#' # Example with SPEI for the black pine site Sezana. The residual
#' # chronology is built as in Izmir et al. (2024).
#' SEZ_chron <- dplR::chron(dplR::detrend(SEZ_rwl, method = "Spline",
#'                                        nyrs = 32, f = 0.5),
#'                          prewhiten = TRUE)
#' SEZ_res <- SEZ_chron[!is.na(SEZ_chron$res), "res", drop = FALSE]
#'
#' monthly_wb <- water_balance(
#'   prec = data_transform(SEZ_daily_psum, format = "monthly",
#'                         monthly_aggregate_function = "sum"),
#'   tmin = data_transform(SEZ_daily_tmin, format = "monthly",
#'                         monthly_aggregate_function = "mean"),
#'   tmax = data_transform(SEZ_daily_tmax, format = "monthly",
#'                         monthly_aggregate_function = "mean"),
#'   lat = 45.71)
#'
#' example_spei <- monthly_response(response = SEZ_res,
#'   env_data = monthly_wb, method = "cor", row_names_subset = TRUE,
#'   previous_year = TRUE, aggregate_function = "SPEI",
#'   remove_insignificant = FALSE)
#'
#' # summary(example_spei)
#' # plot(example_spei, type = 2)
#'
#' }

monthly_response <- function(response, env_data, method = "cor",
                             metric = "r.squared", cor_method = "pearson",
                             previous_year = FALSE, number_previous_years = NULL, neurons = 1,
                             lower_limit = 1, upper_limit = 12,
                             fixed_width = 0, brnn_smooth = TRUE,
                             remove_insignificant = TRUE,
                             remove_insignificant_method = "alpha",
                             perm_n = 999,
                             cluster_threshold = "nominal",
                             alpha = .05, row_names_subset = FALSE,
                             reference_window = "start",
                             aggregate_function = 'mean',
                             quantile_prob = 0.5,
                             temporal_stability_check = "sequential", k = 2,
                             k_running_window = 30, cross_validation_type = "blocked",
                             subset_years = NULL,
                             ylimits = NULL, seed = NULL, tidy_env_data = FALSE,
                             boot = FALSE, boot_n = 1000, boot_ci_type = "norm",
                             boot_conf_int = 0.95,
                             month_interval = NULL,
                             dc_method = NULL,
                             dc_nyrs = NULL,
                             dc_f = 0.5,
                             dc_difference = FALSE,
                             dc_span = "cv",
                             dc_bass = 0,
                             dc_scale = NULL,
                             cor_na_use = "everything"
) {

  ##############################################################################
  # 1 month interval is organized

  months_per_year <- 12L

  if (reference_window == "middle"){

    stop(paste0("reference_window should be 'start' or 'end'. 'middle' reference_window is not implemented for the monthly_response"))
  }

  # previous_year is the master switch for previous-year analyses.
  # If previous_year = FALSE, number_previous_years is ignored.
  if (!is.logical(previous_year) ||
      length(previous_year) != 1 ||
      is.na(previous_year)) {
    stop("previous_year must be either TRUE or FALSE.")
  }

  if (isFALSE(previous_year)) {

    if (!is.null(number_previous_years) &&
        is.numeric(number_previous_years) &&
        length(number_previous_years) == 1 &&
        !is.na(number_previous_years) &&
        number_previous_years > 0) {
      warning(paste0("number_previous_years is ignored because ",
                     "previous_year = FALSE."))
    }

    number_previous_years <- 0L

  } else {

    # Backward compatibility:
    # previous_year = TRUE still means one previous year unless
    # number_previous_years is explicitly supplied.
    if (is.null(number_previous_years)) {
      number_previous_years <- 1L
    }

    if (!is.numeric(number_previous_years) ||
        length(number_previous_years) != 1 ||
        is.na(number_previous_years) ||
        number_previous_years %% 1 != 0 ||
        number_previous_years < 1 ||
        number_previous_years > 5) {
      stop(paste0("When previous_year = TRUE, number_previous_years must be ",
                  "a single integer between 1 and 5."))
    }

    number_previous_years <- as.integer(number_previous_years)
  }

  year_block_width <- months_per_year * (number_previous_years + 1L)

  # Default interval:
  # current year only: c(1, 12)
  # previous-year analysis: c(-1, 12), i.e. from January of the earliest
  # included previous year to December of the current year.
  if (is.null(month_interval)) {
    month_interval <- if (isTRUE(previous_year)) c(-1, months_per_year) else c(1, months_per_year)
  }

  if (!is.numeric(month_interval) ||
      length(month_interval) != 2 ||
      any(is.na(month_interval))) {
    stop("month_interval must be a numeric vector of length 2.")
  }

  offset_start <- month_interval[1]
  offset_end <- month_interval[2]

  if (offset_start == 0 || offset_end == 0) {
    stop("month_interval cannot contain 0. Use negative values for previous-year months and positive values for current-year months.")
  }

  # Negative values in month_interval would imply previous-year climate data.
  # If previous_year = FALSE, previous-year data are not used. In this case,
  # ignore month_interval and use the full current year.
  if (isFALSE(previous_year) && (offset_start < 0 || offset_end < 0)) {

    warning(paste0("Negative values were supplied in month_interval, but ",
                   "previous_year = FALSE. The month_interval argument is ignored ",
                   "and the analysis is performed for the current year only ",
                   "using month_interval = c(1, 12)."))

    month_interval <- c(1, months_per_year)
    offset_start <- month_interval[1]
    offset_end <- month_interval[2]
  }

  # If both limits are positive while previous_year = TRUE, previous years are
  # not actually included in the selected interval. In that case, keep the
  # analysis current-year only.
  if (offset_start > 0 && offset_end > 0 && isTRUE(previous_year)) {

    number_previous_years <- 0L
    previous_year <- FALSE
    year_block_width <- months_per_year

    warning(paste0("Previous years are not included in selected month_interval. ",
                   "The analysis is performed for the current year only."))
  }

  # Convert the user-facing month_interval to column positions in the lagged matrix.
  # For number_previous_years = 2, the internal matrix is:
  # columns 1:12     = year t - 2
  # columns 13:24    = year t - 1
  # columns 25:36    = year t
  if (offset_start < 0 && offset_end < 0) {

    offset_start <- abs(offset_start)
    offset_end <- abs(offset_end)

  } else if (offset_start < 0 && offset_end > 0) {

    offset_start <- abs(offset_start)
    offset_end <- offset_end + months_per_year * number_previous_years

  } else if (offset_start > 0 && offset_end > 0) {

    # Current-year-only interval; already in current-year coordinates.
    offset_start <- offset_start
    offset_end <- offset_end

  } else {

    stop("month_interval must not run from current-year months back to previous-year months.")
  }

  if (offset_start > offset_end) {
    stop("month_interval is invalid after conversion. The start of the interval is after the end.")
  }

  if (offset_start < 1 || offset_end > year_block_width) {
    stop(paste0("month_interval is outside the available climate sequence. ",
                "For number_previous_years = ", number_previous_years,
                ", the available width is ", year_block_width, " months."))
  }

  # Calculate the max_window allowed
  max_window <- offset_end - offset_start + 1

  # If max_window is smaller than upper_limit, upper_limit must be reduced
  if (upper_limit > max_window){

    upper_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The upper_limit is outside your month_interval and",
                     " therefore reduced to the maximum allowed: ", max_window, "."))
    }
  }

  if (lower_limit > max_window){
    lower_limit <- max_window

    if (fixed_width == 0){
      warning(paste0("The lower_limit is outside your month_interval and",
                     " therefore reduced to the minimum allowed: ", max_window, "."))
    }
  }

  # Also correction for fixed_window approach - can only be ERROR
  if (fixed_width > max_window){

    stop(paste0("The selected fixed_width is outside your month_interval. ",
                "Decrease the fixed_width argument to at least: ", max_window, "."))
  }

  # offset_end is converted to the number of months that must be skipped at the
  # end of the full environmental matrix.
  offset_end <- year_block_width - offset_end

  ################################################################################

  if (!is.null(seed)) {
    set.seed(seed)
  }


  # Studentized intervals need a variance estimate for every bootstrap
  # replicate. The statistics computed here return a single value per
  # replicate, so boot.ci() cannot build a studentized interval and previously
  # failed in a way that produced an uninformative error.

  allowed_boot_ci_types <- c("norm", "perc", "basic", "bca")

  if (!(boot_ci_type %in% allowed_boot_ci_types)) {
    stop(paste0("boot_ci_type is '", boot_ci_type, "'. Instead it should be ",
                "one of: ", paste(allowed_boot_ci_types, collapse = ", "),
                ". Studentized intervals are not available because they ",
                "require a variance estimate for each bootstrap replicate."))
  }
  allowed_significance_methods <- c("alpha", "boot", "fdr_BY", "fdr_BH",
                                    "maxstat", "sidak", "cluster")

  if (!(remove_insignificant_method %in% allowed_significance_methods)) {
    stop(paste0("remove_insignificant_method is '", remove_insignificant_method,
                "'. Instead it should be one of: ",
                paste(allowed_significance_methods, collapse = ", "), "."))
  }

  # The bootstrap filter needs bootstrap intervals to exist.
  if (isTRUE(remove_insignificant) &&
      remove_insignificant_method == "boot" && isFALSE(boot)) {

    warning(paste0("remove_insignificant_method = 'boot' requires boot = TRUE. ",
                   "No significance filter is applied."))

    remove_insignificant <- FALSE
  }

  # The max-statistic permutation test re-correlates a shuffled response with
  # every aggregated window. For "cor" this is a vectorised operation, and for
  # a single-predictor "lm" the R squared equals the squared correlation, so the
  # same null distribution applies. For "brnn", and for "lm" with more than one
  # predictor, no such shortcut exists and the test would require refitting the
  # model at every window in every permutation.

  # Cluster-extent inference needs the same permutation machinery as the
  # max-statistic test, plus a threshold at which clusters are formed.
  if (isTRUE(remove_insignificant) &&
      remove_insignificant_method == "cluster") {

    if (!is.numeric(perm_n) || length(perm_n) != 1 || is.na(perm_n) ||
        perm_n < 99) {
      stop("perm_n must be a single number of at least 99.")
    }

    if (!(identical(cluster_threshold, "nominal") ||
          (is.numeric(cluster_threshold) && length(cluster_threshold) == 1 &&
           !is.na(cluster_threshold) && cluster_threshold > 0 &&
           cluster_threshold < 1))) {
      stop(paste0("cluster_threshold must be \"nominal\" or a single number ",
                  "between 0 and 1 giving the absolute correlation at which ",
                  "clusters are formed."))
    }

    if (method != "cor") {
      stop(paste0("remove_insignificant_method = 'cluster' is available only ",
                  "for method = 'cor', because the permuted statistic has to ",
                  "be laid out on the same surface as the observed one. Use ",
                  "'sidak', 'fdr_BY', 'fdr_BH' or 'alpha' instead."))
    }
  }

  if (isTRUE(remove_insignificant) &&
      remove_insignificant_method == "maxstat") {

    if (!is.numeric(perm_n) || length(perm_n) != 1 || is.na(perm_n) ||
        perm_n < 99) {
      stop("perm_n must be a single number of at least 99.")
    }

    if (method == "brnn") {
      stop(paste0("remove_insignificant_method = 'maxstat' is not available for ",
                  "method = 'brnn', because the permutation would require ",
                  "refitting the network at every window. Use 'fdr_BY', ",
                  "'fdr_BH' or 'alpha' instead."))
    }

    if (method == "lm" && ncol(response) > 1) {
      stop(paste0("remove_insignificant_method = 'maxstat' is only available for ",
                  "method = 'lm' with a single predictor. Use 'fdr_BY', ",
                  "'fdr_BH' or 'alpha' instead."))
    }
  }

  if (fixed_width != 0){
    lower_limit = 1
    upper_limit = 12
  }


  # Defining global variables
  median <- NULL
  proxy <- NULL
  optimized_return <- NULL
  transfer_f <- NULL
  journal_theme <- NULL
  CV <- NULL
  Period <- NULL
  Years <- NULL
  yearABC <- NULL
  RMSE <- NULL
  RE <- NULL
  CE <- NULL
  DE <- NULL
  d <- NULL

  temporal_matrix_lower <- NULL
  temporal_matrix_upper <- NULL

  # Check selected aggregation function
  allowed_aggregate_functions <- c("mean", "median", "sum", "quantile", "SPEI")

  if (!(aggregate_function %in% allowed_aggregate_functions)) {
    stop(paste0(
      "aggregate_function is '", aggregate_function,
      "'. Instead it should be one of: ",
      paste(allowed_aggregate_functions, collapse = ", "), "."
    ))
  }

  # Check quantile probability
  if (!is.numeric(quantile_prob) ||
      length(quantile_prob) != 1 ||
      is.na(quantile_prob) ||
      quantile_prob < 0 ||
      quantile_prob > 1) {
    stop("quantile_prob must be a single numeric value between 0 and 1.")
  }

  # Internal helper for aggregating monthly data across rows. data is the
  # climate data frame and cols are the columns of the window.
  aggregate_monthly_window <- function(data, cols) {

    # SPEI is standardised against all years of env_data, which are kept
    # in spei_ref, so the values do not depend on the analysed years.
    if (aggregate_function == "SPEI") {
      return(spei_window_values(spei_ref, cols, years = row.names(data)))
    }

    x <- data.frame(data[, cols, drop = FALSE])

    if (aggregate_function == "mean") {
      return(rowMeans(x, na.rm = TRUE))
    }

    if (aggregate_function == "median") {
      return(apply(x, 1, median, na.rm = TRUE))
    }

    # Missing months are skipped, so a sum adds the values that are
    # available. A window without any value is NA: sum() would return 0,
    # i.e. a dry window that was never observed.
    if (aggregate_function == "sum") {
      out <- apply(x, 1, sum, na.rm = TRUE)
      out[rowSums(!is.na(x)) == 0] <- NA
      return(out)
    }

    if (aggregate_function == "quantile") {
      return(apply(
        x, 1, quantile,
        probs = quantile_prob,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      ))
    }
  }




  # Internal helper for robust BRNN model fitting.
  # BRNN can occasionally be numerically unstable for weak or noisy
  # relationships. In such cases, the calculation should return NA rather than
  # stopping the entire moving-window analysis.
  safe_brnn_metrics <- function(temporal_df,
                                formula,
                                observed,
                                neurons,
                                tol = 1e-6) {

    out <- list(
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      model = NULL,
      predictions = NULL
    )

    temporal_df <- data.frame(temporal_df)
    observed <- as.numeric(observed)

    if (nrow(temporal_df) != length(observed) ||
        length(observed) < 3 ||
        all(is.na(observed))) {
      return(out)
    }

    fit <- try({
      capture.output(
        tmp_model <- brnn(formula,
                          data = temporal_df,
                          neurons = neurons,
                          tol = tol)
      )
      tmp_model
    }, silent = TRUE)

    if (inherits(fit, "try-error") || is.null(fit)) {
      return(out)
    }

    predictions <- try(as.numeric(predict.brnn(fit, temporal_df)),
                       silent = TRUE)

    if (inherits(predictions, "try-error") ||
        length(predictions) != length(observed) ||
        any(!is.finite(predictions))) {
      return(out)
    }

    ss_total <- sum((observed - mean(observed, na.rm = TRUE)) ^ 2,
                    na.rm = TRUE)
    ss_res <- sum((observed - predictions) ^ 2, na.rm = TRUE)

    if (!is.finite(ss_total) ||
        !is.finite(ss_res) ||
        ss_total <= 0) {
      return(out)
    }

    r_squared <- 1 - ss_res / ss_total

    if (!is.finite(r_squared)) {
      r_squared <- NA_real_
    }

    n_obs <- length(observed)
    n_predictors <- ncol(temporal_df) - 1L
    denominator <- n_obs - n_predictors - 1L

    adj_r_squared <- if (is.finite(r_squared) && denominator > 0) {
      1 - ((1 - r_squared) * (n_obs - 1)) / denominator
    } else {
      NA_real_
    }

    if (!is.finite(adj_r_squared)) {
      adj_r_squared <- NA_real_
    }

    out$r.squared <- r_squared
    out$adj.r.squared <- adj_r_squared
    out$model <- fit
    out$predictions <- predictions

    out
  }

  safe_brnn_boot_metrics <- function(temporal_df,
                                     neurons,
                                     boot_n,
                                     boot_conf_int,
                                     boot_ci_type) {

    out <- list(
      r.squared = NA_real_,
      adj.r.squared = NA_real_,
      r.squared.lower = NA_real_,
      r.squared.upper = NA_real_,
      adj.r.squared.lower = NA_real_,
      adj.r.squared.upper = NA_real_
    )

    calc <- try(
      boot(data = temporal_df,
           statistic = boot_f_brnn,
           R = boot_n,
           brnn.formula = "x ~ .",
           neurons = neurons),
      silent = TRUE
    )

    if (inherits(calc, "try-error") ||
        is.null(calc$t) ||
        ncol(calc$t) < 2) {
      return(out)
    }

    r_squared <- as.numeric(calc$t0[1])
    adj_r_squared <- as.numeric(calc$t0[2])

    if (is.finite(r_squared)) {
      out$r.squared <- r_squared
    }

    if (is.finite(adj_r_squared)) {
      out$adj.r.squared <- adj_r_squared
    }

    ci_int_r_squared <- try(
      boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1),
      silent = TRUE
    )

    ci_int_adj_r_squared <- try(
      boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2),
      silent = TRUE
    )

    if (!inherits(ci_int_r_squared, "try-error")) {

      if (boot_ci_type == "norm") {
        out$r.squared.lower <- ci_int_r_squared$norm[2]
        out$r.squared.upper <- ci_int_r_squared$norm[3]
      } else if (boot_ci_type == "perc") {
        out$r.squared.lower <- ci_int_r_squared$perc[4]
        out$r.squared.upper <- ci_int_r_squared$perc[5]
      } else if (boot_ci_type == "stud") {
        out$r.squared.lower <- ci_int_r_squared$student[4]
        out$r.squared.upper <- ci_int_r_squared$student[5]
      } else if (boot_ci_type == "basic") {
        out$r.squared.lower <- ci_int_r_squared$basic[4]
        out$r.squared.upper <- ci_int_r_squared$basic[5]
      } else if (boot_ci_type == "bca") {
        out$r.squared.lower <- ci_int_r_squared$bca[4]
        out$r.squared.upper <- ci_int_r_squared$bca[5]
      } else {
        stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")
      }
    }

    if (!inherits(ci_int_adj_r_squared, "try-error")) {

      if (boot_ci_type == "norm") {
        out$adj.r.squared.lower <- ci_int_adj_r_squared$norm[2]
        out$adj.r.squared.upper <- ci_int_adj_r_squared$norm[3]
      } else if (boot_ci_type == "perc") {
        out$adj.r.squared.lower <- ci_int_adj_r_squared$perc[4]
        out$adj.r.squared.upper <- ci_int_adj_r_squared$perc[5]
      } else if (boot_ci_type == "stud") {
        out$adj.r.squared.lower <- ci_int_adj_r_squared$student[4]
        out$adj.r.squared.upper <- ci_int_adj_r_squared$student[5]
      } else if (boot_ci_type == "basic") {
        out$adj.r.squared.lower <- ci_int_adj_r_squared$basic[4]
        out$adj.r.squared.upper <- ci_int_adj_r_squared$basic[5]
      } else if (boot_ci_type == "bca") {
        out$adj.r.squared.lower <- ci_int_adj_r_squared$bca[4]
        out$adj.r.squared.upper <- ci_int_adj_r_squared$bca[5]
      } else {
        stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")
      }
    }

    # Clean non-finite CI values
    numeric_names <- names(out)
    for (tmp_name in numeric_names) {
      if (!is.finite(out[[tmp_name]])) {
        out[[tmp_name]] <- NA_real_
      }
    }

    out
  }


  # Internal helper for robust cross-validation output.
  # calculate_metrics() returns two rows per fold: calibration and validation.
  # If a BRNN fold fails, we must still return two rows of NA values; otherwise
  # cbind(Years = m2, m1) fails because m2 has 2 * k rows while m1 has fewer.
  empty_cv_metrics <- function(template = NULL) {

    if (!is.null(template)) {

      template <- as.data.frame(template)
      out <- template[rep(NA_integer_, 2), , drop = FALSE]
      row.names(out) <- NULL

      for (tmp_name in names(out)) {

        if (is.factor(out[[tmp_name]])) {
          out[[tmp_name]] <- factor(rep(NA_character_, 2),
                                    levels = levels(out[[tmp_name]]))
        } else if (is.character(out[[tmp_name]])) {
          out[[tmp_name]] <- rep(NA_character_, 2)
        } else {
          out[[tmp_name]] <- rep(NA_real_, 2)
        }
      }

      return(out)
    }

    data.frame(
      cor = c(NA_real_, NA_real_),
      RMSE = c(NA_real_, NA_real_),
      RRSE = c(NA_real_, NA_real_),
      d = c(NA_real_, NA_real_),
      RE = c(NA_real_, NA_real_),
      CE = c(NA_real_, NA_real_),
      DE = c(NA_real_, NA_real_),
      bias = c(NA_real_, NA_real_),
      check.names = FALSE
    )
  }

  normalize_cv_metric_list <- function(metric_list) {

    template <- NULL

    for (tmp_item in metric_list) {

      tmp_df <- try(as.data.frame(tmp_item), silent = TRUE)

      if (!inherits(tmp_df, "try-error") &&
          nrow(tmp_df) == 2 &&
          ncol(tmp_df) > 1) {
        template <- tmp_df
        break
      }
    }

    if (is.null(template)) {
      template <- empty_cv_metrics()
    }

    for (tmp_i in seq_along(metric_list)) {

      tmp_df <- try(as.data.frame(metric_list[[tmp_i]]), silent = TRUE)

      if (inherits(tmp_df, "try-error") ||
          nrow(tmp_df) != 2 ||
          ncol(tmp_df) != ncol(template) ||
          !identical(names(tmp_df), names(template))) {

        metric_list[[tmp_i]] <- empty_cv_metrics(template)

      } else {

        metric_list[[tmp_i]] <- tmp_df
      }
    }

    metric_list
  }



  # Internal helper for constructing a multi-year climate matrix.
  # Rows are current response years. Columns are ordered from the oldest previous
  # year to the current year.
  build_lagged_env_data <- function(x, n_previous_years) {

    x <- data.frame(x, check.names = FALSE)

    if (n_previous_years == 0L) {
      return(x)
    }

    env_years <- suppressWarnings(as.integer(row.names(x)))

    if (any(is.na(env_years))) {
      stop("For previous-year analyses, row.names of env_data must be calendar years.")
    }

    if (any(duplicated(env_years))) {
      stop("Duplicated years are present in row.names(env_data).")
    }

    row.names(x) <- as.character(env_years)
    x <- x[order(env_years), , drop = FALSE]
    env_years <- as.integer(row.names(x))

    current_years <- env_years[
      vapply(env_years, function(y) {
        all((y - n_previous_years):y %in% env_years)
      }, logical(1))
    ]

    if (length(current_years) == 0L) {
      stop(paste0("No years in env_data have the required ",
                  n_previous_years, " previous year(s)."))
    }

    lag_blocks <- lapply(seq(n_previous_years, 0L), function(lag_i) {

      block <- x[as.character(current_years - lag_i), , drop = FALSE]
      row.names(block) <- as.character(current_years)
      colnames(block) <- paste0("Y", -lag_i, "_", colnames(x))

      block
    })

    out <- do.call(cbind, lag_blocks)
    row.names(out) <- as.character(current_years)

    data.frame(out, check.names = FALSE)
  }

  # If there is a column name samp.depth in response data frame, warning is given
  if ("samp.depth" %in% colnames(response)){

    samp.depth_index <- grep("samp.depth", colnames(response))
    response <- response[, -samp.depth_index,F]

    warning("Removed the samp.depth from response data frame")
  }

  # If there is more than 2 columns in response data frame, give a warning
  if (ncol(response) > 1){
    warning(paste0("Your response data frame has more than 1 column! Are you doing a multiproxy research?",
                   " If so, OK. If not, check your response data frame!"))
  }


  # If you insert daily data into monthly response, you get an error
  if ((ncol(env_data) != 12) && (tidy_env_data == FALSE)){
    stop(paste0("You must insert env_data with 12 columns (months)! Instead, you have ", ncol(env_data),
                " columns!"))
  }

  # If env_data is given in tidy version, transformation is needed
  if (tidy_env_data == TRUE){

    n_col_tidy_DF <- ncol(env_data)
    colnames_tidy_DF <- colnames(env_data)

    if (ncol(env_data) != 3){
      stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                 "env_data should have 3 columns, but it has", n_col_tidy_DF, "instead!"))
    }

    if (colnames_tidy_DF[1] != "Year"){
      stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                 "The first column name of the env_data should be 'Year', but it is",
                 colnames_tidy_DF[1], "instead!"))
    }

    if (colnames_tidy_DF[2] != "Month"){
      stop(paste("env_data was inserted in tidy version (tidy_env_data is set to TRUE).",
                 "The second column name of the env_data should be 'DOY', but it is",
                 colnames_tidy_DF[2], "instead!"))
    }

    value_variable = colnames(env_data)[3]
    env_data <- dcast(env_data, Year~Month, value.var = value_variable)
    env_data <- years_to_rownames(env_data, "Year")



  }


  # PART 1 - general data arrangements, warnings and stops
  # Both bojects (response and env_data) are converted to data frames
  response <- data.frame(response)
  env_data <- data.frame(env_data)

  # Here we save the original response data that will be used later
  response_original <- response

  # Previous-year analyses require alignment by calendar year names. This avoids
  # accidental row-position mismatches after adding lagged climate years.
  if (number_previous_years > 0L && row_names_subset == FALSE) {
    row_names_subset <- TRUE
    warning(paste0("Previous-year analyses require alignment by year names. ",
                   "row_names_subset is set to TRUE."))
  }

  # SPEI is standardised against all years of env_data, before they are
  # lagged and matched with the response (see spei_reference()).
  env_data_current_years <- env_data

  # Build current-year or multi-year environmental matrix.
  # For number_previous_years = 0 this returns env_data unchanged.
  env_data <- build_lagged_env_data(env_data, number_previous_years)

  spei_ref <- NULL
  if (aggregate_function == "SPEI") {
    spei_ref <- spei_reference(env_data_current_years,
                               n_previous_years = number_previous_years,
                               daily = FALSE)
  }
  rm(env_data_current_years)

  # This object is used later for optimized_return_all. It should represent the
  # full lagged climate matrix before subset_years is applied.
  env_data_original <- env_data

  # For metric calculations, both objects need to have the same length,
  # with the exception, when row_names_subset is set to TRUE
  # Stop message in case both data frames do not have the same length
  if (nrow(response) != nrow(env_data) & row_names_subset == FALSE)
    stop("Length of env_data and response records differ")

  # Stop message if fixed_width is not between 0 and the full monthly sequence
  if (fixed_width < 0 | fixed_width > year_block_width)
    stop(paste0("fixed_width should be between 0 and ", year_block_width))

  # Stop in case of method == "cor" and ncol(proxies) > 1
  # Correlations could be calculated only for one variable
  if (method == "cor" & ncol(response) > 1)
    stop(paste("More than 1 variable in response data frame not suitable ",
               "for 'cor' method. Use 'lm' or 'brnn'"))

  # Rules for selected window limits
  if (lower_limit > upper_limit)
    stop("lower_limit can not be higher than upper_limit!")

  if (lower_limit > year_block_width | lower_limit < 1)
    stop(paste0("lower_limit out of bounds! It should be between 1 and ",
                year_block_width))

  if (upper_limit > year_block_width | upper_limit < 1)
    stop(paste0("upper_limit out of bounds! It should be between 1 and ",
                year_block_width))

  # ---- climate detrending (dc_* arguments) ---------------------------------
  # Detrending is applied to the aggregated climate series of each moving
  # window, i.e. to a vector with one value per year, exactly as in
  # daily_response(). dc_method = NULL or "none" disables it. "SLD" is the
  # historical dendroTools option. The remaining options follow
  # dplR::detrend().
  spei_check_detrending(aggregate_function, dc_method, dc_difference)

  dc_settings <- dc_check_settings(dc_method = dc_method,
                                   dc_nyrs = dc_nyrs,
                                   dc_f = dc_f,
                                   dc_difference = dc_difference,
                                   dc_span = dc_span,
                                   dc_bass = dc_bass,
                                   dc_scale = dc_scale,
                                   env_list = list(env_data = env_data))

  dc_method <- dc_settings$dc_method
  dc_scale <- dc_settings$dc_scale
  detrending_summary <- dc_settings$summary

  if (!is.null(dc_method)) {

    # Prewhitening returns one missing residual for each lag of the fitted AR
    # model, so every detrended window starts with NA values. With the default
    # cor_na_use = 'everything' this would return NA for the whole matrix.
    # 'pairwise.complete.obs' is used rather than 'complete.obs': for one
    # window it gives the same value, but it returns NA instead of an error
    # for a window without complete years (ratio detrending of a window
    # containing a zero), and the permutation filters correlate all windows
    # at once, where 'complete.obs' would drop a year from every window as
    # soon as it is missing in one of them.
    if (dc_method == "Ar" && method == "cor" &&
        !(cor_na_use %in% c("complete.obs", "na.or.complete",
                            "pairwise.complete.obs"))) {

      warning(paste0("dc_method = 'Ar' introduces missing values at the start ",
                     "of every detrended window. cor_na_use was changed from '",
                     cor_na_use, "' to 'pairwise.complete.obs' so that ",
                     "correlations can be calculated."))

      cor_na_use <- "pairwise.complete.obs"
    }
  }

  # Internal helper. Detrends the aggregated climate series of one window. With
  # dc_method = NULL the series is returned untouched, so that results without
  # detrending are exactly those of earlier versions. The years are used to
  # fit the dplR methods in chronological order.
  dc_apply <- dc_make_apply(dc_settings)


  # Warn users in case of missing values. The original threshold was 8 missing
  # months for one 12-month climate year, so it is scaled when multiple
  # previous years are included.
  missing_month_threshold <- 8 * (number_previous_years + 1L)
  env_temp <- env_data[row.names(env_data) %in% row.names(response),]

  # Subset of years
  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    subset_seq <- seq(lower_subset, upper_subset)
    env_temp <- subset(env_temp, row.names(env_temp) %in% subset_seq)
  }

  na_problem <- data.frame(na_sum = rowSums(is.na(env_temp)))
  na_problem <- na_problem[na_problem$na_sum > missing_month_threshold, , F]
  problematic_years <- paste0(row.names(na_problem), sep = "", collapse=", ")

  if (nrow(na_problem) > 0){

    warning(paste0("Problematic years with missing values are present: ", problematic_years))

  }

  # Previous-year climate data have already been added by
  # build_lagged_env_data().

  # If row_names_subset == TRUE, data is subseted and ordered based on matching
  # row.names. Additionally, number of characters in row.names is checked.
  # There should be at least three characters (assuming years before 100 will
  # never be analysed, there is no such environmental data available)
  if (row_names_subset == TRUE & nchar(row.names(env_data)[1]) >= 3){

    ncol_response <- ncol(response)

    colnames_response <- colnames(response)

    env_data$yearABC <- row.names(env_data)
    response$yearABC <- row.names(response)

    temporal_data <- merge(response, env_data, by = "yearABC")

    response <- data.frame(temporal_data[, c(2:(1 + ncol_response))],
                           row.names = temporal_data$yearABC)
    colnames(response) <- colnames_response

    env_data <- data.frame(temporal_data[, c((1 + ncol_response + 1):
                                               ncol(temporal_data))],
                           row.names = temporal_data$yearABC)
  }

  # if row.names of env_data and the response data frames are not equal,
  # warning is given.
  if (sum(row.names(env_data) == row.names(response)) != nrow(env_data)) {
    warning("row.names between env_data and response do not match!")
  }

  # If row_names_subset == TRUE, but row.names does not appear to be years,
  # error is given.
  if (row_names_subset == TRUE & nchar(row.names(env_data)[1]) < 3){
    stop(paste("row.names does not appear to be years!",
               "At least three characters needed!"))
  }

  # Subset of years
  if (!is.null(subset_years)){
    lower_subset <- subset_years[1]
    upper_subset <- subset_years[2]

    if (lower_subset > upper_subset){
      stop("Change the order of elements in the subset_years argument! First element should be lower than the second!")
    }

    subset_seq <- seq(lower_subset, upper_subset)

    if (any(!(subset_seq %in% row.names(response)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  " in the response data frame. Change the subset_years argument"))
    }

    if (any(!(subset_seq %in% row.names(env_data)))){

      stop(paste0("Undefined columns selected. Subset years don't exist",
                  " in the env_data data frame. Change the subset_years argument"))
    }

    response <- subset(response, row.names(response) %in% subset_seq)
    env_data <- subset(env_data, row.names(env_data) %in% subset_seq)
  }

  # ---- years with gaps long enough to empty a whole window -----------------
  # aggregate_monthly_window() skips missing values, so a window is missing
  # for a given year only when every month it covers is missing. Scattered
  # missing months are absorbed by the aggregation, so anyNA(env_data) is far
  # too broad a test. A year is only a problem when it contains a run of
  # missing months at least as long as the shortest window that is evaluated.
  # Only the analysed years are checked, so a gap in a year outside them, for
  # example a recent year without tree-ring data, raises no warning.
  min_window_width <- if (fixed_width != 0) fixed_width else lower_limit

  years_with_window_gaps <- window_gap_years(env_data, min_window_width)

  # SPEI does not skip missing values: a window that includes a missing month
  # has no SPEI in that year, so all years with missing months that enter the
  # analysis are reported instead. The two warnings below, written for
  # aggregations that skip missing values, are then not needed.
  if (aggregate_function == "SPEI") {

    years_with_window_gaps <- character(0)
    spei_gaps <- spei_gap_years(env_data, months_per_year)

    if (length(spei_gaps) > 0) {

      spei_switch <- method == "cor" &&
        !(cor_na_use %in% c("complete.obs", "na.or.complete",
                            "pairwise.complete.obs"))

      warning(spei_gap_message("env_data", spei_gaps, unit = "months",
                               switched = spei_switch,
                               old_na_use = cor_na_use,
                               detrending = !is.null(dc_method)))

      if (spei_switch) {
        cor_na_use <- "pairwise.complete.obs"
      }
    }
  }

  # Where a window is empty for some years, detrend_climate() fits on the
  # non-missing years and reinserts the gaps afterwards.
  if (!is.null(dc_method) && length(years_with_window_gaps) > 0) {
    warning(dc_gap_message("env_data", years_with_window_gaps,
                           min_window_width, unit = "months"))
  }

  # Missing values are not handled automatically. Tell the user explicitly,
  # so that an all-NA surface is never a silent surprise. Only gaps long
  # enough to empty a window matter, and only for correlations: lm and brnn
  # leave out incomplete years themselves.
  if (method == "cor" && length(years_with_window_gaps) > 0 &&
      !(cor_na_use %in% c("complete.obs", "na.or.complete",
                          "pairwise.complete.obs"))) {

    warning(na_gap_message("env_data", years_with_window_gaps,
                           min_window_width, unit = "months",
                           na_arg = "cor_na_use", na_value = cor_na_use))
  }

  # NA values are not allowed and must be removed from response data.frame
  # exception if cor_na_us accounts for missing values
  if (sum(is.na(response)) > 0 & !(cor_na_use %in% c("complete.obs", "na.or.complete", "pairwise.complete.obs"))){

    prob_year <- row.names(response[is.na(response), , drop = F])

    stop(paste0("NA is not allowed in response data frame. ",
                "Problematic year is ", prob_year))

  }

  # ---- shared window-coverage gate (method-agnostic) -----------------
  # A window is usable only if it has at least `min_complete` complete
  # response/predictor pairs. Used by cor, lm and brnn so all methods
  # return NA (not a silent partial fit) on fully/mostly missing windows.
  min_complete <- 4L
  window_is_usable <- function(x_vec, y_vec, min_n = min_complete) {
    sum(stats::complete.cases(x_vec, y_vec)) >= min_n
  }

  # PART 2 - Based on the selected function arguments, different chunks of code
  # will be used. For demonstration:
  # A) Chunks are used if fixed.withd != 0
  # A.1 method == "cor"
  # A.2 method == "lm"
  # A.3 method == "brnn"

  # B) Chunks are used if fixed.withd == 0
  # B.1 method == "cor"
  # B.2 method == "lm"
  # B.3 method == "brnn"

  # this is a list for climate and and holder for saving mm
  list_climate <- list()
  mm <- 1

  # Vectorised bootstrap. The resample indices are drawn once and shared by
  # every window, and each replicate correlation is obtained from matrix
  # products rather than from a resampling loop. Used only for method = "cor"
  # with interval types computable directly from the replicates; "bca" needs
  # jackknife estimates and keeps the boot() route.
  # bca intervals need at least as many replicates as resampled years
  bca_boot_n_warning(boot, boot_ci_type, boot_n, n_years = nrow(response))

  boot_fast_obj <- NULL
  if (boot == TRUE && method == "cor" &&
      boot_ci_type %in% c("norm", "perc", "basic")) {

    boot_fast_obj <- boot_fast_prep(y = response[, 1],
                                    boot_n = boot_n,
                                    cor_method = cor_method,
                                    cor_na_use = cor_na_use)
  }

  # A.1 method = "cor"
  if (fixed_width != 0 & method == "cor") {
    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }
    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be stored in this matrix


    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (fixed_width != max_window){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width - offset_end - offset_start + 1),
                             style = 3)}

    }
    b = 0

    # An iterating loop. In each itteration x is calculated and represents
    # response (dependent) variable. X is a moving average. Window width of
    # a moving window is fixed_width. Next, statistical metric is calculated
    # based on a selected method (cor, lm or brnn). Calculation is stored in
    # temporal matrix.
    for (j in (0 + offset_start -1): (ncol(env_data) - max((fixed_width + offset_end), offset_end))) {

      b = b + 1

      x <- aggregate_monthly_window(
        env_data, (1 + j):(j + fixed_width)
      )

      # Optional detrending of the aggregated climate series.
      x <- matrix(dc_apply(x, years = row.names(env_data)),
                  nrow = nrow(env_data), ncol = 1)

      x_list <- x
      colnames(x_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x_list) <- row.names(env_data)
      list_climate[[mm]] <- x_list
      mm = mm + 1


      if (boot == FALSE){

        temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method, use = cor_na_use)
        temporal_lower <- NA
        temporal_upper <- NA

      } else if (boot == TRUE){


        if (!is.null(boot_fast_obj)) {

          boot_t <- boot_fast_replicates(boot_fast_obj, x[, 1])
          boot_t0 <- suppressWarnings(
            cor(response[, 1], x[, 1], method = cor_method, use = cor_na_use))

          temporal_correlation <- boot_t0

          boot_ci_vals <- boot_fast_ci(boot_t, boot_t0,
                                       conf = boot_conf_int,
                                       type = boot_ci_type)
          temporal_lower <- boot_ci_vals[1]
          temporal_upper <- boot_ci_vals[2]
          ci_int <- NULL

        } else {

        temp_df_boot <- cbind(response[, 1], x[, 1])
        calc <- boot(temp_df_boot, boot_f_cor, cor.type = cor_method,
                     cor.use = cor_na_use, R = boot_n)

        temporal_correlation <- as.numeric(calc$t0[1])

        ci_int <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type), silent = TRUE)

        if (class(ci_int)[[1]] == "try-error"){

          temporal_lower <- NA
          temporal_upper <- NA

        } else {

          if (boot_ci_type == "norm"){

            temporal_lower <- ci_int$norm[2]
            temporal_upper <- ci_int$norm[3]

          } else if (boot_ci_type == "perc"){

            temporal_lower <- ci_int$perc[4]
            temporal_upper <- ci_int$perc[5]

          } else if (boot_ci_type == "stud") {

            temporal_lower <- ci_int$student[4]
            temporal_upper <- ci_int$student[5]

          } else if (boot_ci_type == "basic") {

            temporal_lower <- ci_int$basic[4]
            temporal_upper <- ci_int$basic[5]

          } else if (boot_ci_type == "bca") {

            temporal_lower <- ci_int$bca[4]
            temporal_upper <- ci_int$bca[5]

          } else {

            stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

          }
        }
        }

      } else {
        print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }

      # Each calculation is printed. Reason: usually it takes several minutes
      # to go through all loops and therefore, users might think that R is
      # not responding. But if each calculation is printed, user could be
      # confident, that R is responding.

      if (reference_window == 'start'){
        temporal_matrix[1, j + 1] <- temporal_correlation
        temporal_matrix_lower[1, j + 1] <- temporal_lower
        temporal_matrix_upper[1, j + 1] <- temporal_upper

      } else if (reference_window == 'end'){
        temporal_matrix[1, j + fixed_width] <- temporal_correlation
        temporal_matrix_lower[1, j + fixed_width] <- temporal_lower
        temporal_matrix_upper[1, j + fixed_width] <- temporal_upper

      } else if (reference_window == 'middle'){

        temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_correlation
        temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_lower
        temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_upper
      }

      if(interactive()){

        if (fixed_width != max_window){setTxtProgressBar(pb, b)}

      }
    }

    if(interactive()){

      if (fixed_width != max_window){close(pb)}

    }

    # temporal_matrix is given rownames and colnames. Rownames represent a
    # window width used fot calculations. Colnames represent the position of
    # moving window in a original env_data data frame.
    row.names(temporal_matrix) <- fixed_width
    row.names(temporal_matrix_lower) <- fixed_width
    row.names(temporal_matrix_upper) <- fixed_width

    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # A.2 method == "lm"
  # For a description see A.1
  if (fixed_width != 0 & method == "lm") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (fixed_width != max_window){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width - offset_end - offset_start + 1),
                             style = 3)}

    }

    b = 0

    for (j in (0 + offset_start -1): (ncol(env_data) - max((fixed_width + offset_end), offset_end))) {

      b = b + 1

      x <- aggregate_monthly_window(
        env_data, (1 + j):(j + fixed_width)
      )

      # Optional detrending of the aggregated climate series.
      x <- matrix(dc_apply(x, years = row.names(env_data)),
                  nrow = nrow(env_data), ncol = 1)

      x_list <- x
      colnames(x_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x_list) <- row.names(env_data)
      list_climate[[mm]] <- x_list
      mm = mm + 1

      if (!window_is_usable(x, response[, 1])) {

        temporal_r_squared <- NA
        temporal_adj_r_squared <- NA
        temporal_r_squared_lower <- NA
        temporal_r_squared_upper <- NA
        temporal_adj_r_squared_lower <- NA
        temporal_adj_r_squared_upper <- NA

      } else if (boot == FALSE){

        temporal_df <- data.frame(cbind(x, response))
        temporal_model <- lm(x ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        temporal_r_squared <- temporal_summary$r.squared
        temporal_adj_r_squared <- temporal_summary$adj.r.squared

        temporal_r_squared_lower<- NA
        temporal_r_squared_upper<- NA
        temporal_adj_r_squared_lower <- NA
        temporal_adj_r_squared_upper <- NA

      } else if (boot == TRUE){

        temporal_df <- data.frame(cbind(x, response))
        calc <- boot(data = temporal_df, statistic = boot_f_lm, R = boot_n, lm.formula = "x ~ .")

        # na.rm is used here, unlike the correlation path, because an NA
        # replicate means a model fit failed to converge rather than that
        # data were missing
        temporal_r_squared <- as.numeric(calc$t0[1])
        temporal_adj_r_squared <- as.numeric(calc$t0[2])

        ci_int_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1), silent = TRUE)
        ci_int_adj_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2), silent = TRUE)


        if (class(ci_int_r_squared)[[1]] == "try-error"){

          temporal_r_squared_lower<- NA
          temporal_r_squared_upper<- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else {

          if (boot_ci_type == "norm"){

            temporal_r_squared_lower <- ci_int_r_squared$norm[2]
            temporal_r_squared_upper <- ci_int_r_squared$norm[3]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$norm[2]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$norm[3]

          } else if (boot_ci_type == "perc"){

            temporal_r_squared_lower <- ci_int_r_squared$perc[4]
            temporal_r_squared_upper <- ci_int_r_squared$perc[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$perc[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$perc[5]

          } else if (boot_ci_type == "stud") {

            temporal_r_squared_lower <- ci_int_r_squared$studen[4]
            temporal_r_squared_upper <- ci_int_r_squared$studen[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$studen[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$studen[5]

          } else if (boot_ci_type == "basic") {

            temporal_r_squared_lower <- ci_int_r_squared$basic[4]
            temporal_r_squared_upper <- ci_int_r_squared$basic[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$basic[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$basic[5]

          } else if (boot_ci_type == "bca") {

            temporal_r_squared_lower <- ci_int_r_squared$bca[4]
            temporal_r_squared_upper <- ci_int_r_squared$bca[5]
            temporal_adj_r_squared_lower <- ci_int_adj_r_squared$bca[4]
            temporal_adj_r_squared_upper <- ci_int_adj_r_squared$bca[5]

          } else {

            stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

          }

        }

      } else {
        stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }

      if (metric == "r.squared"){

        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_r_squared_upper

        } else if (reference_window == 'end') {

          temporal_matrix[1, j + fixed_width] <- temporal_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_upper

        }
      }

      if (metric == "adj.r.squared"){
        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_adj_r_squared_upper

        } else if (reference_window == 'end'){

          temporal_matrix[1, j + fixed_width] <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_adj_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_upper
        }
      }

      if(interactive()){

        if (fixed_width != max_window){setTxtProgressBar(pb, b)}

      }
    }

    if(interactive()){

      if (fixed_width != max_window){close(pb)}

    }

    row.names(temporal_matrix) <- fixed_width
    row.names(temporal_matrix_lower) <- fixed_width
    row.names(temporal_matrix_upper) <- fixed_width

    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # A.3 method == "brnn"
  # For a description see A.1
  if (fixed_width != 0 & method == "brnn") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = (ncol(env_data) - fixed_width) + 1)
    } else if (reference_window == 'end') {
      temporal_matrix <- matrix(NA, nrow = 1, ncol = (ncol(env_data)))
    } else if (reference_window == 'middle') {
      temporal_matrix <- matrix(NA, nrow = 1,
                                ncol = round2((ncol(env_data) - fixed_width +
                                                 1 + fixed_width/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (fixed_width != max_window){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (ncol(env_data) - fixed_width - offset_end - offset_start + 1),
                             style = 3)}
    }

    b = 0

    for (j in (0 + offset_start -1): (ncol(env_data) - max((fixed_width + offset_end), offset_end))) {

      b = b + 1

      x <- aggregate_monthly_window(
        env_data, (1 + j):(j + fixed_width)
      )

      # Optional detrending of the aggregated climate series.
      x <- matrix(dc_apply(x, years = row.names(env_data)),
                  nrow = nrow(env_data), ncol = 1)

      x_list <- x
      colnames(x_list) <- paste0(j + 1, "_" ,j + fixed_width)
      row.names(x_list) <- row.names(env_data)
      list_climate[[mm]] <- x_list
      mm = mm + 1

      if (boot == FALSE){

        temporal_df <- data.frame(cbind(x, response))

        brnn_result <- safe_brnn_metrics(
          temporal_df = temporal_df,
          formula = x ~ .,
          observed = temporal_df[, 1],
          neurons = neurons
        )

        temporal_r_squared <- brnn_result$r.squared
        temporal_adj_r_squared <- brnn_result$adj.r.squared

        temporal_r_squared_lower <- NA
        temporal_r_squared_upper <- NA
        temporal_adj_r_squared_lower <- NA
        temporal_adj_r_squared_upper <- NA

      } else if (boot == TRUE){

        temporal_df <- data.frame(cbind(x, response))

        brnn_boot_result <- safe_brnn_boot_metrics(
          temporal_df = temporal_df,
          neurons = neurons,
          boot_n = boot_n,
          boot_conf_int = boot_conf_int,
          boot_ci_type = boot_ci_type
        )

        temporal_r_squared <- brnn_boot_result$r.squared
        temporal_adj_r_squared <- brnn_boot_result$adj.r.squared

        temporal_r_squared_lower <- brnn_boot_result$r.squared.lower
        temporal_r_squared_upper <- brnn_boot_result$r.squared.upper
        temporal_adj_r_squared_lower <- brnn_boot_result$adj.r.squared.lower
        temporal_adj_r_squared_upper <- brnn_boot_result$adj.r.squared.upper

      } else {

        stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))

      }



      if (metric == "r.squared"){

        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_r_squared_upper

        } else if (reference_window == 'end') {

          temporal_matrix[1, j + fixed_width] <- temporal_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_r_squared_upper

        }
      }

      if (metric == "adj.r.squared"){
        if (reference_window == 'start'){

          temporal_matrix[1, j + 1]  <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + 1]  <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + 1]  <- temporal_adj_r_squared_upper

        } else if (reference_window == 'end'){

          temporal_matrix[1, j + fixed_width] <- temporal_adj_r_squared
          temporal_matrix_lower[1, j + fixed_width] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, j + fixed_width] <- temporal_adj_r_squared_upper

        } else if (reference_window == 'middle'){

          temporal_matrix[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared
          temporal_matrix_lower[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_lower
          temporal_matrix_upper[1, round2(j + 1 + fixed_width/2, 0)] <- temporal_adj_r_squared_upper
        }
      }
      if(interactive()){

        if (fixed_width != max_window){setTxtProgressBar(pb, b)}

      }

    }

    if(interactive()){

      if (fixed_width != max_window){close(pb)}

    }

    row.names(temporal_matrix) <- fixed_width
    row.names(temporal_matrix_lower) <- fixed_width
    row.names(temporal_matrix_upper) <- fixed_width

    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames

  }

  # B fixed_width == 0, in this case, lower_limit and upper_limit arguments
  # will be used to define window width of a moving window.
  # B.1 method == "cor"

  if (fixed_width == 0 & method == "cor") {

    # This is an empty matrix, currently filled with NA's
    # Latter, calculations will be stored in this matrix

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix


    # An iterating double loop: 1 outer loop) iterating from lower_limit :
    # upper_limit defines windo.width used for a moving window. 2) inner loop
    # defines the starting position of a moving window.
    # In each itteration, x is calculated and represents a response (dependent)
    # variable. x is a moving average, based on rowMeans/apply function.
    # Next, statistical metric is calculated based on a selected method (cor,
    # lm or brnn). Calculation is stored in temporal matrix in a proper place.
    # The position of stored calculation is informative later used for
    # indiciating optimal values.

    if (upper_limit != lower_limit){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)

      }
    }

    b = 0


    for (K in lower_limit:upper_limit) {

      b = b + 1

      for (j in (0 + offset_start -1): (ncol(env_data) - max((K + offset_end), offset_end))) {

        x <- aggregate_monthly_window(
          env_data, (1 + j):(j + K)
        )

        # Optional detrending of the aggregated climate series.
        x <- matrix(dc_apply(x, years = row.names(env_data)),
                    nrow = nrow(env_data), ncol = 1)

        x_list <- x
        colnames(x_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x_list) <- row.names(env_data)
        list_climate[[mm]] <- x_list
        mm = mm + 1

        if (boot == FALSE){
          temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method, use = cor_na_use)
          temporal_lower <- NA
          temporal_upper <- NA
        } else if (boot == TRUE){
          if (!is.null(boot_fast_obj)) {

            boot_t <- boot_fast_replicates(boot_fast_obj, x[, 1])
            boot_t0 <- suppressWarnings(
              cor(response[, 1], x[, 1], method = cor_method, use = cor_na_use))

            temporal_correlation <- boot_t0

            boot_ci_vals <- boot_fast_ci(boot_t, boot_t0,
                                         conf = boot_conf_int,
                                         type = boot_ci_type)
            temporal_lower <- boot_ci_vals[1]
            temporal_upper <- boot_ci_vals[2]
            ci_int <- NULL

          } else {

          temp_df_boot <- cbind(response[, 1], x[, 1])
          calc <- boot(temp_df_boot, boot_f_cor, cor.type = cor_method,
                       cor.use = cor_na_use, R = boot_n)

          temporal_correlation <- as.numeric(calc$t0[1])

          ci_int <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type), silent = TRUE)

          if (class(ci_int)[[1]] == "try-error"){

            temporal_lower <- NA
            temporal_upper <- NA

          } else {

            if (boot_ci_type == "norm"){

              temporal_lower <- ci_int$norm[2]
              temporal_upper <- ci_int$norm[3]

            } else if (boot_ci_type == "perc"){

              temporal_lower <- ci_int$perc[4]
              temporal_upper <- ci_int$perc[5]

            } else if (boot_ci_type == "stud") {

              temporal_lower <- ci_int$student[4]
              temporal_upper <- ci_int$student[5]

            } else if (boot_ci_type == "basic") {

              temporal_lower <- ci_int$basic[4]
              temporal_upper <- ci_int$basic[5]

            } else if (boot_ci_type == "bca") {

              temporal_lower <- ci_int$bca[4]
              temporal_upper <- ci_int$bca[5]

            } else {

              stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

            }

          }
          }

        } else {
          print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
        }

        if (reference_window == 'start'){
          temporal_matrix[(K - lower_limit) + 1, j + 1] <- temporal_correlation
          temporal_matrix_lower[(K - lower_limit) + 1, j + 1] <- temporal_lower
          temporal_matrix_upper[(K - lower_limit) + 1, j + 1] <- temporal_upper
        } else if (reference_window == 'end'){
          temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_correlation
          temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_lower
          temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_upper
        } else if (reference_window == 'middle'){
          temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_correlation
          temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_lower
          temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_upper
        }

      }

      if(interactive()){

        if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}
      }
    }

    if(interactive()){

      if (upper_limit != lower_limit){close(pb)}

    }

    # temporal_matrix is given rownames and colnames. Rownames represent a
    # window width used fot calculations. Colnames represent the position of
    # moving window in a original env_data data frame.
    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit,
                                       by = 1))
    row.names(temporal_matrix) <- temporal_rownames
    row.names(temporal_matrix_lower) <- temporal_rownames
    row.names(temporal_matrix_upper) <- temporal_rownames


    temporal_colnames <- as.vector(seq(from = 1,
                                       to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # B.2 method == "lm"
  # For a description see B.1
  if (fixed_width == 0 & method == "lm") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (upper_limit != lower_limit){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)

      }

    }

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1

      for (j in (0 + offset_start -1): (ncol(env_data) - max((K + offset_end), offset_end))) {

        x <- aggregate_monthly_window(
          env_data, (1 + j):(j + K)
        )

        # Optional detrending of the aggregated climate series.
        x <- matrix(dc_apply(x, years = row.names(env_data)),
                    nrow = nrow(env_data), ncol = 1)

        x_list <- x
        colnames(x_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x_list) <- row.names(env_data)
        list_climate[[mm]] <- x_list
        mm = mm + 1

        if (!window_is_usable(x, response[, 1])) {

          temporal_r_squared <- NA
          temporal_adj_r_squared <- NA
          temporal_r_squared_lower <- NA
          temporal_r_squared_upper <- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else if (boot == FALSE){

          temporal_df <- data.frame(cbind(x, response))
          temporal_model <- lm(x ~ ., data = temporal_df)
          temporal_summary <- summary(temporal_model)
          temporal_r_squared <- temporal_summary$r.squared
          temporal_adj_r_squared <- temporal_summary$adj.r.squared

          temporal_r_squared_lower<- NA
          temporal_r_squared_upper<- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else if (boot == TRUE){

          temporal_df <- data.frame(cbind(x, response))
          calc <- boot(data = temporal_df, statistic = boot_f_lm, R = boot_n, lm.formula = "x ~ .")

          # na.rm is used here, unlike the correlation path, because an NA
          # replicate means a model fit failed to converge rather than that
          # data were missing
          temporal_r_squared <- as.numeric(calc$t0[1])
          temporal_adj_r_squared <- as.numeric(calc$t0[2])

          ci_int_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 1), silent = TRUE)
          ci_int_adj_r_squared <- try(boot.ci(calc, conf = boot_conf_int, type = boot_ci_type, index = 2), silent = TRUE)

          if (class(ci_int_r_squared)[[1]] == "try-error"){

            temporal_r_squared_lower<- NA
            temporal_r_squared_upper<- NA
            temporal_adj_r_squared_lower <- NA
            temporal_adj_r_squared_upper <- NA

          } else {

            if (boot_ci_type == "norm"){

              temporal_r_squared_lower <- ci_int_r_squared$norm[2]
              temporal_r_squared_upper <- ci_int_r_squared$norm[3]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$norm[2]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$norm[3]

            } else if (boot_ci_type == "perc"){

              temporal_r_squared_lower <- ci_int_r_squared$perc[4]
              temporal_r_squared_upper <- ci_int_r_squared$perc[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$perc[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$perc[5]

            } else if (boot_ci_type == "stud") {

              temporal_r_squared_lower <- ci_int_r_squared$studen[4]
              temporal_r_squared_upper <- ci_int_r_squared$studen[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$studen[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$studen[5]

            } else if (boot_ci_type == "basic") {

              temporal_r_squared_lower <- ci_int_r_squared$basic[4]
              temporal_r_squared_upper <- ci_int_r_squared$basic[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$basic[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$basic[5]

            } else if (boot_ci_type == "bca") {

              temporal_r_squared_lower <- ci_int_r_squared$bca[4]
              temporal_r_squared_upper <- ci_int_r_squared$bca[5]
              temporal_adj_r_squared_lower <- ci_int_adj_r_squared$bca[4]
              temporal_adj_r_squared_upper <- ci_int_adj_r_squared$bca[5]

            } else {

              stop("boot_ci_type should be 'norm', 'perc', 'stud', 'basic' or 'bca'")

            }

          }

        } else {
          stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))
        }

        if (metric == "r.squared"){

          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_upper
          }

        }

        if (metric == "adj.r.squared"){


          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_upper
          }

        }
      }

      if(interactive()){

        if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}

      }
    }

    if(interactive()){

      if (upper_limit != lower_limit){close(pb)}

    }

    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit, by = 1))
    row.names(temporal_matrix) <- temporal_rownames
    row.names(temporal_matrix_lower) <- temporal_rownames
    row.names(temporal_matrix_upper) <- temporal_rownames

    temporal_colnames <- as.vector(seq(from = 1, to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # B.3 method == "brnn"
  # For a description see B.1
  if (fixed_width == 0 & method == "brnn") {

    if (reference_window == 'start'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data) - lower_limit) + 1)
    } else if (reference_window == 'end'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = (ncol(env_data)))
    } else if (reference_window == 'middle'){
      temporal_matrix <- matrix(NA, nrow = (upper_limit - lower_limit + 1),
                                ncol = round2((ncol(env_data) - lower_limit +
                                                 1 + lower_limit/2 ),0))
    }

    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix

    if (upper_limit != lower_limit){

      if(interactive()){

        pb <- txtProgressBar(min = 0, max = (upper_limit - lower_limit), style = 3)

      }
    }

    b = 0

    for (K in lower_limit:upper_limit) {

      b = b + 1


      for (j in (0 + offset_start -1): (ncol(env_data) - max((K + offset_end), offset_end))) {

        x <- aggregate_monthly_window(
          env_data, (1 + j):(j + K)
        )

        # Optional detrending of the aggregated climate series.
        x <- matrix(dc_apply(x, years = row.names(env_data)),
                    nrow = nrow(env_data), ncol = 1)

        x_list <- x
        colnames(x_list) <- paste0(j + 1, "_" ,j + K)
        row.names(x_list) <- row.names(env_data)
        list_climate[[mm]] <- x_list
        mm = mm + 1

        if (boot == FALSE){

          temporal_df <- data.frame(cbind(x, response))

          brnn_result <- safe_brnn_metrics(
            temporal_df = temporal_df,
            formula = x ~ .,
            observed = temporal_df[, 1],
            neurons = neurons
          )

          temporal_r_squared <- brnn_result$r.squared
          temporal_adj_r_squared <- brnn_result$adj.r.squared

          temporal_r_squared_lower <- NA
          temporal_r_squared_upper <- NA
          temporal_adj_r_squared_lower <- NA
          temporal_adj_r_squared_upper <- NA

        } else if (boot == TRUE){

          temporal_df <- data.frame(cbind(x, response))

          brnn_boot_result <- safe_brnn_boot_metrics(
            temporal_df = temporal_df,
            neurons = neurons,
            boot_n = boot_n,
            boot_conf_int = boot_conf_int,
            boot_ci_type = boot_ci_type
          )

          temporal_r_squared <- brnn_boot_result$r.squared
          temporal_adj_r_squared <- brnn_boot_result$adj.r.squared

          temporal_r_squared_lower <- brnn_boot_result$r.squared.lower
          temporal_r_squared_upper <- brnn_boot_result$r.squared.upper
          temporal_adj_r_squared_lower <- brnn_boot_result$adj.r.squared.lower
          temporal_adj_r_squared_upper <- brnn_boot_result$adj.r.squared.upper

        } else {

          stop(paste0("boot should be TRUE or FALSE, instead it is ", boot))

        }

        if (metric == "r.squared"){

          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_r_squared_upper
          }

        }

        if (metric == "adj.r.squared"){


          if (reference_window == 'start'){
            temporal_matrix[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + 1]  <-  temporal_adj_r_squared_upper

          } else if (reference_window == 'end') {

            temporal_matrix[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, j + K] <- temporal_adj_r_squared_upper

          } else if (reference_window == 'middle'){

            temporal_matrix[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared
            temporal_matrix_lower[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_lower
            temporal_matrix_upper[(K - lower_limit) + 1, round2(j + 1 + K/2, 0)] <- temporal_adj_r_squared_upper
          }

        }

      }

      if(interactive()){

        if (upper_limit != lower_limit){setTxtProgressBar(pb, b)}

      }
    }

    if(interactive()){

      if (upper_limit != lower_limit){close(pb)}

    }

    temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit, by = 1))
    row.names(temporal_matrix) <- temporal_rownames
    row.names(temporal_matrix_lower) <- temporal_rownames
    row.names(temporal_matrix_upper) <- temporal_rownames

    temporal_colnames <- as.vector(seq(from = 1, to = ncol(temporal_matrix), by = 1))
    colnames(temporal_matrix) <- temporal_colnames
    colnames(temporal_matrix_lower) <- temporal_colnames
    colnames(temporal_matrix_upper) <- temporal_colnames
  }

  # PART 3: smoothing function, if brnn method is used. It turnes out, that
  # brnn sometimes (1 - 3 % of calculations) fails to construct a realistic
  # result. In those cases, much lower r.squared (or adj.r.squared) are
  # calculated. smooth_matrix function removes unrealistic calculations and
  # replace them with mean of values in a window 3 x 3. Maximum value
  # is not affected by any means.
  # [i - 1, j - 1], [i - 1, j], [i - 1, j + 1]
  # [    i, j - 1], [    i, j], [    i, j + 1]
  # [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]
  if (method == "brnn" & brnn_smooth == TRUE){
    temporal_matrix <- smooth_matrix(temporal_matrix, factor_drop = 0.7,
                                     repeats = 2)
  }

  # To enhance the visualisation, insignificant values are removed if
  # remove_insignificant == TRUE. The criterion is set by
  # remove_insignificant_method:
  #   "alpha"  per-window significance threshold; does not account for how
  #            many windows were screened
  #   "boot"   bootstrap confidence interval excludes zero (per-window)
  #   "fdr_BY" false discovery rate, Benjamini-Yekutieli, valid under
  #            arbitrary dependence between windows
  #   "fdr_BH" false discovery rate, Benjamini-Hochberg, valid under positive
  #            dependence and less conservative than BY
  significance_filter <- list(
    applied = FALSE,
    method = NA_character_,
    criterion = NA_character_,
    removed_cells = NA_integer_,
    retained_cells = NA_integer_,
    total_cells = length(temporal_matrix)
  )

  if (isTRUE(remove_insignificant)) {

    n_obs <- nrow(response)
    keep  <- NULL

    # 1 Plain alpha threshold
    if (remove_insignificant_method == "alpha") {

      critical_threshold_cor <- critical_r(n_obs, alpha = alpha,
                                           cor_method = cor_method)
      critical_threshold_cor2 <- critical_threshold_cor ^ 2

      if (method == "cor") {

        keep <- abs(temporal_matrix) >= abs(critical_threshold_cor)

      } else if (method == "lm" | method == "brnn") {

        keep <- abs(temporal_matrix) >= abs(critical_threshold_cor2)
      }

      significance_filter$applied <- TRUE
      significance_filter$method <- "alpha"
      significance_filter$criterion <- paste0("|r| >= ",
                                              round(critical_threshold_cor, 3),
                                              " (alpha = ", alpha, ")")

      significance_filter$critical_r <- critical_threshold_cor


      # 2 Bootstrap confidence interval
    } else if (remove_insignificant_method == "boot") {

      if (is.null(temporal_matrix_lower) ||
          is.null(temporal_matrix_upper) ||
          !all(dim(temporal_matrix) == dim(temporal_matrix_lower)) ||
          !all(dim(temporal_matrix) == dim(temporal_matrix_upper))) {

        warning(paste0("remove_insignificant_method = 'boot' was requested, but ",
                       "bootstrap confidence-interval matrices are missing or ",
                       "do not match the calculation matrix. The significance ",
                       "filter is ignored."))

      } else {

        keep <- (temporal_matrix_lower > 0 & temporal_matrix_upper > 0) |
          (temporal_matrix_lower < 0 & temporal_matrix_upper < 0)

        significance_filter$applied <- TRUE
        significance_filter$method <- "boot"
        significance_filter$criterion <- "Bootstrap confidence interval excludes zero"
      }

      # 3 False discovery rate
    } else if (remove_insignificant_method == "fdr_BY" |
               remove_insignificant_method == "fdr_BH") {

      fdr_type <- sub("^fdr_", "", remove_insignificant_method)

      if (method == "cor") {

        r_surface <- temporal_matrix

      } else if (method == "lm" | method == "brnn") {

        r_surface <- sqrt(abs(temporal_matrix))
      }

      p_surface <- p_from_r(r_surface, n_obs,
                            cor_method = cor_method)

      p_adjusted <- stats::p.adjust(as.vector(p_surface), method = fdr_type)
      p_adjusted <- matrix(p_adjusted, nrow = nrow(p_surface))

      keep <- p_adjusted < alpha

      significance_filter$applied <- TRUE
      significance_filter$method <- remove_insignificant_method
      significance_filter$criterion <- paste0("Adjusted p < ", alpha,
                                              " (", fdr_type, ")")

      # 4 Max-statistic permutation
    } else if (remove_insignificant_method == "maxstat") {

      climate_matrix <- do.call(cbind, list_climate)

      if (is.null(climate_matrix) || ncol(climate_matrix) < 1) {

        warning(paste0("remove_insignificant_method = 'maxstat' was requested, ",
                       "but the matrix of aggregated climate windows is not ",
                       "available. The significance filter is ignored."))

      } else {

        shared_years <- intersect(row.names(response), row.names(climate_matrix))

        if (length(shared_years) < 10) {

          warning(paste0("remove_insignificant_method = 'maxstat' was requested, ",
                         "but fewer than 10 years are shared between the response ",
                         "and the aggregated climate windows. The significance ",
                         "filter is ignored."))

        } else {

          y_perm <- response[shared_years, 1]
          x_perm <- climate_matrix[shared_years, , drop = FALSE]

          # Null distribution of the largest absolute correlation obtained
          # anywhere in the search when the response is unrelated to climate.
          null_max <- rep(NA_real_, perm_n)

          for (perm_i in seq_len(perm_n)) {

            perm_cor <- suppressWarnings(
              stats::cor(sample(y_perm), x_perm,
                         method = cor_method, use = cor_na_use))

            null_max[perm_i] <- max(abs(perm_cor), na.rm = TRUE)
          }

          critical_threshold_perm <- as.numeric(
            stats::quantile(null_max, probs = 1 - alpha, na.rm = TRUE))

          # For "lm" the stored statistic is R squared, so the threshold is
          # squared as well. This is exact for a single predictor.
          if (method == "cor") {

            observed_max <- max(abs(temporal_matrix), na.rm = TRUE)
            keep <- abs(temporal_matrix) >= critical_threshold_perm

          } else if (method == "lm") {

            observed_max <- sqrt(max(abs(temporal_matrix), na.rm = TRUE))
            keep <- abs(temporal_matrix) >= critical_threshold_perm ^ 2
          }

          p_global <- (sum(null_max >= observed_max, na.rm = TRUE) + 1) /
            (perm_n + 1)

          significance_filter$applied <- TRUE
          significance_filter$method <- "maxstat"
          significance_filter$criterion <- paste0("|r| >= ",
                                                  round(critical_threshold_perm, 3),
                                                  " (max-statistic permutation, ",
                                                  perm_n, " permutations, alpha = ",
                                                  alpha, ")")
          significance_filter$critical_r <- critical_threshold_perm
          significance_filter$observed_max_r <- observed_max
          significance_filter$p_global <- p_global
          significance_filter$perm_n <- perm_n
          significance_filter$n_windows_searched <- ncol(x_perm)
        }
      }

      # 5 Sidak correction on the effective number of tests
    } else if (remove_insignificant_method == "sidak") {

      climate_matrix <- do.call(cbind, list_climate)

      if (is.null(climate_matrix) || ncol(climate_matrix) < 2) {

        warning(paste0("remove_insignificant_method = 'sidak' was requested, ",
                       "but the matrix of aggregated climate windows is not ",
                       "available. The significance filter is ignored."))

      } else {

        shared_years <- intersect(row.names(response), row.names(climate_matrix))

        if (length(shared_years) < 10) {

          warning(paste0("remove_insignificant_method = 'sidak' was requested, ",
                         "but fewer than 10 years are shared between the ",
                         "response and the aggregated climate windows. The ",
                         "significance filter is ignored."))

        } else {

          n_obs_sidak <- length(shared_years)

          x_sidak <- climate_matrix[shared_years, , drop = FALSE]

          # Missing years are left out rather than whole windows, so that the
          # effective number of tests describes every screened window, see
          # sidak_window_matrix(). Complete data are not affected.
          x_sidak <- sidak_window_matrix(x_sidak)
          x_sidak <- x_sidak[, apply(x_sidak, 2, stats::sd) > 0, drop = FALSE]

          if (ncol(x_sidak) < 2) {

            warning(paste0("remove_insignificant_method = 'sidak' was ",
                           "requested, but fewer than two complete climate ",
                           "windows remain. The significance filter is ",
                           "ignored."))

          } else {

            # Effective number of independent tests, following Li and Ji
            # (2005). The correlation matrix among the screened windows has
            # rank at most the number of years, so its non-zero eigenvalues
            # are obtained from the much smaller years-by-years matrix
            # instead of the windows-by-windows one.
            z_sidak <- scale(x_sidak)
            gram <- tcrossprod(z_sidak) / (nrow(z_sidak) - 1)

            eig <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
            eig <- eig[eig > 1e-10]

            m_eff <- sum((eig >= 1) + (eig - floor(eig)))

            alpha_sidak <- 1 - (1 - alpha) ^ (1 / m_eff)
            critical_threshold_sidak <- critical_r(n_obs_sidak,
                                                   alpha = alpha_sidak,
                                                   cor_method = cor_method)

            if (method == "cor") {

              keep <- abs(temporal_matrix) >= critical_threshold_sidak

            } else if (method == "lm" | method == "brnn") {

              keep <- abs(temporal_matrix) >= critical_threshold_sidak ^ 2
            }

            significance_filter$applied <- TRUE
            significance_filter$method <- "sidak"
            significance_filter$criterion <- paste0("|r| >= ",
                                                    round(critical_threshold_sidak, 3),
                                                    " (Sidak on ", round(m_eff, 1),
                                                    " effective tests, alpha = ",
                                                    alpha, ")")
            significance_filter$critical_r <- critical_threshold_sidak
            significance_filter$effective_tests <- m_eff
            significance_filter$alpha_per_test <- alpha_sidak
            significance_filter$n_windows_searched <- ncol(x_sidak)
          }
        }
      }

      # 6 Cluster-extent inference
    } else if (remove_insignificant_method == "cluster") {

      climate_matrix <- do.call(cbind, list_climate)

      if (is.null(climate_matrix) || ncol(climate_matrix) < 1) {

        warning(paste0("remove_insignificant_method = 'cluster' was requested, ",
                       "but the matrix of aggregated climate windows is not ",
                       "available. The significance filter is ignored."))

      } else {

        shared_years <- intersect(row.names(response), row.names(climate_matrix))

        if (length(shared_years) < 10) {

          warning(paste0("remove_insignificant_method = 'cluster' was ",
                         "requested, but fewer than 10 years are shared ",
                         "between the response and the aggregated climate ",
                         "windows. The significance filter is ignored."))

        } else {

          y_perm <- response[shared_years, 1]
          x_perm <- climate_matrix[shared_years, , drop = FALSE]
          n_obs_cl <- length(shared_years)

          # Unlike the max-statistic test, which needs only the largest
          # permuted correlation, cluster inference needs the permuted values
          # laid out on the same surface as the observed ones. Each aggregated
          # window is named "start_end", which gives its width and position;
          # the mapping is checked against the observed surface before use.
          win_nm <- colnames(x_perm)
          parts <- strsplit(win_nm, "_", fixed = TRUE)
          w_start <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 1)))
          w_end <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
          w_width <- w_end - w_start + 1

          surf_widths <- suppressWarnings(as.numeric(row.names(temporal_matrix)))

          if (nrow(temporal_matrix) == 1) {
            row_idx <- rep(1L, length(w_width))
          } else {
            row_idx <- match(w_width, surf_widths)
          }

          col_idx <- if (reference_window == "start") {
            w_start
          } else if (reference_window == "end") {
            w_end
          } else {
            round2(w_start + w_width / 2, 0)
          }

          usable <- !is.na(row_idx) & !is.na(col_idx) &
            row_idx >= 1 & row_idx <= nrow(temporal_matrix) &
            col_idx >= 1 & col_idx <= ncol(temporal_matrix)

          map_ok <- FALSE
          if (any(usable)) {
            obs_cor <- suppressWarnings(
              as.numeric(stats::cor(y_perm, x_perm[, usable, drop = FALSE],
                                    method = cor_method, use = cor_na_use)))
            on_surface <- temporal_matrix[cbind(row_idx[usable],
                                                col_idx[usable])]
            both <- is.finite(obs_cor) & is.finite(on_surface)
            map_ok <- sum(both) > 0 &&
              max(abs(obs_cor[both] - on_surface[both])) < 1e-8
          }

          if (!map_ok) {

            warning(paste0("remove_insignificant_method = 'cluster' was ",
                           "requested, but the aggregated windows could not be ",
                           "matched to their positions on the response ",
                           "surface. The significance filter is ignored."))

          } else {

            idx_cells <- cbind(row_idx[usable], col_idx[usable])
            x_use <- x_perm[, usable, drop = FALSE]

            form_thr <- if (identical(cluster_threshold, "nominal")) {
              critical_r(n_obs_cl, alpha = alpha, cor_method = cor_method)
            } else {
              cluster_threshold
            }

            null_clust <- rep(NA_real_, perm_n)
            null_max <- rep(NA_real_, perm_n)

            perm_surface <- matrix(NA_real_, nrow(temporal_matrix),
                                   ncol(temporal_matrix))

            for (perm_i in seq_len(perm_n)) {

              perm_cor <- suppressWarnings(
                as.numeric(stats::cor(sample(y_perm), x_use,
                                      method = cor_method, use = cor_na_use)))

              perm_surface[] <- NA_real_
              perm_surface[idx_cells] <- perm_cor

              null_max[perm_i] <- max(abs(perm_cor), na.rm = TRUE)
              null_clust[perm_i] <- max(cluster_fast_masses(perm_surface,
                                                            form_thr))
            }

            crit_mass <- as.numeric(
              stats::quantile(null_clust, probs = 1 - alpha, na.rm = TRUE))
            crit_max_cl <- as.numeric(
              stats::quantile(null_max, probs = 1 - alpha, na.rm = TRUE))

            obs_mask <- !is.na(temporal_matrix) &
              abs(temporal_matrix) > form_thr
            cc <- cluster_fast_label(obs_mask)

            keep <- matrix(FALSE, nrow(temporal_matrix),
                           ncol(temporal_matrix))
            n_retained <- 0L
            obs_masses <- numeric(0)

            if (cc$n > 0L) {
              for (cl in seq_len(cc$n)) {
                in_cl <- cc$labels == cl
                mass <- sum(abs(temporal_matrix[in_cl]), na.rm = TRUE)
                obs_masses <- c(obs_masses, mass)
                if (mass > crit_mass) {
                  keep[in_cl] <- TRUE
                  n_retained <- n_retained + 1L
                }
              }
            }

            significance_filter$applied <- TRUE
            significance_filter$method <- "cluster"
            significance_filter$criterion <- paste0("cluster mass > ",
                                                    round(crit_mass, 1),
                                                    ", clusters formed at |r| > ",
                                                    round(form_thr, 3), " (",
                                                    perm_n,
                                                    " permutations, alpha = ",
                                                    alpha, ")")
            significance_filter$forming_threshold <- form_thr
            significance_filter$cluster_mass_critical <- crit_mass
            significance_filter$n_clusters_formed <- cc$n
            significance_filter$n_clusters_retained <- n_retained
            significance_filter$largest_cluster_mass <-
              if (length(obs_masses)) max(obs_masses) else 0
            significance_filter$maxstat_critical_r <- crit_max_cl
            significance_filter$perm_n <- perm_n
            significance_filter$n_windows_searched <- ncol(x_use)
          }
        }
      }
    }

    # Apply whichever mask was produced
    if (!is.null(keep)) {

      keep[is.na(keep)] <- FALSE

      temporal_matrix[!keep] <- NA

      # The same mask is applied to the confidence-interval matrices, so that
      # a filtered surface is not accompanied by intervals for windows the
      # filter removed.
      if (exists("temporal_matrix_lower") && !is.null(temporal_matrix_lower) &&
          all(dim(temporal_matrix_lower) == dim(temporal_matrix))) {
        temporal_matrix_lower[!keep] <- NA
        temporal_matrix_upper[!keep] <- NA
      }

      significance_filter$removed_cells  <- sum(!keep, na.rm = TRUE)
      significance_filter$retained_cells <- sum(keep, na.rm = TRUE)
      significance_filter$total_cells    <- length(keep)
    }
  }


  if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){

    if (isTRUE(significance_filter$applied) &&
        identical(significance_filter$method, "maxstat")) {

      warning(paste0("No window exceeded the max-statistic permutation ",
                     "threshold (global p = ",
                     signif(significance_filter$p_global, 3),
                     "). This indicates that the strongest correlation found ",
                     "anywhere in the search is no larger than expected when ",
                     "the response is unrelated to climate."))

    } else if (isTRUE(significance_filter$applied)) {
      warning(paste0("All calculations were removed by the significance filter (",
                     significance_filter$method, "). Change alpha or ",
                     "remove_insignificant_method."))
    } else {
      warning("All calculations are insignificant! Please change the alpha argument.")
    }

    final_list <- list(calculations = temporal_matrix,
                       method = method,
                       metric = cor_method,
                       analysed_period = NA,
                       optimized_return = NA,
                       optimized_return_all = NA,
                       transfer_function = NA, temporal_stability = NA,
                       cross_validation = NA,
                       plot_heatmap = NA,
                       plot_extreme = NA,
                       type = "monthly",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       significance_filter = significance_filter,
                       aggregated_climate = NA,
                       detrending = detrending_summary)

    class(final_list) <- 'dmrs'

    return(final_list)

  } else {



    ########################################################################
    # PART 4: Final list is being created and returned as a function output#
    ########################################################################

    # The first three elements of the final list are already created: calculated
    # values, method and metric used.
    # Here we create the fourth element: the optimal sequence of days that
    # returns the best selected statistical metric. We name it optimal_return.

    # In case of negative correlations, different strategy is applied.
    # For more detailed description see plot_extreme()

    if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){
      stop("All calculations are insignificant! Change the alpha argument!")
    }

    overall_max <- max(temporal_matrix, na.rm = TRUE)
    overall_min <- min(temporal_matrix, na.rm = TRUE)

    # absolute vales of overall_maximum and overall_minimum are compared and
    # one of the following two if functions is used
    # There are unimportant warnings produced:
    # no non-missing arguments to max; returning -Inf

    if ((abs(overall_max) >= abs(overall_min)) == TRUE) {

      # maximum value is located. Row indeces are needed to query information
      # about the window width used to calculate the maximum. Column name is
      # needed to query the starting day.
      max_result <- suppressWarnings(which.max(apply(temporal_matrix,
                                                     MARGIN = 2, max,
                                                     na.rm = TRUE)))
      plot_column <- max_result
      max_index <- which.max(temporal_matrix[, names(max_result)])
      row_index <- row.names(temporal_matrix)[max_index]
    }

    if ((abs(overall_max) < abs(overall_min)) == TRUE) {

      min_result <- suppressWarnings(which.min(apply(temporal_matrix,
                                                     MARGIN = 2, min,
                                                     na.rm = TRUE)))
      plot_column <- min_result
      min_index <- which.min(temporal_matrix[, names(min_result)])
      row_index <- row.names(temporal_matrix)[min_index]
    }

    # The fourth return element is being created: rowMeans/ apply of optimal sequence:
    # So, here we consider more options, based on the reference_winow
    # 1. reference window = "start"
    if (reference_window == 'start'){


      dataf <- data.frame(aggregate_monthly_window(
        env_data,
        as.numeric(plot_column):
          (as.numeric(plot_column) + as.numeric(row_index) - 1)
      ))

      # if detrending was applied, should also be applied here
      dataf[, 1] <- dc_apply(dataf[, 1], years = row.names(dataf))

      dataf_full <- cbind(response, dataf)
      colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      dataf_original <- data.frame(aggregate_monthly_window(
        env_data_original,
        as.numeric(plot_column):
          (as.numeric(plot_column) + as.numeric(row_index) - 1)
      ))

      dataf_full_original <- dataf_original

      dataf_full_original[, 1] <- dc_apply(dataf_full_original[, 1],
                                           years = row.names(dataf_full_original))

      colnames(dataf_full_original) <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      # Additional check: (we should get the same metric as before in the loop)
      if (method == "lm" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$r.squared
      }

      if (method == "lm" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$adj.r.squared
      }

      if (method == "brnn" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))

        brnn_result <- safe_brnn_metrics(
          temporal_df = temporal_df,
          formula = Optimized.rowNames ~ .,
          observed = temporal_df[, 1],
          neurons = neurons
        )

        optimized_result <- brnn_result$r.squared
      }

      if (method == "brnn" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))

        brnn_result <- safe_brnn_metrics(
          temporal_df = temporal_df,
          formula = Optimized.rowNames ~ .,
          observed = temporal_df[, 1],
          neurons = neurons
        )

        optimized_result <- brnn_result$adj.r.squared
      }

      if (method == "cor"){
        optimized_result <- cor(dataf, response, method = cor_method, use = cor_na_use)
      }

      # Just give a nicer colname
      colnames(dataf) <- "Optimized return"

    }

    # Option 2, reference window = "end"
    if (reference_window == 'end'){

      dataf <- data.frame(aggregate_monthly_window(
        env_data,
        (as.numeric(plot_column) - as.numeric(row_index) + 1):
          as.numeric(plot_column)
      ))

      # if detrending was applied, should also be applied here
      dataf[, 1] <- dc_apply(dataf[, 1], years = row.names(dataf))

      dataf_full <- cbind(response, dataf)
      colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      dataf_original <- data.frame(aggregate_monthly_window(
        env_data_original,
        (as.numeric(plot_column) - as.numeric(row_index) + 1):
          as.numeric(plot_column)
      ))

      dataf_full_original <- dataf_original

      dataf_full_original[, 1] <- dc_apply(dataf_full_original[, 1],
                                           years = row.names(dataf_full_original))

      colnames(dataf_full_original) <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      # Additional check: (we should get the same metric as before in the loop)
      if (method == "lm" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$r.squared
      }

      if (method == "lm" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))
        temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
        temporal_summary <- summary(temporal_model)
        optimized_result <- temporal_summary$adj.r.squared
      }

      if (method == "brnn" & metric == "r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))

        brnn_result <- safe_brnn_metrics(
          temporal_df = temporal_df,
          formula = Optimized.rowNames ~ .,
          observed = temporal_df[, 1],
          neurons = neurons
        )

        optimized_result <- brnn_result$r.squared
      }

      if (method == "brnn" & metric == "adj.r.squared"){
        temporal_df <- data.frame(cbind(dataf, response))

        brnn_result <- safe_brnn_metrics(
          temporal_df = temporal_df,
          formula = Optimized.rowNames ~ .,
          observed = temporal_df[, 1],
          neurons = neurons
        )

        optimized_result <- brnn_result$adj.r.squared
      }

      if (method == "cor"){
        optimized_result <- cor(dataf, response, method = cor_method, use = cor_na_use)
      }

      # Just give a nicer colname
      colnames(dataf) <- "Optimized return"

    }

    # 1. reference window = "middle"
    if (reference_window == 'middle'){

      if (as.numeric(row_index)%%2 == 0){
        adjustment_1 = 0
        adjustment_2 = 1
      } else {
        adjustment_1 = 1
        adjustment_2 = 2
      }

      dataf <- data.frame(aggregate_monthly_window(
        env_data,
        (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
          (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2)
      ))

      # if detrending was applied, should also be applied here
      dataf[, 1] <- dc_apply(dataf[, 1], years = row.names(dataf))

      dataf_full <- cbind(response, dataf)
      colnames(dataf_full)[ncol(dataf_full)] <- "Optimized_return"
      colnames(dataf) <- "Optimized.rowNames"

      ## Once again, the same procedure, to get the optimal sequence, but this time for whole data, not only
      # for the analysed period.

      dataf_original <- data.frame(aggregate_monthly_window(
        env_data_original,
        (round2((as.numeric(plot_column) - as.numeric(row_index) / 2)) - adjustment_1):
          (round2((as.numeric(plot_column) + as.numeric(row_index) / 2)) - adjustment_2)
      ))

      dataf_full_original <- dataf_original

      dataf_full_original[, 1] <- dc_apply(dataf_full_original[, 1],
                                           years = row.names(dataf_full_original))

    }

    colnames(dataf_full_original) <- "Optimized_return"
    colnames(dataf) <- "Optimized.rowNames"



    # Additional check: (we should get the same metric as before in the loop)
    if (method == "lm" & metric == "r.squared"){
      temporal_df <- data.frame(cbind(dataf, response))
      temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
      temporal_summary <- summary(temporal_model)
      optimized_result <- temporal_summary$r.squared
    }

    if (method == "lm" & metric == "adj.r.squared"){
      temporal_df <- data.frame(cbind(dataf, response))
      temporal_model <- lm(Optimized.rowNames ~ ., data = temporal_df)
      temporal_summary <- summary(temporal_model)
      optimized_result <- temporal_summary$adj.r.squared
    }

    if (method == "brnn" & metric == "r.squared"){
      temporal_df <- data.frame(cbind(dataf, response))

      brnn_result <- safe_brnn_metrics(
        temporal_df = temporal_df,
        formula = Optimized.rowNames ~ .,
        observed = temporal_df[, 1],
        neurons = neurons
      )

      optimized_result <- brnn_result$r.squared
    }

    if (method == "brnn" & metric == "adj.r.squared"){
      temporal_df <- data.frame(cbind(dataf, response))

      brnn_result <- safe_brnn_metrics(
        temporal_df = temporal_df,
        formula = Optimized.rowNames ~ .,
        observed = temporal_df[, 1],
        neurons = neurons
      )

      optimized_result <- brnn_result$adj.r.squared
    }

    if (method == "cor"){
      optimized_result <- cor(dataf, response, method = cor_method, use = cor_na_use)
    }

    # Just give a nicer colname
    colnames(dataf) <- "Optimized return"

  }


  # Element 5
  # Here we create the fifth element of the final list: Analysed period in the
  # form of min(year) - max(year), e.g. 1950 - 2015
  min_env_data <- min(as.numeric(row.names(env_data)))
  min_response <- min(as.numeric(row.names(response)))

  max_env_data <- max(as.numeric(row.names(env_data)))
  max_response <- max(as.numeric(row.names(response)))

  min_together <- min(min_env_data, min_response)
  max_together <- min(max_env_data, max_response)


  analysed_period <- paste(as.character(min_together),
                           as.character(max_together),
                           sep = " - ")
  if (nchar(analysed_period) < 9) {
    analysed_period <- NA
  }

  # Here, the transfer function is being created
  transfer_data = data.frame(proxy = response[,1], optimized_return =dataf[,1])
  lm_model = lm(optimized_return ~ proxy, data = transfer_data)

  brnn_transfer_result <- safe_brnn_metrics(
    temporal_df = transfer_data,
    formula = optimized_return ~ proxy,
    observed = transfer_data$optimized_return,
    neurons = neurons
  )

  brnn_model <- brnn_transfer_result$model

  full_range = data.frame(proxy = seq(from = min(response[,1], na.rm = TRUE), to = max(response[,1], na.rm = TRUE), length.out = 100))

  if (method == "lm" | method == "cor"){
    full_range$transfer_f = predict(lm_model, full_range)
  }

  if (method == "brnn"){
    if (!is.null(brnn_model)) {
      transfer_predictions <- try(as.numeric(predict.brnn(brnn_model, full_range)),
                                  silent = TRUE)

      if (inherits(transfer_predictions, "try-error") ||
          length(transfer_predictions) != nrow(full_range) ||
          any(!is.finite(transfer_predictions))) {
        full_range$transfer_f <- NA_real_
      } else {
        full_range$transfer_f <- transfer_predictions
      }
    } else {
      full_range$transfer_f <- NA_real_
    }
  }

  # String for titles
  if (method == "cor"){
    title_string <- "Correlation Coefficients"
  } else if (method == "lm"){
    title_string <- "Linear Regression"
  } else if (method == "brnn"){
    title_string <- "ANN With Bayesian Regularization"
  } else (print("The selection of method is not correct"))

  # The definition of theme
  journal_theme <- theme_bw() +
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 18), text = element_text(size = 18),
          plot.title = element_text(size = 16,  face = "bold"))

  p1 <- ggplot(transfer_data, aes(proxy, optimized_return)) +
    geom_point() +
    geom_line(aes(proxy, transfer_f), full_range) +
    journal_theme +
    ggtitle(paste("Analysed Period:", analysed_period, "\nMethod:", title_string))

  # If there is more than one independent variable in the model,
  # transfer function is not given, since we should return a 3d model
  if (ncol(response) > 1){
    p1 <- "No transfer function is created for two or more response variables"
  }


  #######################################################################
  ############## The temporal stability of optimized_return #############
  #######################################################################

  dataset = data.frame(optimized_return =dataf[,1], proxy = response)

  empty_list = list()
  empty_list_period = list()
  empty_list_significance = list()

  temporal_stability <- data.frame()

  # 1. Progressive stability check
  if (temporal_stability_check == "progressive"){
    foldi <- seq(1:k)
    #foldi <- paste("fold_", foldi)
    folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

    for (m in 1:k){
      #Segement your data by fold using the which() function
      trainIndexes <- which(folds <= m, arr.ind = TRUE)
      dataset_temp <- dataset[trainIndexes, ]
      MAKS <- max(as.numeric(row.names(dataset_temp)))
      MIN <- min(as.numeric(row.names(dataset_temp)))
      empty_list_period[[m]] <- paste(MIN, "-", MAKS)

      if (method == "cor"){
        calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method, use = cor_na_use)
        sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exact = FALSE)$p.value
        empty_list_significance[[m]] <- sig
        empty_list[[m]] <- calculation
        colname = "correlation"
      } else if (method == "lm" & metric == "r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        colname = "r.squared"
        empty_list[[m]] <- summary(MLR)$r.squared
        empty_list_significance[[m]] <- NA
      } else if (method == "lm" & metric == "adj.r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        empty_list[[m]] <- summary(MLR)$adj.r.squared
        empty_list_significance[[m]] <- NA
        colname = "adj.r.squared"
      } else if (method == "brnn" & metric == "r.squared"){

        brnn_result <- safe_brnn_metrics(
          temporal_df = dataset_temp,
          formula = optimized_return ~ .,
          observed = dataset_temp[, 1],
          neurons = neurons
        )

        empty_list[[m]] <- brnn_result$r.squared
        empty_list_significance[[m]] <- NA
        colname = "r.squared"

      } else if (method == "brnn" & metric == "adj.r.squared"){

        brnn_result <- safe_brnn_metrics(
          temporal_df = dataset_temp,
          formula = optimized_return ~ .,
          observed = dataset_temp[, 1],
          neurons = neurons
        )

        empty_list[[m]] <- brnn_result$adj.r.squared
        empty_list_significance[[m]] <- NA
        colname = "adj.r.squared"
      }
    }
    m1 <- do.call(rbind, empty_list)
    m2 <- do.call(rbind, empty_list_period)
    m3 <- do.call(rbind, empty_list_significance)

    temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(m3, 4), nsmall = 3)))
    colnames(temporal_stability) <-c("Period", colname, "p value")
    temporal_stability
  }

  # 2. Sequential stability check
  if (temporal_stability_check == "sequential"){
    foldi <- seq(1:k)
    #foldi <- paste("fold_", foldi)
    folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)

    for (m in 1:k){
      #Segement your data by fold using the which() function
      trainIndexes <- which(folds == m, arr.ind = TRUE)
      dataset_temp <- dataset[trainIndexes, ]

      MAKS <- max(as.numeric(row.names(dataset_temp)))
      MIN <- min(as.numeric(row.names(dataset_temp)))
      empty_list_period[[m]] <- paste(MIN, "-", MAKS)

      if (method == "cor"){
        calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method, use = cor_na_use)
        sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exact = FALSE)$p.value
        empty_list[[m]] <- calculation
        empty_list_significance[[m]] <- sig
        colname = "correlation"

      } else if (method == "lm" & metric == "r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        colname = "r.squared"
        empty_list[[m]] <- summary(MLR)$r.squared
        empty_list_significance[[m]] <- NA
      } else if (method == "lm" & metric == "adj.r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        empty_list[[m]] <- summary(MLR)$adj.r.squared
        empty_list_significance[[m]] <- NA
        colname = "adj.r.squared"
      } else if (method == "brnn" & metric == "r.squared"){

        brnn_result <- safe_brnn_metrics(
          temporal_df = dataset_temp,
          formula = optimized_return ~ .,
          observed = dataset_temp[, 1],
          neurons = neurons
        )

        empty_list[[m]] <- brnn_result$r.squared
        empty_list_significance[[m]] <- NA
        colname = "r.squared"

      } else if (method == "brnn" & metric == "adj.r.squared"){

        brnn_result <- safe_brnn_metrics(
          temporal_df = dataset_temp,
          formula = optimized_return ~ .,
          observed = dataset_temp[, 1],
          neurons = neurons
        )

        empty_list[[m]] <- brnn_result$adj.r.squared
        empty_list_significance[[m]] <- NA
        colname = "adj.r.squared"
      }
    }
    m1 <- do.call(rbind, empty_list)
    m2 <- do.call(rbind, empty_list_period)
    m3 <- do.call(rbind, empty_list_significance)

    temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(m3, 4), nsmall = 3)))
    colnames(temporal_stability) <-c("Period", colname, "p value")
    temporal_stability

  }


  # 3. running_window
  if (temporal_stability_check == "running_window"){

    k_end <- nrow(dataset)
    place_list <- 1
    empty_list_datasets <- list()
    empty_list_period <- list()

    empty_list <- list()
    empty_list_significance <- list()

    place_list = 1

    if (k_end < k_running_window){
      stop("k_running_window is less than the number of analysed years. Reduce the argument k_running_window")
    }

    for (w in 0:(k_end - k_running_window)){

      dataset_temp <- dataset[(1+w):(k_running_window + w),]
      empty_list_datasets[[place_list]] <- dataset_temp

      MAKS <- max(as.numeric(row.names(dataset_temp)))
      MIN <- min(as.numeric(row.names(dataset_temp)))
      empty_list_period[[place_list]] <- paste(MIN, "-", MAKS)

      place_list <- place_list + 1


    }

    for (m in 1:length(empty_list_datasets)){

      dataset_temp <- empty_list_datasets[[m]]

      if (method == "cor"){
        calculation <- cor(dataset_temp[,1], dataset_temp[,2], method = cor_method, use = cor_na_use)
        sig <- cor.test(dataset_temp[,1], dataset_temp[,2], method = cor_method, exact = FALSE)$p.value
        empty_list[[m]] <- calculation
        empty_list_significance[[m]] <- sig
        colname = "correlation"

      } else if (method == "lm" & metric == "r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        colname = "r.squared"
        empty_list[[m]] <- summary(MLR)$r.squared
        empty_list_significance[[m]] <- NA
      } else if (method == "lm" & metric == "adj.r.squared"){
        MLR <- lm(optimized_return ~ ., data = dataset_temp)
        empty_list[[m]] <- summary(MLR)$adj.r.squared
        empty_list_significance[[m]] <- NA
        colname = "adj.r.squared"
      } else if (method == "brnn" & metric == "r.squared"){

        brnn_result <- safe_brnn_metrics(
          temporal_df = dataset_temp,
          formula = optimized_return ~ .,
          observed = dataset_temp[, 1],
          neurons = neurons
        )

        empty_list[[m]] <- brnn_result$r.squared
        empty_list_significance[[m]] <- NA
        colname = "r.squared"

      } else if (method == "brnn" & metric == "adj.r.squared"){

        brnn_result <- safe_brnn_metrics(
          temporal_df = dataset_temp,
          formula = optimized_return ~ .,
          observed = dataset_temp[, 1],
          neurons = neurons
        )

        empty_list[[m]] <- brnn_result$adj.r.squared
        empty_list_significance[[m]] <- NA
        colname = "adj.r.squared"
      }
    }
    m1 <- do.call(rbind, empty_list)
    m2 <- do.call(rbind, empty_list_period)
    m3 <- do.call(rbind, empty_list_significance)

    temporal_stability <- data.frame(cbind(m2, format(round(m1, 3), nsmall = 3), format(round(as.numeric(m3), digits = 3), nsmall = 3)))
    colnames(temporal_stability) <-c("Period", colname, "p value")
    temporal_stability
  }




  #########################################################################
  ################## Out of sample estimates ##############################
  #########################################################################

  dataset = data.frame(optimized_return =dataf[,1], proxy = response)

  # Years without a value, such as the first years of a window prewhitened
  # with dc_method = "Ar", can neither be predicted nor used for calibration,
  # and they made every fold that contained them return NA. They are left
  # out before the folds are formed. Complete data are not affected.
  dataset <- dataset[stats::complete.cases(dataset), , drop = FALSE]

  empty_list = list()
  empty_list_period = list()

  if (cross_validation_type == "blocked"){
    dataset <- dataset
  } else if (cross_validation_type == "randomized"){
    dataset <- dataset[sample(nrow(dataset)), ]
  } else (stop(paste("The cross_validation_type is not selected correctly! It is ", cross_validation_type,
                     ". It should be 'blocked' or 'randomized'!", sep = "")))

  foldi <- seq(1:k)
  #foldi <- paste("fold_", foldi)
  folds <- cut(seq(1, nrow(dataset)), breaks = k, labels = FALSE)
  bl = 1

  for (m in 1:k){
    #Segement your data by fold using the which() function
    testIndexes <- which(folds == m, arr.ind = TRUE)
    test <- dataset[testIndexes, ]
    train <- dataset[-testIndexes, ]

    MAKS <- max(as.numeric(row.names(train)))
    MIN <- min(as.numeric(row.names(train)))
    empty_list_period[[bl]] <- paste(MIN, "-", MAKS)
    bl <- bl + 1

    MAKS <- max(as.numeric(row.names(test)))
    MIN <- min(as.numeric(row.names(test)))
    empty_list_period[[bl]] <- paste(MIN, "-", MAKS)
    bl <- bl + 1

    if (method == "lm" | method == "cor"){
      MLR <- lm(optimized_return ~ ., data = train)
      train_predicted <- predict(MLR, train)
      test_predicted <- predict(MLR, test)
      train_observed <- train[, 1]
      test_observed <- test[, 1]
      calculations <- calculate_metrics(train_predicted, test_predicted,
                                        train_observed, test_observed, test = test,
                                        formula = optimized_return ~ ., digits = 15)

      empty_list[[m]] <- calculations
    }

    if (method == "brnn"){

      brnn_cv_result <- safe_brnn_metrics(
        temporal_df = train,
        formula = optimized_return ~ .,
        observed = train[, 1],
        neurons = neurons
      )

      BRNN <- brnn_cv_result$model

      if (!is.null(BRNN)){

        train_predicted <- try(as.numeric(predict(BRNN, train)), silent = TRUE)
        test_predicted <- try(as.numeric(predict(BRNN, test)), silent = TRUE)

        if (!inherits(train_predicted, "try-error") &&
            !inherits(test_predicted, "try-error") &&
            length(train_predicted) == nrow(train) &&
            length(test_predicted) == nrow(test) &&
            all(is.finite(train_predicted)) &&
            all(is.finite(test_predicted))) {

          train_observed <- train[, 1]
          test_observed <- test[, 1]

          calculations <- calculate_metrics(train_predicted, test_predicted,
                                            train_observed, test_observed,
                                            digits = 15,
                                            test = test,
                                            formula = optimized_return ~ .)

          empty_list[[m]] <- calculations

        } else {

          empty_list[[m]] <- empty_cv_metrics()
        }

      } else {

        empty_list[[m]] <- empty_cv_metrics()
      }
    }


  }

  empty_list <- normalize_cv_metric_list(empty_list)

  m1 <- do.call(rbind, empty_list)
  # m1 <- m1[, -c(3, 4, 7)]
  m2 <- do.call(rbind, empty_list_period)

  if (nrow(as.data.frame(m1)) != nrow(as.data.frame(m2))) {
    stop(paste0("Internal cross-validation error: metric table has ",
                nrow(as.data.frame(m1)), " rows but period table has ",
                nrow(as.data.frame(m2)), " rows. This usually indicates ",
                "that a model-fitting branch did not return both calibration ",
                "and validation rows."))
  }

  cross_validation <- cbind(Years = m2, m1)
  cross_validation$Period <- c("Calibration", "Validation")
  cross_validation$CV <- rep(1:k, each = 2)
  row.names(cross_validation) <- NULL

  if (cross_validation_type == "blocked"){
    cross_validation <- dplyr::select(cross_validation, CV, Period, Years, cor, RMSE, RRSE, d, RE, CE, DE)
  }
  if (cross_validation_type == "randomized"){
    cross_validation <- dplyr::select(cross_validation, CV, Period, cor, RMSE, RRSE, d, RE, CE, DE)
  }

  ################################################################
  #### Here the final list is being filled with six elements #####
  ################################################################

  # When metohod == "cor", different final_list is created
  if (method == "lm" | method == "brnn") {
    final_list <- list(calculations = temporal_matrix, method = method,
                       metric = metric, analysed_period = analysed_period,
                       optimized_return = dataf_full,
                       optimized_return_all = dataf_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = cross_validation,
                       aggregated_climate = do.call(cbind, list_climate),
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       significance_filter = significance_filter,
                       detrending = detrending_summary)
  }

  if (method == "cor"){
    final_list <- list(calculations = temporal_matrix, method = method,
                       metric = cor_method, analysed_period = analysed_period,
                       optimized_return = dataf_full,
                       optimized_return_all = dataf_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = cross_validation,
                       aggregated_climate = do.call(cbind, list_climate),
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       significance_filter = significance_filter,
                       detrending = detrending_summary)
  }



  plot_heatmapA <- plot_heatmap(final_list, reference_window = reference_window, type = "monthly")
  plot_extremeA <- plot_extreme(final_list, ylimits = ylimits, reference_window = reference_window, type = "monthly")

  # Here, for the sake of simplicity, we create final list again
  if (method == "lm" | method == "brnn") {
    final_list <- list(calculations = temporal_matrix, method = method,
                       metric = metric, analysed_period = analysed_period,
                       optimized_return = dataf_full,
                       optimized_return_all = dataf_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = cross_validation,
                       plot_heatmap = plot_heatmapA,
                       plot_extreme = plot_extremeA,
                       type = "monthly",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate = do.call(cbind, list_climate),
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       significance_filter = significance_filter,
                       detrending = detrending_summary)
  }

  if (method == "cor"){

    final_list <- list(calculations = temporal_matrix, method = method,
                       metric = cor_method, analysed_period = analysed_period,
                       optimized_return = dataf_full,
                       optimized_return_all = dataf_full_original,
                       transfer_function = p1, temporal_stability = temporal_stability,
                       cross_validation = cross_validation,
                       plot_heatmap = plot_heatmapA,
                       plot_extreme = plot_extremeA,
                       type = "monthly",
                       reference_window = reference_window,
                       boot_lower = temporal_matrix_lower,
                       boot_upper = temporal_matrix_upper,
                       aggregated_climate = do.call(cbind, list_climate),
                       previous_year = previous_year,
                       number_previous_years = number_previous_years,
                       significance_filter = significance_filter,
                       detrending = detrending_summary)
  }

  class(final_list) <- 'dmrs'

  return(final_list)

}
