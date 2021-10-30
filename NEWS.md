# dendroTools 1.2.1
Minor correction for ggplot2 object -> guide = FALSE --> guide = 'none'


# dendroTools 1.2.0
There is a major upgrade related to new arguments that now allow you to specify the time of interest to be used for calculating correlation coefficients or other statistical metrics.
For example, if you are using daily_response() or daily_response_seascorr() and want to calculate correlations only from the beginning of last October to the current end of September, use the argument daily_interval =c(-274, 273). Note that a negative sign indicates the previous year's doy, while a positive sign indicates the current year doy. In normal (non-leap ) years, October 1 represents doy 274, and September 30 represents doy 273.
If you use monthly_response() or monthly_response_seascorr() and want to calculate correlations only from the beginning of the previous October to the current end of September, for example, use the argument monthly_interval =c(-10, 9). Note that a negative sign indicates the month of the previous year, while a positive sign indicates the month of the current year.
Other functions have been changed only slightly, so they will work even if new arguments are used.


# dendroTools 1.1.4
* new error is implemented in case of wrong subset_years
* minor corrections to partial correlation functions

# dendroTools 1.1.3
* deleted example in vignette Examples_with_daily_climate_data which caused problems with CRAN check

# dendroTools 1.1.2
* a bug corrected in the data_transform()
* minor bug is removed which happened in a case of nonsignificant calculations or only one significant correlation

# dendroTools 1.1.1
* dendroTools version is changed to 1.1.1
* Correctcions realted to testhat changes.
* imports is updated: oce(>= 1.2-0)

# dendroTools 1.0.8
* Minor corrections after new publication
* Sys.setlocale("LC_TIME", "C") is now added, so optimal dates will always show in English language
* The vignette titles are changed

# dendroTools 1.0.7.
* Bootstrapping is now available for all metrics, monthly and daily functions
* Seacorr is now available also for monthly data, including bootstrapping
* swit272 daily precipitation dataset is added
* swit272 monthly temperature dataset is removed
* KNMI_daily_transform() is now called data_transform(). The entire function is updated and now also enables transformation of daily data into monthly.

# dendroTools 1.0.6.
* Bootstrapping of correlation coefficients is introduced. To do so, use the argument boot in daily_response() and monthly_response(). Bootstrapping is currently available only for correlation coefficients. In future version, it will also be available for model fitting. 
* Package version is changed to 1.0.6

# dendroTools 1.0.5.
* new function: monthly_response_seascorr()
* titles of plot_extreme, plot_specific and plot_heatmap are updated
* for daily_response() and monthly_response(), it is now possible to define method for correlation coefficient: "pearson", "kendall", "spearman". To do so, use , cor_method argument
* daily_response output list is class "dmrs"
* Using newly defined class "dmrs", new summary function is defined
* Package version is changed to 1.0.5

# dendroTools 1.0.4.
* New data is available for examples demonstration: swit272_monthly_temperatures
* New function is available: daily_response_seascorr(). This function implements partial correlations in the analysis of daily response functions
* Package version is changed to 1.0.4

# dendroTools 1.0.3.
* Package version is changed to 1.0.3
* There is new function introduced: KNMI_daily_transform() which transforms daily data obtained from KNMI Climate explorer into data frame suitable for daily_response(). 
* New TRW chronology is included, swit272
* New function: monthly_response()
* New data: LJ_monthly_temperatures
* New data: LJ_monthly_precipitation

# dendroTools 1.0.2.
* Package version is changed to 1.0.2
* There have been many issues with RWeka package, which depends on rJava. Therefore, functions from RWeka are replaced with other functions. For random forest model we now use randomForest package
* There are six new output elements, all of them are residual diagnostic plots for calibration, holdout and edge data.
* Due to many problems related to rJava in R, MT and BMT methods from RWeka are now joined and replaced by cubist function from Cubist r package. Cubist fits a model tree or ensemble of model trees, if committees argument is more than 1.
