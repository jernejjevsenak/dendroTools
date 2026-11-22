
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dendroTools

dendroTools provides dendroclimatological methods for studying linear
and nonlinear relationships between tree-ring proxies and daily or
monthly climate data in R.

## Main functions

- *daily_response()* and *monthly_response()* aggregate daily or monthly
  climate data over moving windows of all widths and positions,
  optionally extending into up to five previous years. For each window,
  they calculate the correlation (Pearson, Spearman or Kendall) with a
  tree-ring proxy, or fit linear models or neural networks with Bayesian
  regularization (brnn) with one or more proxies, and find the optimal
  window, i.e. the window with the strongest relationship.
- *daily_response_seascorr()* and *monthly_response_seascorr()*
  calculate partial correlations between a tree-ring proxy and a primary
  climate variable, controlling for a second variable, following
  seascorr ([Meko et
  al. 2011](https://doi.org/10.1016/j.cageo.2011.01.013)). They
  partition the variance explained in the optimal window between the two
  variables, and with `pcor_control = TRUE` they also return the partial
  correlations of the second variable.
- *water_balance()* calculates the daily or monthly climatic water
  balance, i.e. precipitation minus potential evapotranspiration from
  the Hargreaves equation. With `aggregate_function = "SPEI"`, the
  response functions turn it into the Standardized
  Precipitation-Evapotranspiration Index (SPEI) of every window.
- *compare_methods()* compares multiple linear regression (MLR),
  artificial neural networks with Bayesian regularization (BRNN),
  (ensembles of) model trees (MT) and random forests (RF) for climate
  reconstruction. The performance metrics are the correlation
  coefficient (r), root mean squared error (RMSE), root relative squared
  error (RRSE), index of agreement (d), reduction of error (RE),
  coefficient of efficiency (CE), detrended efficiency (DE) and mean
  bias.
- *data_transform()* turns daily data with a date and a value column
  into the wide format used by the response functions,
  *glimpse_daily_data()* shows daily data and their gaps, and
  *summary()* and *plot()* report and visualise the results.

The aggregated climate data can be detrended before the analysis, by
simple linear detrending or with the spline, mean, autoregressive and
Friedman methods of dplR. Windows that are not significant can be
removed with a per-window test or with corrections for multiple testing
(false discovery rate, Šidák correction, and max-statistic and
cluster-extent permutation tests). Bootstrap confidence intervals are
available in all four response functions.

The package includes example tree-ring and climate data, among them the
raw ring widths and E-OBS daily climate data of black pine from Sežana,
Slovenia (`SEZ_rwl`, `SEZ_daily_tmin`, `SEZ_daily_tavg`,
`SEZ_daily_tmax` and `SEZ_daily_psum`).

## Example

``` r
library(dendroTools)

data(data_MVA)               # mean vessel area chronology, Mlače, Slovenia
data(LJ_daily_temperatures)  # daily mean temperatures, Ljubljana

result <- daily_response(response = data_MVA, env_data = LJ_daily_temperatures,
                         method = "cor", lower_limit = 21, upper_limit = 90,
                         row_names_subset = TRUE, remove_insignificant = TRUE,
                         remove_insignificant_method = "fdr_BH")

summary(result)         # optimal window: 14 March - 12 May, r = 0.767
plot(result, type = 1)  # correlations of all windows with the optimal width
plot(result, type = 2)  # correlations of all windows
```

## Vignettes

- *Examples_with_daily_climate_data*: introduction to *daily_response()*
- *Partial_correlations_seascorr*: partial correlations with
  *daily_response_seascorr()* and *monthly_response_seascorr()*
- *Detrending_climate_data*: detrending the aggregated climate data
- *SPEI_drought_index*: climatic water balance and SPEI, with the data
  from Sežana
- *Compare_different_regression_methods*: *compare_methods()*
- *Modify_plots_with_ggplot*: customising the plots

Open them with `browseVignettes("dendroTools")` or, for example,
`vignette("SPEI_drought_index", package = "dendroTools")`.

## Installation

``` r
install.packages("dendroTools") # from CRAN

# development version from GitHub, with vignettes
# install.packages("devtools")
devtools::install_github("jernejjevsenak/dendroTools", build_vignettes = TRUE)
```

## Citation

Please cite dendroTools as:

Jevšenak J. (2020) New features in the dendroTools R package:
Bootstrapped and partial correlation coefficients for monthly and daily
climate data. Dendrochronologia 63, 125753.
<https://doi.org/10.1016/j.dendro.2020.125753>

The methods are also described in:

- Jevšenak J., Levanič T. (2018) dendroTools: R package for studying
  linear and nonlinear responses between tree-rings and daily
  environmental data. Dendrochronologia 48, 32–39.
  <https://doi.org/10.1016/j.dendro.2018.01.005>
- Jevšenak J., Levanič T., Džeroski S. (2018) Comparison of an optimal
  regression method for climate reconstruction with the
  compare_methods() function from the dendroTools R package.
  Dendrochronologia 52, 96–104.
  <https://doi.org/10.1016/j.dendro.2018.10.001>

## Authors

- **Jernej Jevšenak**

Please report bugs and suggestions at
<https://github.com/jernejjevsenak/dendroTools/issues>.
