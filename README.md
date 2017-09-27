
<!-- README.md is generated from README.Rmd. Please edit that file -->
dendroTools
===========

The core purpose of the dendroTools package is to introduce novel dendroclimatological methods to study linear and nonlinear relationship between daily climate data and tree-ring sequences. There are two core functions. The first core function is daily\_response(), which finds the optimal sequence of days that are linearly or nonlinearly related to a tree-ring proxy records. The second core function is compare\_methods() that calculates several performance measures for train and test data for different regression method.

To use daily\_response(), two data frames are required, one with daily climate data, e.g. temperatures; and one with tree-ring proxy records. Example data is provided, so users can see, how data frames should be organized. The daily\_response() calculates all possible values of a selected statistical measure between response variables and daily environmental data. Calculations are based on a moving window, which runs through daily environmental data and calculates moving averages.

compare\_methods() calcualtes performance measures for multiple linear regression (MLR), artificial neural networks with Bayesian regularization training algorithm (ANN), M5P model trees (MT), model trees with bagging (BMT) and random forest of regression trees (RF). Calculated performance measures are correlation coefficient, root mean squared error (RMSE), root relative squared error (RSSE), index of agreement (d), reduction of error (RE), coefficient of efficiency (CE) nad mean bias.

Installation
------------

You can install dendroTools from github with:

``` r
# install.packages("devtools")
devtools::install_github("jernejjevsenak/dendroTools")
```

Examples
--------

This is a basic example which shows you how to use the daily\_response():

``` r
## basic example code
library(dendroTools)
data(daily_temperatures_example) 
data(example_proxies_1)
result1 <- daily_response(response = example_proxies_1, env_data = daily_temperatures_example, 
                            method = "lm", measure = "r.squared", lower_limit = 90, upper_limit = 150)
```

This function is computationally intensive and it takes a while to calculate all possible values. Especially, if nonlinear "brnn" method is used. Each calculated value is printed, therefore user can be sure, that algorithm is still running. The return of this function is a list with three elements: @calculations, @method, @measure. The return is organized in a way, that can be used by three plotting functions: plot\_extreme(), plot\_specific() and plot\_heatmap(). Function plot\_extreme() graphs a line plot of a row with the highest calculated measure. It indicates the sequence of days, that give the best value of selected statistical measure. With plot\_specific(), measures with selected window width are plotted. Function plot\_heatmap() is a visual representation of calculated values.

``` r
plot_extreme(result1, title = TRUE)
```

![](README-plot%20examples-1.png)

``` r
plot_specific(result1, window_width = 100, title = TRUE)
```

![](README-plot%20examples-2.png)

``` r
plot_heatmap(result1)
```

![](README-plot%20examples-3.png)

Here is an example of using the compare\_methods(). The core idea is to compare different regression methods and select the best model for climate reconstruction. Users should spend some time to understand what are the tunning parameters of specific method. Basically, methods could be run with default parameters, but usually it turns out that optimization of parameters gives better results. To optimize parameters, users might want to explore the caret package.

``` r
## basic example with default parameter values
library(dendroTools)
data(example_dataset_1)
experiment_1 <- compare_methods(formula = MVA~.,
dataset = example_dataset_1, k = 3, repeats = 100)
```

To see the outputs, use the following code.

``` r
experiment_1[[1]] # See a data frame results
#>    Measure Period MLR_M MLR_SD MLR_AR MLR_S1 ANN_M ANN_SD ANN_AR ANN_S1
#> 5        r    cal 0.676  0.006   3.88   0.00 0.677  0.006   3.25   0.00
#> 6        r    val 0.642  0.030   1.48   0.55 0.642  0.030   2.13   0.45
#> 9     RMSE    cal 0.389  0.003   3.61   0.00 0.388  0.003   3.58   0.00
#> 10    RMSE    val 0.419  0.014   1.49   0.54 0.419  0.014   2.16   0.45
#> 11    RSSE    cal 0.734  0.005   3.62   0.00 0.734  0.005   3.68   0.00
#> 12    RSSE    val 0.813  0.044   1.54   0.50 0.813  0.044   2.18   0.46
#> 3        d    cal 0.782  0.006   2.66   0.00 0.775  0.006   4.95   0.00
#> 4        d    val 0.740  0.022   1.04   0.96 0.732  0.022   3.10   0.03
#> 7       RE    cal 0.490  0.029   3.67   0.00 0.491  0.029   3.66   0.00
#> 8       RE    val 0.372  0.070   1.63   0.46 0.372  0.070   2.09   0.50
#> 1       CE    cal 0.459  0.008   3.66   0.00 0.460  0.008   3.64   0.00
#> 2       CE    val 0.326  0.084   1.61   0.47 0.326  0.085   2.12   0.48
#>     MT_M MT_SD MT_AR MT_S1 BMT_M BMT_SD BMT_AR BMT_S1  RF_M RF_SD RF_AR
#> 5  0.676 0.006  3.98  0.00 0.687  0.007   2.01   0.00 0.959 0.003  1.00
#> 6  0.636 0.041  1.84  0.51 0.631  0.035   3.71   0.01 0.534 0.057  4.98
#> 9  0.389 0.003  3.77  0.00 0.384  0.003   2.01   0.00 0.171 0.005  1.00
#> 10 0.421 0.018  1.86  0.49 0.424  0.015   3.63   0.02 0.469 0.023  4.99
#> 11 0.734 0.006  3.78  0.00 0.725  0.006   2.01   0.00 0.323 0.009  1.00
#> 12 0.818 0.047  1.90  0.45 0.822  0.044   3.56   0.04 0.912 0.058  4.98
#> 3  0.782 0.006  2.80  0.00 0.784  0.007   2.71   0.00 0.967 0.002  1.00
#> 4  0.735 0.030  1.50  0.82 0.727  0.025   3.55   0.01 0.673 0.038  4.95
#> 7  0.490 0.029  3.82  0.00 0.503  0.029   2.01   0.00 0.901 0.008  1.00
#> 8  0.364 0.075  1.98  0.42 0.358  0.071   3.49   0.04 0.208 0.101  4.95
#> 1  0.459 0.008  3.81  0.00 0.473  0.009   2.01   0.00 0.895 0.006  1.00
#> 2  0.318 0.089  1.96  0.43 0.311  0.085   3.52   0.04 0.150 0.124  4.95
#>    RF_S1
#> 5   1.00
#> 6   0.00
#> 9   1.00
#> 10  0.00
#> 11  1.00
#> 12  0.00
#> 3   1.00
#> 4   0.00
#> 7   1.00
#> 8   0.01
#> 1   1.00
#> 2   0.01
experiment_1[[2]] # See a ggplot of mean bias for validation data
```

![](README-example%202a-1.png)

``` r
## An example with tunned parameter values
library(dendroTools)
data(example_dataset_1)
experiment_2 <- compare_methods(formula = MVA~.,
dataset = example_dataset_1, k = 3, repeats = 100, neurons = 1,
MT_M = 4, MT_N = FALSE, MT_U = FALSE, MT_R = FALSE, BMT_P = 100,
BMT_I = 100, BMT_M = 4, BMT_N = FALSE, BMT_U = FALSE, BMT_R = FALSE,
RF_P = 100, RF_I = 100, RF_depth= 0, multiply = 5)
```

Aagin, to see the outputs, use the following code.

``` r
experiment_2[[1]] # See a data frame results
#>    Measure Period MLR_M MLR_SD MLR_AR MLR_S1 ANN_M ANN_SD ANN_AR ANN_S1
#> 5        r    cal 0.676  0.006   3.88   0.00 0.677  0.006   3.25   0.00
#> 6        r    val 0.642  0.030   1.48   0.55 0.642  0.030   2.13   0.45
#> 9     RMSE    cal 0.389  0.003   3.61   0.00 0.388  0.003   3.58   0.00
#> 10    RMSE    val 0.419  0.014   1.49   0.54 0.419  0.014   2.16   0.45
#> 11    RSSE    cal 0.734  0.005   3.62   0.00 0.734  0.005   3.68   0.00
#> 12    RSSE    val 0.813  0.044   1.54   0.50 0.813  0.044   2.18   0.46
#> 3        d    cal 0.782  0.006   2.66   0.00 0.775  0.006   4.95   0.00
#> 4        d    val 0.740  0.022   1.04   0.96 0.732  0.022   3.10   0.03
#> 7       RE    cal 0.490  0.029   3.67   0.00 0.491  0.029   3.66   0.00
#> 8       RE    val 0.372  0.070   1.63   0.46 0.372  0.070   2.09   0.50
#> 1       CE    cal 0.459  0.008   3.66   0.00 0.460  0.008   3.64   0.00
#> 2       CE    val 0.326  0.084   1.61   0.47 0.326  0.085   2.12   0.48
#>     MT_M MT_SD MT_AR MT_S1 BMT_M BMT_SD BMT_AR BMT_S1  RF_M RF_SD RF_AR
#> 5  0.676 0.006  3.98  0.00 0.687  0.007   2.01   0.00 0.959 0.003  1.00
#> 6  0.636 0.041  1.84  0.51 0.631  0.035   3.71   0.01 0.534 0.057  4.98
#> 9  0.389 0.003  3.77  0.00 0.384  0.003   2.01   0.00 0.171 0.005  1.00
#> 10 0.421 0.018  1.86  0.49 0.424  0.015   3.63   0.02 0.469 0.023  4.99
#> 11 0.734 0.006  3.78  0.00 0.725  0.006   2.01   0.00 0.323 0.009  1.00
#> 12 0.818 0.047  1.90  0.45 0.822  0.044   3.56   0.04 0.912 0.058  4.98
#> 3  0.782 0.006  2.80  0.00 0.784  0.007   2.71   0.00 0.967 0.002  1.00
#> 4  0.735 0.030  1.50  0.82 0.727  0.025   3.55   0.01 0.673 0.038  4.95
#> 7  0.490 0.029  3.82  0.00 0.503  0.029   2.01   0.00 0.901 0.008  1.00
#> 8  0.364 0.075  1.98  0.42 0.358  0.071   3.49   0.04 0.208 0.101  4.95
#> 1  0.459 0.008  3.81  0.00 0.473  0.009   2.01   0.00 0.895 0.006  1.00
#> 2  0.318 0.089  1.96  0.43 0.311  0.085   3.52   0.04 0.150 0.124  4.95
#>    RF_S1
#> 5   1.00
#> 6   0.00
#> 9   1.00
#> 10  0.00
#> 11  1.00
#> 12  0.00
#> 3   1.00
#> 4   0.00
#> 7   1.00
#> 8   0.01
#> 1   1.00
#> 2   0.01
experiment_2[[2]] # See a ggplot of mean bias for validation data
```

![](README-example%203a-1.png)

Authors
-------

-   **Jernej Jevšenak**

Collaborators
-------------

-   **Prof. dr. Tom Levanič**
