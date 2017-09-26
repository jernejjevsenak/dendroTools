
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

Here is an example on using the compare\_methods(). The core idea is to compare different regression methods and select the best model for climate reconstruction. Users should spend some time to understand what are the tunning parameters of specific method. Basically, methods could be run with default parameters, but usually it turns out that optimization of parameters gives better results. To optimize parameters, users might want to explore the caret package.

``` r
## basic example with default parameter values
library(dendroTools)
data(example_dataset_1)
experiment_1 <- compare_methods(formula = MVA~.,
dataset = example_dataset_1, k = 3, repeats = 10)
```

To see the outputs, use the following code.

``` r
experiment_1[[1]] # See a data frame results
#>    Measure Period MLR_M MLR_SD MLR_AR MLR_S1 ANN_M ANN_SD ANN_AR ANN_S1
#> 5        r    cal 0.675  0.005    4.0    0.0 0.676  0.005    3.3    0.0
#> 6        r    val 0.646  0.024    1.4    0.6 0.646  0.023    2.1    0.4
#> 9     RMSE    cal 0.389  0.002    3.6    0.0 0.389  0.002    3.5    0.0
#> 10    RMSE    val 0.417  0.009    1.5    0.5 0.417  0.009    1.9    0.5
#> 11    RSSE    cal 0.734  0.004    3.5    0.0 0.734  0.004    3.8    0.0
#> 12    RSSE    val 0.807  0.023    1.4    0.6 0.806  0.024    2.0    0.4
#> 3        d    cal 0.781  0.005    2.8    0.0 0.773  0.005    5.0    0.0
#> 4        d    val 0.742  0.017    1.0    1.0 0.733  0.018    3.0    0.0
#> 7       RE    cal 0.481  0.016    3.7    0.0 0.481  0.016    3.8    0.0
#> 8       RE    val 0.374  0.038    1.7    0.4 0.375  0.039    1.7    0.6
#> 1       CE    cal 0.458  0.006    3.8    0.0 0.459  0.006    3.6    0.0
#> 2       CE    val 0.337  0.045    1.7    0.4 0.340  0.045    1.7    0.6
#>     MT_M MT_SD MT_AR MT_S1 BMT_M BMT_SD BMT_AR BMT_S1  RF_M RF_SD RF_AR
#> 5  0.676 0.007   3.9   0.0 0.686  0.006    2.0      0 0.959 0.002     1
#> 6  0.638 0.030   1.9   0.5 0.638  0.024    3.8      0 0.543 0.044     5
#> 9  0.389 0.002   3.5   0.0 0.384  0.003    2.0      0 0.173 0.004     1
#> 10 0.421 0.014   2.0   0.4 0.421  0.009    3.8      0 0.462 0.019     5
#> 11 0.734 0.005   3.5   0.0 0.725  0.005    2.0      0 0.326 0.009     1
#> 12 0.813 0.027   2.0   0.4 0.813  0.023    3.8      0 0.891 0.029     5
#> 3  0.781 0.005   2.7   0.0 0.782  0.006    2.7      0 0.966 0.002     1
#> 4  0.737 0.023   1.6   0.8 0.730  0.019    3.5      0 0.682 0.031     5
#> 7  0.482 0.014   3.7   0.0 0.494  0.013    2.0      0 0.898 0.005     1
#> 8  0.364 0.047   2.2   0.3 0.365  0.038    3.6      0 0.239 0.060     5
#> 1  0.459 0.007   3.8   0.0 0.472  0.008    2.0      0 0.894 0.006     1
#> 2  0.328 0.049   2.2   0.3 0.329  0.043    3.6      0 0.196 0.054     5
#>    RF_S1
#> 5      1
#> 6      0
#> 9      1
#> 10     0
#> 11     1
#> 12     0
#> 3      1
#> 4      0
#> 7      1
#> 8      0
#> 1      1
#> 2      0
experiment_1[[2]] # See a ggplot of mean bias for validation data
```

![](README-example%202a-1.png)

``` r
## basic example with tunned parameter values
library(dendroTools)
data(example_dataset_1)
experiment_2 <- compare_methods(formula = MVA~.,
dataset = example_dataset_1, k = 3, repeats = 50, neurons = 1,
MT_M = 4, MT_N = FALSE, MT_U = FALSE, MT_R = FALSE, BMT_P = 100,
BMT_I = 100, BMT_M = 4, BMT_N = FALSE, BMT_U = FALSE, BMT_R = FALSE,
RF_P = 100, RF_I = 100, RF_depth= 0, multiply = 5)
```

Aain, to see the outputs, use the following code.

``` r
experiment_2[[1]] # See a data frame results
#>    Measure Period MLR_M MLR_SD MLR_AR MLR_S1 ANN_M ANN_SD ANN_AR ANN_S1
#> 5        r    cal 0.677  0.006   3.94   0.00 0.678  0.006   3.16   0.00
#> 6        r    val 0.638  0.037   1.48   0.52 0.638  0.036   2.18   0.46
#> 9     RMSE    cal 0.388  0.003   3.64   0.00 0.388  0.003   3.48   0.00
#> 10    RMSE    val 0.422  0.015   1.50   0.50 0.423  0.015   2.08   0.50
#> 11    RSSE    cal 0.733  0.005   3.64   0.00 0.732  0.006   3.58   0.00
#> 12    RSSE    val 0.822  0.049   1.56   0.46 0.822  0.050   2.08   0.50
#> 3        d    cal 0.783  0.006   2.72   0.00 0.776  0.006   4.94   0.00
#> 4        d    val 0.735  0.025   1.06   0.94 0.727  0.025   3.00   0.06
#> 7       RE    cal 0.491  0.030   3.70   0.00 0.492  0.030   3.56   0.00
#> 8       RE    val 0.356  0.080   1.60   0.44 0.357  0.082   2.06   0.52
#> 1       CE    cal 0.461  0.008   3.74   0.00 0.462  0.008   3.48   0.00
#> 2       CE    val 0.309  0.099   1.60   0.44 0.309  0.101   2.06   0.52
#>     MT_M MT_SD MT_AR MT_S1 BMT_M BMT_SD BMT_AR BMT_S1  RF_M RF_SD RF_AR
#> 5  0.676 0.006  4.14  0.00 0.689  0.007   2.00   0.00 0.960 0.003  1.00
#> 6  0.628 0.049  1.98  0.50 0.626  0.042   3.64   0.02 0.532 0.060  4.96
#> 9  0.388 0.003  3.94  0.00 0.383  0.004   2.00   0.00 0.170 0.005  1.00
#> 10 0.427 0.020  2.02  0.46 0.428  0.016   3.62   0.02 0.472 0.023  4.98
#> 11 0.733 0.005  3.94  0.00 0.723  0.006   2.00   0.00 0.322 0.008  1.00
#> 12 0.829 0.053  2.08  0.40 0.831  0.049   3.54   0.04 0.921 0.061  4.98
#> 3  0.782 0.006  2.96  0.00 0.785  0.006   2.62   0.00 0.967 0.002  1.00
#> 4  0.727 0.034  1.74  0.74 0.721  0.028   3.48   0.00 0.668 0.039  4.94
#> 7  0.490 0.030  3.98  0.00 0.505  0.029   2.00   0.00 0.902 0.008  1.00
#> 8  0.345 0.087  2.10  0.40 0.342  0.081   3.50   0.04 0.189 0.114  4.98
#> 1  0.460 0.008  4.02  0.00 0.475  0.009   2.00   0.00 0.896 0.005  1.00
#> 2  0.297 0.103  2.10  0.40 0.294  0.098   3.50   0.04 0.129 0.142  4.98
#>    RF_S1
#> 5      1
#> 6      0
#> 9      1
#> 10     0
#> 11     1
#> 12     0
#> 3      1
#> 4      0
#> 7      1
#> 8      0
#> 1      1
#> 2      0
experiment_2[[2]] # See a ggplot of mean bias for validation data
```

![](README-example%203a-1.png)

Authors
-------

-   **Jernej Jevšenak**

Collaborators
-------------

-   **Prof. dr. Tom Levanič**
