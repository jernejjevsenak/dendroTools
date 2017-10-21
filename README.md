
<!-- README.md is generated from README.Rmd. Please edit that file -->
dendroTools
===========

The core purpose of the dendroTools package is to introduce novel dendroclimatological methods to study linear and nonlinear relationships between climate data and tree-ring sequences. There are two core functions. The first core function is daily\_response(), which finds the optimal sequence of days that are linearly or nonlinearly related to a tree-ring proxy records. The second core function is compare\_methods() which calculates several performance measures for train and test data for different regression methods.

To use daily\_response(), two data frames are required, one with daily climate data, e.g. mean daily temperature; and one with tree-ring proxy records. Example data is provided, so users can see, how data frames should be organized. The daily\_response() calculates all possible values of a selected statistical measure between response variable(s) and daily environmental data. Calculations are based on a moving window, which runs through daily environmental data and calculates moving averages.

compare\_methods() calcualtes performance measures for multiple linear regression (MLR), artificial neural networks with Bayesian regularization training algorithm (ANN), M5P model trees (MT), model trees with bagging (BMT) and random forest of regression trees (RF). Calculated performance measures are correlation coefficient (r), root mean squared error (RMSE), root relative squared error (RSSE), index of agreement (d), reduction of error (RE), coefficient of efficiency (CE) nad mean bias.

Installation
------------

You can install dendroTools using:

``` r
# install.packages("devtools")
devtools::install_github("jernejjevsenak/dendroTools") # current version under development

install.packages("devtools") # from CRAN
```

Examples
--------

This is a basic example which demonstrates how to use the daily\_response():

``` r
## basic example code
library(dendroTools)
data(LJ_daily_temperatures) 
data(example_proxies_1)
result1 <- daily_response(response = example_proxies_1, env_data = daily_temperatures_example, 
                            method = "lm", measure = "r.squared", lower_limit = 30, upper_limit = 270, remove_insignificant = TRUE, row_names_subset = TRUE)
```

This function is computationally intensive and it takes a while to calculate all possible values. Especially, if nonlinear "brnn" method is used. While calculating, progress bar shows the percentage of calculated values.

The return of this function is a list with four elements: @calculations, @method, @measure and @optimized\_result. The return is organized in a way, that can be used by three plotting functions: plot\_extreme(), plot\_specific() and plot\_heatmap(). Function plot\_extreme() graphs a line plot of a row with the highest calculated measure. It indicates the sequence of days, that give the best value of selected statistical measure. With plot\_specific(), measures with selected window width are plotted. Function plot\_heatmap() is a visual representation of calculated values.

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

Here is an example of using the compare\_methods(). The core idea is to compare different regression methods and select the best model for climate reconstruction. Users should spend some time to understand what are the tunning parameters of specific method. To boost the performance of regression methods, set the parameter use\_caret = TRUE, to optimize some parameters automatically.

``` r
library(dendroTools)
data(example_dataset_1)
experiment_1 <- compare_methods(formula = MVA ~ ., dataset = example_dataset_1, k = 3, repeats = 100, 
                                use_caret = TRUE)
#> Warning: package 'caret' was built under R version 3.4.2
#> Loading required package: lattice
#> Loading required package: ggplot2
#> Warning: package 'brnn' was built under R version 3.4.1
#> Loading required package: Formula
#> Warning: package 'Formula' was built under R version 3.4.1
#> randomForest 4.6-12
#> Type rfNews() to see new features/changes/bug fixes.
#> 
#> Attaching package: 'randomForest'
#> The following object is masked from 'package:ggplot2':
#> 
#>     margin
#> Warning: package 'bindrcpp' was built under R version 3.4.1
```

To see the outputs, use the following code.

``` r
experiment_1[[1]] # See a data frame results of mean and standard deviation for different methods
#>    Measure Period Mean MLR Std MLR Mean ANN Std ANN Mean MT Std MT
#> 1        r    cal    0.721   0.042    0.729   0.045   0.721  0.045
#> 2        r    val    0.691   0.096    0.692   0.098   0.677  0.098
#> 3     RMSE    cal    0.401   0.027    0.396   0.028   0.401  0.028
#> 4     RMSE    val    0.439   0.057    0.440   0.059   0.448  0.063
#> 5     RSSE    cal    0.690   0.043    0.682   0.049   0.690  0.047
#> 6     RSSE    val    0.783   0.119    0.785   0.124   0.797  0.120
#> 7        d    cal    0.819   0.034    0.819   0.038   0.818  0.036
#> 8        d    val    0.776   0.056    0.770   0.056   0.765  0.061
#> 9       RE    cal    0.563   0.065    0.573   0.069   0.562  0.069
#> 10      RE    val    0.439   0.162    0.435   0.174   0.418  0.166
#> 11      CE    cal    0.522   0.060    0.533   0.066   0.522  0.064
#> 12      CE    val    0.373   0.203    0.369   0.217   0.351  0.205
#>    Mean BMT Std BMT Mean RF Std RF
#> 1     0.736   0.043   0.931  0.010
#> 2     0.681   0.097   0.615  0.119
#> 3     0.392   0.027   0.237  0.017
#> 4     0.444   0.058   0.474  0.068
#> 5     0.675   0.046   0.409  0.029
#> 6     0.791   0.116   0.843  0.127
#> 7     0.824   0.035   0.943  0.010
#> 8     0.766   0.056   0.714  0.077
#> 9     0.581   0.067   0.846  0.026
#> 10    0.427   0.160   0.350  0.184
#> 11    0.542   0.063   0.832  0.023
#> 12    0.361   0.198   0.274  0.230
```

``` r
experiment_1[[2]] # See a data frame results of average rank and share of rank 1 for different methods
#>    Measure Period Avg rank MLR Share of rank 1 MLR Avg rank ANN
#> 5        r    cal        3.917               0.000        3.057
#> 6        r    val        2.047               0.307        2.150
#> 9     RMSE    cal        3.850               0.000        3.120
#> 10    RMSE    val        2.053               0.367        2.463
#> 11    RSSE    cal        3.857               0.000        3.117
#> 12    RSSE    val        2.057               0.363        2.463
#> 3        d    cal        3.187               0.000        3.877
#> 4        d    val        1.663               0.520        2.720
#> 7       RE    cal        3.857               0.000        3.117
#> 8       RE    val        2.057               0.363        2.463
#> 1       CE    cal        3.857               0.000        3.120
#> 2       CE    val        2.057               0.363        2.463
#>    Share of rank 1 ANN Avg rank MT Share of rank 1 MT Avg rank BMT
#> 5                0.000       4.043              0.000        2.173
#> 6                0.460       2.487              0.273        3.163
#> 9                0.000       3.990              0.000        2.207
#> 10               0.297       2.537              0.313        3.037
#> 11               0.000       3.997              0.000        2.207
#> 12               0.297       2.540              0.310        3.037
#> 3                0.000       3.413              0.000        2.673
#> 4                0.223       2.210              0.417        3.150
#> 7                0.000       3.997              0.000        2.207
#> 8                0.297       2.540              0.310        3.037
#> 1                0.000       3.997              0.000        2.210
#> 2                0.297       2.540              0.310        3.037
#>    Share of rank 1 BMT Avg rank RF Share of rank 1 RF
#> 5                0.000       1.000              1.000
#> 6                0.110       4.313              0.127
#> 9                0.000       1.000              1.000
#> 10               0.173       4.100              0.167
#> 11               0.000       1.000              1.000
#> 12               0.173       4.100              0.167
#> 3                0.000       1.000              1.000
#> 4                0.157       4.443              0.103
#> 7                0.000       1.000              1.000
#> 8                0.173       4.100              0.167
#> 1                0.000       1.000              1.000
#> 2                0.173       4.100              0.167
```

``` r
experiment_1[[3]] # See a ggplot of mean bias for calibration data
```

![](README-visualise%20element%203-1.png)

``` r
experiment_1[[4]] # See a ggplot of mean bias for validation data
```

![](README-visualise%20element%204-1.png)

``` r
experiment_1[[5]] # See a data frame with the parameter values used by regression methods
#>    Method   Parameter Value
#> 1     ANN ANN_neurons     1
#> 2      MT        MT_M     4
#> 3      MT        MT_N FALSE
#> 4      MT        MT_U FALSE
#> 5      MT        MT_R FALSE
#> 6     BMT       BMT_P   100
#> 7     BMT       BMT_I   100
#> 8     BMT       BMT_M     4
#> 9     BMT       BMT_N FALSE
#> 10    BMT       BMT_U FALSE
#> 11    BMT       BMT_R FALSE
#> 12     RF     RF_mtry     2
#> 13     RF RF_maxnodes     4
#> 14     RF    RF_ntree   200
```

Authors
-------

-   **Jernej Jevšenak**

Collaborators
-------------

-   **Prof. dr. Tom Levanič**
