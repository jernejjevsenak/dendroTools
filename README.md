
<!-- README.md is generated from README.Rmd. Please edit that file -->
dendroTools
===========

The core purpose of the dendroTools package is to introduce novel dendroclimatological methods to study linear and nonlinear relationship between daily climate data and tree-ring sequences. The core function is daily\_response(), which finds the optimal sequence of days that are linearly or nonlinearly related to a tree-ring proxy records. . To use daily\_response(), two data frames are required, one with daily climate data, e.g. temperatures; and one with tree-ring proxy records. Example data is provided, so users can see, how data frames should be organized. The daily\_response() calculates all possible values of a selected statistical measure between response variables and daily environmental data. Calculations are based on a moving window, which runs through daily environmental data and calculates moving averages.

Installation
------------

You can install dendroTools from github with:

``` r
# install.packages("devtools")
devtools::install_github("jernejjevsenak/dendroTools")
```

Examples
--------

This is a basic example which shows you how to use the package:

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

Authors
-------

-   **Jernej Jevšenak**

Collaborators
-------------

-   **Prof. dr. Tom Levanič**
