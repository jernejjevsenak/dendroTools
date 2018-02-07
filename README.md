
<!-- README.md is generated from README.Rmd. Please edit that file -->
dendroTools
===========

The core purpose of the dendroTools package is to introduce novel dendroclimatological methods to study linear and nonlinear relationships between environmental (climate) data and tree-ring sequences. Currently, there are two core functions. The first core function is the daily\_response(), which finds the optimal sequence of days that are linearly or nonlinearly related to a tree-ring proxy records. The second core function is compare\_methods() which calculates several performance metrics for train and test data for different regression methods.

To use daily\_response(), two data frames are required, one with daily climate data, e.g. mean daily temperature; and one with tree-ring proxy records. Example data is provided, so users can see, how data frames should be organized. The daily\_response() calculates all possible values of a selected statistical metric between response variable(s) and daily environmental data. Calculations are based on a moving window, which runs through daily environmental data and calculates moving averages.

compare\_methods() calcualtes performance metrics for multiple linear regression (MLR), artificial neural networks with Bayesian regularization training algorithm (BRNN), M5P model trees (MT), model trees with bagging (BMT) and random forest of regression trees (RF). Calculated performance metrics are correlation coefficient (r), root mean squared error (RMSE), root relative squared error (RSSE), index of agreement (d), reduction of error (RE), coefficient of efficiency (CE) nad mean bias.

Installation
------------

You can install dendroTools using:

``` r
library("devtools")
devtools::install_github("jernejjevsenak/dendroTools") # current version under development

install.packages("dendroTools") # from CRAN
```

Authors
-------

-   **Jernej Jevšenak**

Collaborators
-------------

-   **Prof. dr. Tom Levanič**
