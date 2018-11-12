# dendroTools 1.0.2.

* Package version is changed to 1.0.2
* There have been many issues with RWeka package, which depends on rJava. Therefore, functions from RWeka are replaced with other functions. For random forest model we now use randomForest package
* There are six new output elements, all of them are residual diagnostic plots for calibration, holdout and edge data.
* Due to many problems related to rJava in R, MT and BMT methods from RWeka are now joined and replaced by cubist function from Cubist r package. Cubist fits a model tree or ensemble of model trees, if committees argument is more than 1.
