# dendroTools 0.0.7.

* There is a new vignette created for the compare_methods() description
* For the optimization phase, you can now choose between RMSE and RSquared metrics for final model selection
* All tuning parameters have new arguments, which allows users to specify a vector of possible values to be tested in the optimization process. 
* caret package was removed from dependencies, since we implemented our own optimization functions
* Random forest model from the randomForest package was replaced with random forest model from the RWeka package 
* The returns argument from the compare_methods() was removed
* For the compare_methods(), the bias distribution plots were changed from density plots to histograms
* Package version is changed to 0.0.7
