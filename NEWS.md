# dendroTools 0.0.4.

* Fifth element is given to the compare_methods() output list. This is a data frame with the final values of parameters used by different regression methods. 
* New function is added: round_df(). This function rounds all numeric columns in a data frame. 
* Datasets are updated
* the method description as part of the title is added to the plot_extreme(), plot_specific() and plot_heatmap()
* iter function was joined with the compare_methods()
* random forests model was previously based on RWeka package. From 0.0.4 version on, randomForest package is used to fit random forest models
* Package version is changed to 0.0.4
* For the compare_methods(), new parameter is introduced: use_caret (boolean). The parameter is used for model tunning by caret package. 
* Progress bar is given to compare_methods() and daily_response()

